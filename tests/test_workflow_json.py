import json
from pathlib import Path

WF = json.loads(Path("dxworkflow.json").read_text())
STAGE_IDS = {s["id"] for s in WF["stages"]}


def test_nine_stages():
    assert STAGE_IDS == {"chr_prefix", "amber", "cobalt", "sage", "purple", "qc_flags",
                         "cnvkit_batch", "cnv_chr_strip", "purple_plotter"}
    assert len(WF["stages"]) == 9


def _links(stage_input):
    for v in stage_input.values():
        link = v.get("$dnanexus_link") if isinstance(v, dict) else None
        if isinstance(link, dict) and "stage" in link:
            yield link["stage"], link["outputField"]


def test_all_stage_links_reference_existing_stages():
    for s in WF["stages"]:
        for src, field in _links(s.get("input", {})):
            assert src in STAGE_IDS, f"{s['id']} links unknown stage {src}"


def test_downstream_bams_come_from_chr_prefix():
    # amber/cobalt/sage/cnvkit_batch must take tumour_bam/bai from the chr_prefix stage
    for sid in ("amber", "cobalt", "sage", "cnvkit_batch"):
        st = next(s for s in WF["stages"] if s["id"] == sid)
        assert st["input"]["tumour_bam"]["$dnanexus_link"] == {"stage": "chr_prefix", "outputField": "output_bam"}, sid
        assert st["input"]["tumour_bai"]["$dnanexus_link"] == {"stage": "chr_prefix", "outputField": "output_bai"}, sid


def test_purity_ploidy_scalars_linked_into_batch():
    batch = next(s for s in WF["stages"] if s["id"] == "cnvkit_batch")
    fields = {f for src, f in _links(batch["input"]) if src == "purple"}
    assert {"purity", "ploidy"} <= fields


def test_cn_reference_is_a_workflow_input_not_a_stage_link():
    batch = next(s for s in WF["stages"] if s["id"] == "cnvkit_batch")
    cn = batch["input"]["cn_reference"]["$dnanexus_link"]
    assert "workflowInputField" in cn   # supplied as a workflow input, never a stage output


def test_plotter_consumes_chr_prefixed_files_not_stripped():
    # igv.js/hg38 needs chr names: plotter must read cnvkit_batch/purple, never cnv_chr_strip
    plot = next(s for s in WF["stages"] if s["id"] == "purple_plotter")
    srcs = {src for src, _ in _links(plot["input"])}
    assert {"purple", "qc_flags", "cnvkit_batch"} <= srcs
    assert "cnv_chr_strip" not in srcs


def test_strip_stage_consumes_cnvkit_and_purple():
    strip = next(s for s in WF["stages"] if s["id"] == "cnv_chr_strip")
    srcs = {src for src, _ in _links(strip["input"])}
    assert {"cnvkit_batch", "purple"} <= srcs


def test_declared_workflow_outputs():
    outs = {o["name"]: o["outputSource"]["$dnanexus_link"] for o in WF["outputs"]}
    # CNV call files are promoted in their chr-stripped form (from cnv_chr_strip)
    for name in ("igv_html", "purple_cnv_somatic_nochr", "purple_cnv_gene_nochr",
                 "cnvkit_call_cns_nochr", "cnvkit_genemetrics_nochr", "cnvkit_cnr_nochr",
                 "qc_report", "purple_tar", "purity", "ploidy"):
        assert name in outs, f"missing workflow output {name}"
    assert outs["purple_cnv_somatic_nochr"]["stage"] == "cnv_chr_strip"
    assert outs["cnvkit_call_cns_nochr"]["stage"] == "cnv_chr_strip"
    assert outs["purple_tar"] == {"stage": "purple", "outputField": "purple_tar"}


def test_workflow_output_sources_are_real_stage_outputs():
    for o in WF["outputs"]:
        link = o["outputSource"]["$dnanexus_link"]
        assert link["stage"] in STAGE_IDS


# Map each stage id to the app dir whose dxapp.json defines its outputSpec.
STAGE_APP = {
    "chr_prefix": "eggd_chr_prefix", "amber": "eggd_cgp-amber", "cobalt": "eggd_cgp-cobalt",
    "sage": "eggd_cgp-sage", "purple": "eggd_cgp-purple", "qc_flags": "eggd_cgp-qc-flags",
    "cnv_chr_strip": "eggd_cnv_chr_strip",
    "cnvkit_batch": "eggd_cgp-cnvkit-batch", "purple_plotter": "eggd_purple_plotter",
}


def _outputspec(stage_id):
    app = json.loads(Path(f"apps/{STAGE_APP[stage_id]}/dxapp.json").read_text())
    return {o["name"]: o["class"] for o in app["outputSpec"]}


def _inputspec(stage_id):
    app = json.loads(Path(f"apps/{STAGE_APP[stage_id]}/dxapp.json").read_text())
    return {i["name"]: i["class"] for i in app["inputSpec"]}


def test_stage_input_names_exist_in_app_inputspec():
    # every input key on a stage must be a real inputSpec field of that stage's app,
    # and a stage->stage file/scalar link must match classes
    for s in WF["stages"]:
        ispec = _inputspec(s["id"])
        for name, v in s.get("input", {}).items():
            base = name.split(".")[-1]          # tolerate stage-qualified keys
            assert base in ispec, f"{s['id']} sets unknown input {name}"
            link = v.get("$dnanexus_link") if isinstance(v, dict) else None
            if isinstance(link, dict) and "stage" in link:
                assert _outputspec(link["stage"])[link["outputField"]] == ispec[base], \
                    f"class mismatch: {s['id']}.{base} <- {link['stage']}.{link['outputField']}"


def test_linked_output_fields_exist_in_app_outputspec():
    # every stage->stage link and every workflow output must name a real app outputSpec field
    for s in WF["stages"]:
        for src, field in _links(s.get("input", {})):
            assert field in _outputspec(src), f"{s['id']} links missing {src}.{field}"
    for o in WF["outputs"]:
        link = o["outputSource"]["$dnanexus_link"]
        spec = _outputspec(link["stage"])
        assert link["outputField"] in spec, f"workflow output {o['name']} links missing field"
        assert spec[link["outputField"]] == o["class"], f"class mismatch for {o['name']}"


def test_purple_promoted_cnv_outputs_exist():
    # guards the hard M6->M9 dependency: PURPLE must declare the CNV TSV outputs
    spec = _outputspec("purple")
    assert {"cnv_somatic_tsv", "cnv_gene_tsv"} <= set(spec)
