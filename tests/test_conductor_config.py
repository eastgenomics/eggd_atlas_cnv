import json
from pathlib import Path

PROVIDED = json.loads(Path("conductor/atlas_cnv_pon_provided.example.json").read_text())
BUILT = json.loads(Path("conductor/atlas_cnv_pon_built.example.json").read_text())


def _by_name(cfg, name):
    return next(v for v in cfg["executables"].values() if v.get("name") == name)


def test_provided_atlas_uses_workflow_input_cn_reference():
    atlas = _by_name(PROVIDED, "eggd_atlas_cnv")
    assert atlas["per_sample"] is True
    assert "cnvkit_cn_reference" in atlas["inputs"]          # workflow-level input name
    assert "stage-cnvkit_batch.cn_reference" not in atlas["inputs"]  # never stage-qualified
    assert "sample_id" in atlas["inputs"]                    # required input present
    assert "input_bai" in atlas["inputs"]                    # separately supplied chr-prefixed BAI


def test_built_pon_step_is_per_run_and_holds():
    pon = _by_name(BUILT, "eggd_cgp-cnvkit-pon")
    assert pon["per_sample"] is False
    assert pon.get("hold") is True          # must hold to gather the coverage array
    cov = _by_name(BUILT, "eggd_cgp-cnvkit-coverage")
    assert cov["analysis"] in pon["depends_on"]


def test_built_coverage_and_chr_prefix_per_sample_and_required_inputs():
    cov = _by_name(BUILT, "eggd_cgp-cnvkit-coverage")
    chrp = _by_name(BUILT, "eggd_chr_prefix")
    assert cov["per_sample"] is True and chrp["per_sample"] is True
    assert "sample_id" in cov["inputs"]                      # coverage needs sample_id
    assert chrp["analysis"] in json.dumps(cov["inputs"])     # coverage consumes prefixed BAM


def test_built_inputs_filter_is_include_regex():
    pon = _by_name(BUILT, "eggd_cgp-cnvkit-pon")
    pats = pon["inputs_filter"]["coverage_files"]
    # include filter: patterns match REAL specimen IDs, not the controls we want to drop
    assert all("NA12878" not in p for p in pats)


def test_built_atlas_links_cn_reference_from_pon_analysis():
    atlas = _by_name(BUILT, "eggd_atlas_cnv")
    pon = _by_name(BUILT, "eggd_cgp-cnvkit-pon")
    assert atlas["per_sample"] is True and pon["analysis"] in atlas["depends_on"]
    assert "sample_id" in atlas["inputs"]
    assert atlas["inputs"]["input_bam"]["$dnanexus_link"]["analysis"] == "analysis_1"
    assert atlas["inputs"]["input_bai"]["$dnanexus_link"]["analysis"] == "analysis_1"
    found = json.dumps(atlas["inputs"])
    assert f'"{pon["analysis"]}"' in found and '"cn_reference"' in found
