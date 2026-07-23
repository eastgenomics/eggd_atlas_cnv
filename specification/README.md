# eggd_atlas_cnv — per-sample somatic CNV workflow (PURPLE + CNVkit + IGV)

`eggd_atlas_cnv` is a DNAnexus **workflow** that takes a chr-prefixed tumour BAM and its
BAI index and produces a combined somatic copy-number call set: it produces a genome-wide PURPLE profile (purity/ploidy-fitted, with an optional ploidy cap), a focal
gene-level CNVkit profile (using either a supplied panel-of-normals or one built once per
run by eggd_conductor outside the workflow), and a single IGV.js HTML report that overlays
both (igv.js loads from a CDN at view time — not fully offline).

It formalises the two previously separate tracks — the PURPLE backbone workflow
(`cnv-backbone-purple-atlas`) and the CNVkit gene atlas (`cnv-genes-atlas`) — into
one locked, versioned, per-sample workflow composed entirely of **DNAnexus apps**
(`org-emee_1`, `aws:eu-central-1`).

## What this workflow does

Given a chr-prefixed tumour BAM and matching BAI (produced by a separately-run
`eggd_chr_prefix`), a sample ID, and PoN configuration, the workflow runs (AMBER+COBALT
are parallel stages, so there are seven DNAnexus app stages: `amber`, `cobalt`, `sage`,
`purple`, `cnvkit_batch`, `cnv_chr_strip`, `purple_plotter` — see DESIGN §3):

1. Runs **AMBER** (BAF per germline site) and **COBALT** (read-depth ratios) in parallel.
2. Runs **SAGE** to produce a somatic SNV/indel VCF (panel mode).
3. Runs **PURPLE** to fit purity, ploidy, and genome-wide copy-number segments —
   applying a **maximum-ploidy cap** when requested, either as a hard cap
   (`max_ploidy`) or conditionally on a purity threshold (e.g. cap ploidy at 2 when
   fitted purity < 0.35). PURPLE emits purity and ploidy as scalar job outputs.
4. Runs **CNVkit batch** (coverage → fix → segment → call → plot → genemetrics)
   against a supplied panel-of-normals (`cnvkit_cn_reference`), feeding PURPLE's purity
   and ploidy into integer CN calling.
5. Runs **eggd_purple_plotter** to overlay PURPLE + AMBER + CNVkit into one
   IGV.js HTML report (igv.js loads from a CDN at view time — not fully offline).
6. Runs **eggd_cnv_chr_strip** to write Ensembl-named (`*.nochr.*`) copies of the CNV call
   files for downstream apps, **retaining the chr-prefixed originals**.

The panel-of-normals is **an input to the workflow, never built inside it**. When a run
needs a fresh PoN, it is built **once per sequencing run** by `eggd_conductor` as a
per-run step and linked into every per-sample workflow job (see "Orchestration" below) —
not rebuilt 48 times.

## Status of this document set

These documents are the **complete design and build specification**. A fresh agent
session (or a human developer) should be able to open this directory, read the files
in order, and build a working, tested workflow without needing the original
conversation.

Read in this order:

1. **README.md** (this file) — orientation and quick start
2. **DESIGN.md** — architecture, the app-conversion strategy, the two-pass ploidy-cap
   state machine, scalar output linking (JBOR), the eggd_conductor per-run PoN wiring, error and test strategy
3. **IMPLEMENTATION.md** — milestone-by-milestone TDD build plan with code sketches
4. **REFERENCE.md** — app I/O catalogue, DNAnexus resource file IDs, workflow-link
   schema + declared outputs, conductor config, env vars, glossary
5. **build-prompt.md** — the ready-to-use build prompt: build rules, per-milestone
   verification table, and the invariants a fresh agent must uphold

## Project layout (target)

```
eggd_atlas_cnv/
├── dxworkflow.json                 ← the locked workflow (built with dx build)
├── atlas_helpers/                  ← pure-Python decision logic (unit-tested here)
│   ├── __init__.py
│   ├── ploidy_gate.py              ← two-pass ploidy-cap decision (imported by purple app)
│   ├── purity.py                   ← parse *.purple.purity.tsv → purity/ploidy/sex
│   └── chr_strip.py                ← strip chr from CNV file chromosome column (imported by cnv_chr_strip)
├── conductor/
│   └── atlas_cnv_pon_provided.example.json / atlas_cnv_pon_built.example.json  ← eggd_conductor executables (one per PoN topology)
├── scripts/
│   ├── resource_ids.env.template   ← all DNAnexus file/app IDs to fill in
│   ├── app_ids.json                ← stage-app-name → built DNAnexus app ID (single source
│   │                                  of truth for app identity; each app builds/publishes
│   │                                  itself from its own per-app GitHub repo, not from here)
│   └── run_e2e.sh                  ← run the workflow on one sample end-to-end
├── tests/
│   ├── test_ploidy_gate.py         ← ploidy-cap decision unit tests
│   ├── test_purity.py              ← purity.tsv parser unit tests
│   ├── test_chr_strip.py           ← chr-strip transform unit tests
│   ├── test_workflow_json.py       ← dxworkflow.json structure + link validity
│   └── test_conductor_config.py    ← conductor executables block structure + link validity
├── pyproject.toml
└── specification/                  ← this directory
```

## Quick start (once built)

```bash
# 1. Init and install dev deps for the pure-Python helpers/tests
python3 -m venv .venv && .venv/bin/pip install -e ".[dev]"

# 2. Log in to DNAnexus (test_workflow_json.py validates dxworkflow.json against the
#    live, deployed apps named in scripts/app_ids.json — no local app source needed)
dx login

# 3. Run the unit + workflow-JSON tests
.venv/bin/pytest tests/ -v

# 4. Each app builds and publishes itself from its own per-app GitHub repo (see
#    scripts/app_ids.json for the eastgenomics/<app-name> repo per stage). Fill in
#    remaining reference file IDs, then select the build project for the workflow itself
cp scripts/resource_ids.env.template scripts/resource_ids.env   # edit IDs
dx select project-XXXX

# 5. Regenerate + build the workflow (references apps by exact app ID, not by name)
python3 scripts/gen_workflow.py
dx build --workflow . --overwrite

# 6. Run the workflow on one sample (chr-prefixed BAM/BAI supplied by a separately-run
#    eggd_chr_prefix; PoN supplied as a workflow input)
dx run workflow-XXXX \
  -iinput_bam="project-XXXX:file-CHR-PREFIXED-BAM"   \
  -iinput_bai="project-XXXX:file-CHR-PREFIXED-BAI"    \
  -isample_id="25330S0047"              \
  -iploidy_cap_purity_threshold=0.35    \
  -icnvkit_cn_reference="project-XXXX:file-PON" \
  --destination "project-XXXX:/atlas_cnv/25330S0047/" --priority high --watch
```

## Orchestration (eggd_conductor)

A whole sequencing run is driven by [`eggd_conductor`](https://github.com/eastgenomics/eggd_conductor),
which reads an assay config `executables` block and launches each step **per run**
(`"per_sample": false`) or **per sample** (`"per_sample": true`), wiring outputs to
downstream inputs with `$dnanexus_link {analysis, stage, field}` and `depends_on`.

The PoN choice is a **config decision, not a workflow branch**:

| Use case | Conductor topology |
|---|---|
| **PoN provided** | Per-sample `eggd_atlas_cnv` only; `cnvkit_cn_reference` is a static `$dnanexus_link {project,id}`. |
| **PoN built** | `eggd_chr_prefix` (per-sample) → `eggd_cgp-cnvkit-coverage` (per-sample, parallel) → `eggd_cgp-cnvkit-pon` (per-run, `depends_on` + `hold`, gathers the coverage array with an **include** `inputs_filter` matching real specimen IDs) → `eggd_atlas_cnv` (per-sample, `cnvkit_cn_reference` ← the per-run `cn_reference`; `input_bam`/`input_bai` ← the same per-sample `chr_prefix` analysis). |

See `conductor/atlas_cnv_pon_{provided,built}.example.json` and REFERENCE §4.

## Target user

A clinical-bioinformatics scientist at CUH / East Genomics running the CGP + backbone
somatic panel. They have DNAnexus access to `org-emee_1` in `aws:eu-central-1`, need a
single reproducible per-sample call that combines genome-wide (PURPLE) and focal
(CNVkit) copy number, and want one HTML deliverable to inspect the result. They should
not need to know the internal tool flags or the two-track history.

## Non-goals

- **Not a batch/cohort runner.** One BAM per invocation. Cohort runs are driven by
  launching the workflow N times (a separate orchestration script), not by scatter.
- **No matched normal.** Tumour-only throughout (AMBER germline-PON sites).
- **No SIGS / CUPpa / CHORD.** Signature, tissue-of-origin and HRD stages from the
  legacy workflow are out of scope for v0.1 (QC-flags treats them as optional inputs,
  left unlinked). They may be re-added as optional stages later.
- **No LOH / cnLOH calling.** Use FACETS separately for copy-neutral LOH.
- **Does not build the PoN.** The panel-of-normals is always a workflow *input*.
  Building it (when needed) is a per-run `eggd_conductor` step, so it runs once per run,
  not once per sample. The workflow contains no PoN-build or reference-resolver stage.
- **No conditional stage skipping in the workflow graph.** DNAnexus native workflows
  cannot branch; the only in-app conditional behaviour is PURPLE's ploidy cap.
- **Not a general HMF-tools wrapper.** Flags are locked to the CGP+backbone panel.
