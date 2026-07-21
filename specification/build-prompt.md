# build-prompt.md — agent build instructions for eggd_atlas_cnv

Copy everything below the line as the opening message to a fresh agent session.

---

You are building **eggd_atlas_cnv**, a per-sample DNAnexus workflow that runs a somatic
CNV pipeline (chr-prefixed BAM/BAI supplied from a separately-run `eggd_chr_prefix`, then
AMBER + COBALT + SAGE → PURPLE → QC-flags, plus CNVkit batch, joined by an IGV.js plotter),
composed of DNAnexus **apps**.

The complete specification lives in `specification/`. **Read all five documents in full,
in order, before writing any code** (this file, `build-prompt.md`, is the fifth — it holds
the invariants you must uphold):

1. `specification/README.md` — orientation, layout, quick start, non-goals, orchestration
2. `specification/DESIGN.md` — architecture, the applet→app conversion, the two-pass
   ploidy-cap state machine (§4.2), scalar output linking, the eggd_conductor per-run PoN
   orchestration (§13), error table, test strategy, limitations
3. `specification/IMPLEMENTATION.md` — the milestone-by-milestone TDD plan (M1–M11) with
   runnable test and code sketches
4. `specification/REFERENCE.md` — app I/O catalogue, DNAnexus resource file IDs, the
   workflow link schema + declared `outputs`, conductor config JSON, glossary
5. `specification/build-prompt.md` (this file) — build rules, per-milestone verification,
   and the invariants below

## Build rules

- **TDD for everything local.** For M1, M2, M9 and M10 write the tests first (they are
  provided in IMPLEMENTATION), run them, confirm they fail (red), then implement until
  green. Do not write implementation before its test.
- **Commit at each milestone** with a concise message (e.g. `M2: ploidy-gate, purity and
  chr_strip helpers + tests`). **No AI attribution** in any commit, file, or comment.
- **DNAnexus milestones (M3–M8, M11) are hand-offs.** They require `dx login`, a build
  project, and real reference data. Prepare the app sources and `dxapp.json`/`code.sh`
  edits fully, run `python3 -c "import json; json.load(...)"` validation locally, but
  STOP before `dx build --app` / `dx run` and hand to the human operator with the exact
  command to run. Do not fabricate app IDs or run outputs. **M9 is local-TDD** (write and
  pass `test_workflow_json.py` first); only its final `dx build --workflow` is a hand-off.
- **Never modify the legacy applet sources in place.** Copy them from
  `/home/wook/Documents/cnv-backbone-purple-atlas/applets/` into `apps/`, then edit the copies.

## Per-milestone verification

| M | Command | Green when |
|---|---|---|
| M1 | `.venv/bin/pytest tests/test_import.py -v` | package imports |
| M2 | `.venv/bin/pytest tests/test_ploidy_gate.py tests/test_purity.py tests/test_chr_strip.py -v` | all green |
| M3–M5 ⚠️ | `dx build --app apps/eggd_cgp-<tool> --overwrite` + smoke | builds; documented output appears |
| M6 ⚠️ | `pytest tests/test_ploidy_gate.py tests/test_purity.py` + `dx build --app apps/eggd_cgp-purple` | builds; conditional re-run + scalar purity/ploidy/sample_sex + cnv_somatic_tsv/cnv_gene_tsv |
| M7 ⚠️ | `dx build --app apps/eggd_cgp-qc-flags --overwrite` | builds; 16-col qc_report |
| M8 ⚠️ | `pytest tests/test_chr_strip.py` + `dx build --app` cnv_chr_strip + coverage/pon/batch/plotter | strip writes `*.nochr.*` keeping originals; five apps build; plotter accepts `.genemetrics.tsv` |
| M9 | `.venv/bin/pytest tests/test_workflow_json.py -v` | 8 stages; chr-prefixed BAM/BAI are workflow inputs; CNV `*_nochr` outputs from cnv_chr_strip; scalar links; cn_reference is a workflow input |
| M10 | `.venv/bin/pytest tests/test_conductor_config.py -v` | per-run PoN topology valid; cn_reference linked from the pon analysis |
| M11 ⚠️ | `bash scripts/run_e2e.sh <sample>` | non-empty `igv_html` |

## Invariants (must hold throughout the entire build)

### The pure helpers stay free of dxpy, network, and unowned I/O
`atlas_helpers/{ploidy_gate,purity}.py` must be importable and testable with no DNAnexus
connection. `ploidy_gate` has no filesystem/network/DNAnexus I/O at all; `purity` only
*reads* the supplied local TSV (no writes, no network). **Verify:**
`grep -rnE "import dxpy|from dxpy|requests|urllib|dx-jobutil" atlas_helpers/` → no output;
and `grep -nE "open\(|write_text|write_bytes" atlas_helpers/ploidy_gate.py` → no output
(`ploidy_gate` performs no file I/O; only `purity`/`chr_strip` read local files, `chr_strip`
also writes its one declared output).

### max_ploidy and ploidy_cap_purity_threshold are mutually exclusive, and the threshold is a fraction
Supplying both is a configuration error, not a silent precedence; and
`ploidy_cap_purity_threshold` must be in `(0, 1]`. **Verify:**
`pytest tests/test_ploidy_gate.py::test_mutually_exclusive tests/test_ploidy_gate.py::test_invalid_threshold_range`
both pass (each raises `PloidyConfigError`).

### The conditional ploidy cap decides on PURPLE's own fitted purity, and re-runs at most once
In CONDITIONAL mode: run PURPLE unbounded (pass 1), read `purity` from
`{s}.purple.purity.tsv`, and only if `purity < threshold` re-run PURPLE **once** with
`-max_ploidy ploidy_cap_value`. Never re-run AMBER/COBALT; never loop more than two PURPLE
passes. **Verify:** M6 smoke log shows exactly one or two `java -jar purple.jar`
invocations, never three; `needs_rerun` boundary tests pass (`0.34→True`, `0.35→False`).

### PURPLE emits purity (float), ploidy (int), sample_sex (string) as typed job outputs
These are the ONLY way purity/ploidy reach CNVkit batch — the workflow cannot parse a
float out of a file (see DESIGN §2.4 for the scalar-output / JBOR mechanism). `ploidy` is
`round(raw)` with a floor of 1 (downstream input class is `int`). **Verify:**
`dx describe <purple-job> --json | jq '.output | {purity, ploidy, sample_sex}'` returns a
float, an int, and a string; and `tests/test_workflow_json.py::test_purity_ploidy_scalars_linked_into_batch` passes.

### The PoN is a workflow input, never built inside the per-sample workflow
The workflow contains no PoN-build, coverage-pooling, or reference-resolver stage.
`cnvkit_batch.cn_reference` is linked from the required workflow input `cnvkit_cn_reference`
(`workflowInputField`), never from a stage output. Building a PoN (when needed) is a
per-run `eggd_conductor` step (`per_sample: false`) that runs once per run and reuses the
existing `eggd_cgp-cnvkit-coverage` + `eggd_cgp-cnvkit-pon` apps. **Verify:**
`tests/test_workflow_json.py::test_cn_reference_is_a_workflow_input_not_a_stage_link`,
`tests/test_conductor_config.py::test_built_pon_step_is_per_run_and_holds`, and
`::test_provided_atlas_uses_workflow_input_cn_reference` pass.

### Every app declares a timeoutPolicy and installs no dependencies at run time
No inline `apt-get install` in any `code.sh`; deps come from `execDepends` or the cnvkit
Docker asset. **Verify:** `grep -RniE '\b(apt|apt-get)\b.*\binstall\b' apps/*/src apps/*/resources` → no output (broader than just `apt-get install`);
`grep -rL timeoutPolicy apps/*/dxapp.json` → no output.

### Tool flags and output names are frozen for converted apps
AMBER/COBALT/SAGE/PURPLE/QC-flags keep every legacy tool flag, reference input, and output
name — downstream links in `dxworkflow.json` and legacy validation depend on them. The
only additions are PURPLE's new ploidy inputs and outputs. **Verify (machine-checkable):**
(a) `tests/test_app_specs.py` asserts each converted app's `dxapp.json` `outputSpec` names
equal the recorded legacy set (plus PURPLE's `purity`/`ploidy`/`sample_sex`/`cnv_somatic_tsv`/
`cnv_gene_tsv`); (b) `grep -c "apt-get\|apt install" apps/<app>/src/code.sh` — removed lines
only; (c) `diff -u <legacy>/src/code.sh apps/<app>/src/code.sh` as supporting evidence that
no tool flag changed. A failing (a) or unexpected flag in (c) blocks the milestone.

### Stage 0 chr_prefix emits exactly one BAM and never starves downstream
`eggd_chr_prefix` is adapted to single-file I/O (`input_bam` → `output_bam`, `output_bai`)
and must **always emit a BAM** — reheadered when needed, or the input passed through when
already chr-prefixed. It must never emit nothing (the public app's skip behaviour would
break the workflow). **Verify:** smoke on an already-chr BAM still produces `output_bam`;
`tests/test_workflow_json.py::test_downstream_bams_come_from_chr_prefix` passes (amber,
cobalt, sage, cnvkit_batch all take their BAM from the `chr_prefix` stage).

### chr-stripping produces copies and never mutates the originals
`eggd_cnv_chr_strip` writes `*.nochr.*` copies with the `chr` prefix removed from the
**chromosome column only** (`chrM`→`MT`); it must not touch gene names/other columns and
must not overwrite its inputs. CNVkit originals stay in `cnvkit_batch`'s outputs, PURPLE
originals in `purple_tar`. The plotter keeps consuming the **chr-prefixed** files (igv.js/
hg38), never the stripped ones. **Verify:** `tests/test_chr_strip.py` passes
(`test_only_chromosome_column_changes`, `chrM→MT`, header-only passthrough); and
`tests/test_workflow_json.py::test_plotter_consumes_chr_prefixed_files_not_stripped` passes.

### Exactly one coverage computation per sample in the workflow
The workflow has no `cnvkit_coverage` stage — `cnvkit_batch` computes its own on-target
coverage. Do not add a second coverage stage. (The per-run PoN build computes coverage
separately via conductor; that is a different, run-level job.) **Verify:** `dxworkflow.json`
has exactly 7 stages and no stage named `cnvkit_coverage`.

## Starting instruction

First, print the target file tree (README §"Project layout") and the contents of
`pyproject.toml` and the two helper test files (`test_ploidy_gate.py`, `test_purity.py`)
from IMPLEMENTATION §3–§4. Then begin at
M1 and proceed milestone by milestone, running each verification command before moving on.
When you reach a ⚠️ DNAnexus milestone, prepare all files, validate JSON locally, and hand
the exact `dx build`/`dx run` command to the operator rather than executing it yourself.
