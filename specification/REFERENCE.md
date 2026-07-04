# REFERENCE — eggd_atlas_cnv

## 1. App I/O catalogue

Each row: input/output name → class (patterns / notes). ⚠️ = new or changed vs legacy.

### 1.0 `eggd_chr_prefix` (folded in, adapted) — stage 0

| Dir | Name | Class | Notes |
|---|---|---|---|
| in | input_bam | file | `*.bam`; raw tumour BAM (Ensembl or already-chr) |
| in | mode | string | fixed `add_chr` (default; not a workflow input) |
| out ⚠️ | output_bam | file | `*.bam`; chr-prefixed (or passthrough) |
| out ⚠️ | output_bai | file | `*.bai`; generated index |

Adapted from public `eggd_chr_prefix` (array I/O → **single-file**; emits input unchanged
when already prefixed). Ubuntu 20.04, samtools asset, `mem1_ssd1_v2_x2`, timeout 2h.

### 1.1 `eggd_cgp-amber` (converted, v1.0.0)

| Dir | Name | Class | Notes |
|---|---|---|---|
| in | tumour_bam | file | chr-prefix GRCh38 BAM |
| in | tumour_bai | file | index |
| in | sample_id | string | output stem |
| in | amber_jar | file | AMBER 4.3-beta.4 JAR |
| in | germline_sites | file | het PON loci `.tsv.gz` |
| out | amber_tar | file | `{s}.amber.tar.gz` |

Instance `mem1_ssd1_v2_x4`; timeout 6h; execDepends: openjdk-21-jre-headless, samtools, tabix.

### 1.2 `eggd_cgp-cobalt` (converted, v1.0.0)

| Dir | Name | Class | Notes |
|---|---|---|---|
| in | tumour_bam, tumour_bai, sample_id | file/file/string | |
| in | cobalt_jar | file | COBALT 3.0-beta.5 |
| in | norm_file | file | target-region normalisation TSV |
| in | diploid_regions | file | `DiploidRegions.38.bed.gz` |
| in | gc_profile | file | `GC_profile.1000bp.38.cnp` |
| in | ref_fasta, ref_fai | file | bgzf FASTA + fai |
| out | cobalt_tar | file | `{s}.cobalt.tar.gz` |

Instance `mem2_ssd1_v2_x4`; timeout 6h. Always `-ref_genome_version 38`.

### 1.3 `eggd_cgp-sage` (converted, v1.0.0)

| Dir | Name | Class | Notes |
|---|---|---|---|
| in | tumour_bam, tumour_bai, sample_id | | |
| in | sage_jar | file | SAGE 5.0-beta.11 |
| in | ref_fasta, ref_fai | file | plain `.fa` + fai |
| in | hotspots_vcf, hotspots_tbi | file | |
| in | panel_bed, hc_bed | file | panel + high-confidence BED |
| in | pon_file | file | SAGE PON `.tsv.gz` |
| in | ensembl_data | file | `ensembl_data.tar.gz` |
| out | somatic_vcf, somatic_vcf_tbi | file | `{s}.sage.somatic.vcf.gz(.tbi)` |

Instance `mem2_ssd1_v2_x16`; timeout 8h. Flags: `-panel_only -high_depth_mode
-skip_msi_jitter -skip_bqr -ref_sample_count 0`.

### 1.4 `eggd_cgp-purple` (converted + extended, v1.0.0) ⚠️

| Dir | Name | Class | Notes |
|---|---|---|---|
| in | sample_id | string | |
| in | purple_jar | file | PURPLE 4.4-beta.12 |
| in | amber_tar, cobalt_tar | file | from AMBER/COBALT |
| in | somatic_vcf, somatic_vcf_tbi | file | **optional**, from SAGE |
| in | gc_profile | file | `GC_profile.1000bp.38.cnp` |
| in | target_regions_bed | file | panel target regions |
| in | ref_fasta, ref_fai | file | bgzf FASTA + fai |
| in | ensembl_data | file | ensembl tarball |
| in ⚠️ | max_ploidy | int | optional; hard `-max_ploidy` cap |
| in ⚠️ | ploidy_cap_purity_threshold | float | optional; purity-conditional cap trigger |
| in ⚠️ | ploidy_cap_value | int | optional, default 2; cap applied on trigger |
| out | purple_tar | file | `{s}.purple.tar.gz` (full output dir) |
| out | purity_tsv | file | `{s}.purple.purity.tsv` |
| out | purity_range_tsv | file | `{s}.purple.purity.range.tsv` |
| out | plots_tar | file | optional; JFreeChart PNGs |
| out ⚠️ | cnv_somatic_tsv | file | `{s}.purple.cnv.somatic.tsv` (segment-level CN calls) |
| out ⚠️ | cnv_gene_tsv | file | `{s}.purple.cnv.gene.tsv` (gene-level CN calls) |
| out ⚠️ | purity | float | fitted purity (final pass) |
| out ⚠️ | ploidy | int | `round(ploidy)`, min 1 (final pass) |
| out ⚠️ | sample_sex | string | `gender` lowercased (male/female/unknown) |

Instance `mem2_ssd1_v2_x4`; timeout 6h. execDepends include R packages + circos (see
DESIGN §4.2 / legacy applet). `max_ploidy` and `ploidy_cap_purity_threshold` are mutually
exclusive (helper raises `PloidyConfigError`).

`*.purple.purity.tsv` columns consumed: `purity`, `ploidy`, `status`, `gender`
(also `wholeGenomeDuplication`, `diploidProportion`, `minPurity`, `maxPurity` by qc-flags).

### 1.5 `eggd_cgp-qc-flags` (converted, v1.0.0)

| Dir | Name | Class | Notes |
|---|---|---|---|
| in | sample_id | string | |
| in | purity_tsv, purity_range_tsv | file | from PURPLE |
| in | sigs_allocation, cup_summary | file | **optional; unlinked in v0.1** |
| out | qc_report | file | `{s}.qc_report.tsv` (16 cols) |

`qc_report.tsv` header (order fixed — plotter reads by name):
`sample_id, flags, purity, ploidy, status, wgd, diploidProportion, ci_width,
bestfit_purity, bestfit_ploidy, bestfit_score, total_snvs, cuppa_top1, cuppa_prob1,
cuppa_top2, cuppa_prob2`. Instance `mem1_ssd1_v2_x2`; timeout 1h.

Flag thresholds: `PURITY_FLOOR` (<0.20), `WGD_SUSPECT` (WGD & <0.35),
`FLAT_GENOME` (diploidProportion>0.85), `WIDE_CI` (maxPurity−minPurity>0.40),
`NO_TUMOR_BESTFIT`, `LOW_SNV_COUNT` (0<SNVs<50).

### 1.6 `eggd_cgp-cnvkit-pon` (reused, v2.0.0) — per-run PoN build (via conductor)

Not a workflow stage. Invoked once per run by `eggd_conductor` when a PoN must be built.

| Dir | Name | Class | Notes |
|---|---|---|---|
| in | coverage_files | array:file | `*.cnn` from per-sample coverage (gathered) |
| in | pon_name | string | optional stem |
| in | fasta, fasta_fai, fasta_gzi | file | optional; GC correction (chr-prefix FASTA) |
| out | cn_reference | file | `*.cnn` → linked into `eggd_atlas_cnv.cnvkit_cn_reference` |
| out | reference_stats | file | summary TSV |

Instance `mem2_ssd1_v2_x4`. (No new reference/resolver app exists — this existing app is
the per-run builder.)

### 1.7 `eggd_cgp-cnvkit-coverage` (reused, v2.0.0) — per-run PoN build input (via conductor)

Not a workflow stage. Run per-sample by conductor to feed the per-run PoN build.
in: tumour_bam, tumour_bai, sample_id, baits, target_avg_size(=250).
out: `coverage_cnn` (`{s}.targetcoverage.cnn`). Instance `mem1_ssd1_v2_x4`.

### 1.8 `eggd_cgp-cnvkit-batch` (reused, v2.0.1)

| Dir | Name | Class | Notes |
|---|---|---|---|
| in | tumour_bam, tumour_bai, sample_id | | |
| in | cn_reference | file | **workflow input** `cnvkit_cn_reference` (supplied, or linked from a per-run PoN) |
| in | baits | file | same BED as PoN |
| in | target_avg_size | int | default 250 |
| in | purity | float | optional; **from PURPLE** (used only if >0.40) |
| in | ploidy | int | optional, default 2; **from PURPLE** |
| in | sample_sex | string | optional; **from PURPLE**; corrects diagram chrX/Y |
| in | drop_low_coverage | boolean | default true |
| out | copy_ratios | file | `{s}.cnr` → plotter `cnvkit_cnr` |
| out | segments | file | `{s}.cns` |
| out | call_segments | file | `{s}.call.cns` → plotter `cnvkit_call_cns` |
| out | scatter_png, diagram_pdf | file | optional |
| out | genemetrics | file | `{s}.genemetrics.tsv` → plotter `cnvkit_genemetrics` (⚠️ ext, §DESIGN 8.1) |

Instance `mem1_ssd1_v2_x4`; timeout 2h.

### 1.9 `eggd_purple_plotter` (reused, v1.0.0)

in: sample_id, purple_tar?, qc_report?, cnvkit_cnr?, cnvkit_call_cns?,
cnvkit_genemetrics?, msi_report?, locus?(=all). out: `igv_html` (`{s}.igv.bars.html`).
Requires ≥1 of purple_tar/cnvkit_cnr. Instance `mem1_ssd1_v2_x2`; timeout 1h. Consumes the
**chr-prefixed** CNVkit/PURPLE files (igv.js/hg38) — never the `eggd_cnv_chr_strip` outputs.

### 1.10 `eggd_cnv_chr_strip` (NEW, v0.1.0) — final stage

| Dir | Name | Class | Notes |
|---|---|---|---|
| in | sample_id | string | output stem |
| in | cnvkit_cnr | file | optional `*.cnr` |
| in | cnvkit_cns | file | optional `*.cns` |
| in | cnvkit_call_cns | file | optional `*.call.cns` |
| in | cnvkit_genemetrics | file | optional `*.genemetrics.tsv` |
| in | purple_cnv_somatic | file | optional `*.purple.cnv.somatic.tsv` |
| in | purple_cnv_gene | file | optional `*.purple.cnv.gene.tsv` |
| out | cnvkit_cnr_nochr | file | `{s}.nochr.cnr` |
| out | cnvkit_cns_nochr | file | `{s}.nochr.cns` |
| out | cnvkit_call_cns_nochr | file | `{s}.nochr.call.cns` |
| out | cnvkit_genemetrics_nochr | file | `{s}.nochr.genemetrics.tsv` |
| out | purple_cnv_somatic_nochr | file | `{s}.purple.cnv.somatic.nochr.tsv` |
| out | purple_cnv_gene_nochr | file | `{s}.purple.cnv.gene.nochr.tsv` |

Strips the `chr` prefix from the `chromosome` column only (`chrM`/`chrMT`→`MT`), leaving
originals untouched. Bundles `atlas_helpers/chr_strip.py`; system `python3`, no execDepends.
Instance `mem1_ssd1_v2_x2`; timeout 1h. All six file inputs and all six `*_nochr` outputs
are `optional: true` — each output is emitted only if its input is supplied (in the workflow
all six inputs are linked, so all six outputs exist).

## 2. DNAnexus resource file IDs

Fill `scripts/resource_ids.env` from this table (verified from the legacy
`cgp-cnv-somatic` `dxworkflow.json` and `cnv-genes-atlas`). Reference genome objects live
in `project-Fkb6Gkj433GVVvj73J7x8KbV` (the `ADMINISTER`/canonical project); resolve every
ID with `dx api file-XXX listProjects` and record as `project:file`.

| Var | Purpose | File ID |
|---|---|---|
| AMBER_JAR | AMBER 4.3-beta.4 | `file-J893948470j1X5zP1gGZYk6G` |
| GERMLINE_SITES | AMBER het PON loci | `file-J88xxvQ4QyV82JKY5PkqXGbv` |
| COBALT_JAR | COBALT 3.0-beta.5 | `file-J893p9Q470j4zY3zzpVBjP11` |
| NORM_FILE | COBALT target-region norm | `file-J89183j4zxxZ5v7Gv5F9BK06` |
| DIPLOID_REGIONS | `DiploidRegions.38.bed.gz` | `file-J88xxvQ4QyV7vyPx294b4Z23` |
| GC_PROFILE | `GC_profile.1000bp.38.cnp` | `file-J88xxvQ4QyVPb8K6VFqX1FKB` |
| COBALT_REF_FASTA | bgzf FASTA (cobalt) | `project-Fkb6Gkj433GVVvj73J7x8KbV:file-GjPxXq84qv8xz3ZV1jFq6z3g` |
| COBALT_REF_FAI | fai | `file-GjPxp3Q4qv8Vk76xFVjfyPgJ` |
| SAGE_JAR | SAGE 5.0-beta.11 | `file-J8F1bFj4gPFvK81Z330Yfk3j` |
| SAGE_REF_FASTA | plain `.fa` | `file-G5xBZvj4yPz1xZfqKKjkG8xQ` |
| SAGE_REF_FAI | fai | `file-G5xF4984yPzJ61ZB1K7gkBXb` |
| HOTSPOTS_VCF | SAGE hotspots | `file-J8F1kQQ4gPFbZQZ1BG7KY62z` |
| HOTSPOTS_TBI | | `file-J8F1kV04gPFz7kZ45XfQPp4K` |
| PANEL_BED | SAGE panel BED | `file-J8F1kv84gPFfQKqYbZVYQqv8` |
| HC_BED | high-confidence BED | `file-J8F1kVj4gPFkBK29yj6P4BgQ` |
| SAGE_PON | SAGE PON `.tsv.gz` | `file-J8F1kX84gPFpJ2kvz13KJ5kP` |
| ENSEMBL_DATA | ensembl tarball | `file-J892FbQ493ZPf7pZzpjX05P6` |
| PURPLE_JAR | PURPLE 4.4-beta.12 | `file-J893gpQ470jJPF18gjJ1zg7Y` |
| PURPLE_TARGET_BED | panel target regions | `file-J88gVF84Y8X123K6JX8jB8Z5` |
| PURPLE_REF_FASTA | bgzf FASTA (purple) | `file-Gb757784XGyY3FPvkPQ74K9z` |
| PURPLE_REF_FAI | fai | `file-Gb7578Q4XGyQ8xfvyxJkBgx7` |
| CNVKIT_BED | CGP dedup BED (20,905 intervals) | `file-J8F94G845FG5ZG55xFzF1fP6` |
| CNVKIT_PON | GC-corrected PoN (41 samples) | `file-J8G49YQ4ZP680j2bpv0F3kFX` |
| CHR_PREFIX_FASTA | chr-prefix FASTA (GC correction) | `project-Fkb6Gkj433GVVvj73J7x8KbV:file-Gb73P8Q4XGyX9QZ0g9qp61xV` |
| CNVKIT_DOCKER | cnvkit 1.0.0 image tarball | `project-Fkb6Gkj433GVVvj73J7x8KbV:file-J8j7Vyj45FG1BbK26JQgQY6q` |

`scripts/resource_ids.env.template`:

```bash
# App IDs (filled after dx build --app)
APP_CHR_PREFIX=app-XXXX
APP_AMBER=app-XXXX
APP_COBALT=app-XXXX
APP_SAGE=app-XXXX
APP_PURPLE=app-XXXX
APP_QC_FLAGS=app-XXXX
APP_CNVKIT_COVERAGE=app-XXXX
APP_CNVKIT_PON=app-XXXX
APP_CNVKIT_BATCH=app-XXXX
APP_PURPLE_PLOTTER=app-XXXX
APP_CNV_CHR_STRIP=app-XXXX
# Reference file IDs (from the table above)
AMBER_JAR=file-...
# ... (all rows) ...
```

## 3. Workflow link schema (`dxworkflow.json`)

Workflow inputs:

| Name | Class | Optional | Note |
|---|---|---|---|
| input_bam | file | no | raw tumour BAM (Ensembl or already-chr); → chr_prefix.input_bam |
| sample_id | string | no | |
| cnvkit_cn_reference | file | **no** | PoN `.cnn`; supplied static, or linked from a per-run `eggd_cgp-cnvkit-pon` job by conductor |
| max_ploidy | int | yes | → purple.max_ploidy |
| ploidy_cap_purity_threshold | float | yes | → purple.ploidy_cap_purity_threshold |
| ploidy_cap_value | int | yes | → purple.ploidy_cap_value (default 2) |

(No `tumour_bai` input — stage 0 generates the index. No `pon_coverage`/`pon_bams` inputs —
PoN building happens once per run in conductor, §4, not inside the workflow.)

Stage link table (`$dnanexus_link`: `workflowInputField` = W, `{stage,outputField}` = S):

| Stage | Input | Source |
|---|---|---|
| chr_prefix | input_bam | W input_bam |
| chr_prefix | mode | fixed `"add_chr"` |
| amber | tumour_bam | S chr_prefix.output_bam |
| amber | tumour_bai | S chr_prefix.output_bai |
| amber | sample_id | W |
| amber | amber_jar / germline_sites | fixed file IDs |
| cobalt | tumour_bam / tumour_bai | S chr_prefix.output_bam / output_bai |
| cobalt | sample_id | W |
| cobalt | cobalt_jar / norm_file / diploid_regions / gc_profile / ref_fasta / ref_fai | fixed |
| sage | tumour_bam / tumour_bai | S chr_prefix.output_bam / output_bai |
| sage | sample_id | W |
| sage | sage_jar / ref_fasta / ref_fai / hotspots_* / panel_bed / hc_bed / pon_file / ensembl_data | fixed |
| purple | sample_id | W |
| purple | amber_tar | S amber.amber_tar |
| purple | cobalt_tar | S cobalt.cobalt_tar |
| purple | somatic_vcf / somatic_vcf_tbi | S sage.somatic_vcf / somatic_vcf_tbi |
| purple | purple_jar / gc_profile / target_regions_bed / ref_fasta / ref_fai / ensembl_data | fixed |
| purple | max_ploidy / ploidy_cap_purity_threshold / ploidy_cap_value | W |
| qc_flags | sample_id | W |
| qc_flags | purity_tsv | S purple.purity_tsv |
| qc_flags | purity_range_tsv | S purple.purity_range_tsv |
| cnvkit_batch | tumour_bam / tumour_bai | S chr_prefix.output_bam / output_bai |
| cnvkit_batch | sample_id | W |
| cnvkit_batch | cn_reference | **W cnvkit_cn_reference** |
| cnvkit_batch | baits | fixed CNVKIT_BED |
| cnvkit_batch | purity | **S purple.purity** |
| cnvkit_batch | ploidy | **S purple.ploidy** |
| cnvkit_batch | sample_sex | **S purple.sample_sex** |
| purple_plotter | sample_id | W |
| purple_plotter | purple_tar | S purple.purple_tar |
| purple_plotter | qc_report | S qc_flags.qc_report |
| purple_plotter | cnvkit_cnr | S cnvkit_batch.copy_ratios |
| purple_plotter | cnvkit_call_cns | S cnvkit_batch.call_segments |
| purple_plotter | cnvkit_genemetrics | S cnvkit_batch.genemetrics |
| cnv_chr_strip | sample_id | W |
| cnv_chr_strip | cnvkit_cnr / cnvkit_cns / cnvkit_call_cns / cnvkit_genemetrics | S cnvkit_batch.copy_ratios / segments / call_segments / genemetrics |
| cnv_chr_strip | purple_cnv_somatic / purple_cnv_gene | S purple.cnv_somatic_tsv / cnv_gene_tsv |

### 3.1 Declared workflow-level outputs

`dxworkflow.json` includes an `outputs` block linking a curated subset of stage outputs to
stable workflow-output names (`outputSource` = `{stage, outputField}`). This is the verified
DNAnexus schema (top-level `outputs` array, each with `name`/`class`/`outputSource`;
workflow-level `inputs` use `workflowInputField`) — see
https://documentation.dnanexus.com/developer/workflows/intro-to-building-workflows
("Locked Workflows"). These are references to stage outputs, not copies. Non-promoted stage
outputs remain accessible as `stage-<id>.<field>`.

| Workflow output | Class | outputSource (stage.field) |
|---|---|---|
| igv_html | file | purple_plotter.igv_html |
| qc_report | file | qc_flags.qc_report |
| purple_cnv_somatic_nochr | file | cnv_chr_strip.purple_cnv_somatic_nochr |
| purple_cnv_gene_nochr | file | cnv_chr_strip.purple_cnv_gene_nochr |
| cnvkit_call_cns_nochr | file | cnv_chr_strip.cnvkit_call_cns_nochr |
| cnvkit_genemetrics_nochr | file | cnv_chr_strip.cnvkit_genemetrics_nochr |
| cnvkit_cnr_nochr | file | cnv_chr_strip.cnvkit_cnr_nochr |
| purple_tar | file | purple.purple_tar |
| purity | float | purple.purity |
| ploidy | int | purple.ploidy |

CNV calls are promoted in **chr-stripped** form (Ensembl names) for downstream apps; the
chr-prefixed originals (`cnvkit_batch.*`, `purple.cnv_*_tsv`, `purple_tar` contents) remain
reachable as `stage-<id>.<field>`.

Output-link syntax in `dxworkflow.json`:

```json
"outputs": [
  { "name": "igv_html",                 "class": "file",
    "outputSource": { "$dnanexus_link": { "stage": "purple_plotter", "outputField": "igv_html" } } },
  { "name": "purple_cnv_somatic_nochr", "class": "file",
    "outputSource": { "$dnanexus_link": { "stage": "cnv_chr_strip", "outputField": "purple_cnv_somatic_nochr" } } },
  { "name": "purity",                   "class": "float",
    "outputSource": { "$dnanexus_link": { "stage": "purple", "outputField": "purity" } } }
]
```

Scalar-link syntax example (typed job output → input). Scalars (`int`/`float`/`string`)
are literal values in the producing job's output hash, linked by JBOR exactly like files;
the downstream input class must match the output class, and the link is validated by
`dx build --workflow`. See DESIGN §2.4 for the full mechanism.

```json
"purity": { "$dnanexus_link": { "stage": "purple", "outputField": "purity" } }
```

## 4. eggd_conductor orchestration (per-run PoN)

The workflow is launched per-sample by `eggd_conductor` from an assay config
`executables` block. Executable keys are app/workflow IDs; each value carries: `analysis`
(`analysis_1`…), `per_sample` (bool), optional `depends_on` (list of analysis labels),
optional `hold` (bool — block conductor until done, required when gathering an output
array), `inputs`, optional `inputs_filter` (regex per input — an **include** filter: only gathered
per-sample outputs whose sample name matches a pattern are kept), and `output_dirs`.
Controls are excluded by giving an include pattern that matches the *wanted* samples;
non-matching samples (e.g. `NA12878`, `Q`-numbers) are simply not gathered.

Input value forms:
- `"INPUT-R1-R2"`, `"INPUT-SAMPLE-NAME"`, `"INPUT-dx_project_name"`, `"INPUT-parent_out_dir"`,
  `"INPUT-SAMPLESHEET"` — conductor-injected tokens.
- `{"$dnanexus_link": {"project": "...", "id": "file-..."}}` — static file.
- `{"$dnanexus_link": {"analysis": "analysis_N", "stage": "stage-...", "field": "..."}}` —
  cross-analysis link (`stage` present only when the source is a *workflow*; omit `stage`
  for a source that is an *app*). Gathers an array when the source analysis is per-sample.

### 4.1 Topology A — PoN provided (single per-sample step)

The atlas workflow adds the `chr` prefix itself (stage 0), so it takes the raw BAM directly.
`cnvkit_cn_reference` is the **workflow-level input** (never a stage-qualified key).
Example file: `conductor/atlas_cnv_pon_provided.example.json`.

```json
{
  "executables": {
    "workflow-ATLAS": {
      "name": "eggd_atlas_cnv", "analysis": "analysis_1", "per_sample": true,
      "process_fastqs": false,
      "inputs": {
        "input_bam": {"$dnanexus_link": {"analysis": "analysis_0", "stage": "stage-align", "field": "bam"}},
        "sample_id": "INPUT-SAMPLE-NAME",
        "cnvkit_cn_reference": {"$dnanexus_link": {"project": "project-XXXX", "id": "file-PON"}},
        "ploidy_cap_purity_threshold": 0.35
      },
      "output_dirs": {"workflow-ATLAS": "/OUT-FOLDER/STAGE-NAME"}
    }
  }
}
```

### 4.2 Topology B — PoN built once per run

Coverage (for the PoN) needs chr-prefixed BAMs, so `eggd_chr_prefix` runs per-sample
*before* coverage; the same prefixed BAM is fed to the atlas workflow (stage 0 passthrough).
Example file: `conductor/atlas_cnv_pon_built.example.json`.

```json
{
  "executables": {
    "app-CHR-PREFIX": {
      "name": "eggd_chr_prefix", "analysis": "analysis_1", "per_sample": true,
      "process_fastqs": false,
      "inputs": {
        "input_bam": {"$dnanexus_link": {"analysis": "analysis_0", "stage": "stage-align", "field": "bam"}},
        "mode": "add_chr"
      },
      "output_dirs": {"app-CHR-PREFIX": "/OUT-FOLDER/APP-NAME"}
    },
    "app-CNVKIT-COVERAGE": {
      "name": "eggd_cgp-cnvkit-coverage", "analysis": "analysis_2", "per_sample": true,
      "process_fastqs": false, "depends_on": ["analysis_1"],
      "inputs": {
        "tumour_bam": {"$dnanexus_link": {"analysis": "analysis_1", "field": "output_bam"}},
        "tumour_bai": {"$dnanexus_link": {"analysis": "analysis_1", "field": "output_bai"}},
        "sample_id": "INPUT-SAMPLE-NAME",
        "baits": {"$dnanexus_link": {"project": "project-XXXX", "id": "file-CNVKIT_BED"}}
      },
      "output_dirs": {"app-CNVKIT-COVERAGE": "/OUT-FOLDER/APP-NAME"}
    },
    "app-CNVKIT-PON": {
      "name": "eggd_cgp-cnvkit-pon", "analysis": "analysis_3", "per_sample": false,
      "process_fastqs": false, "depends_on": ["analysis_2"], "hold": true,
      "inputs": {
        "coverage_files": {"$dnanexus_link": {"analysis": "analysis_2", "field": "coverage_cnn"}},
        "fasta":     {"$dnanexus_link": {"project": "project-Fkb6Gkj433GVVvj73J7x8KbV", "id": "file-CHR_PREFIX_FASTA"}}
      },
      "inputs_filter": {
        "coverage_files": ["^[0-9]+S[0-9]+"]
      },
      "output_dirs": {"app-CNVKIT-PON": "/OUT-FOLDER/APP-NAME"}
    },
    "workflow-ATLAS": {
      "name": "eggd_atlas_cnv", "analysis": "analysis_4", "per_sample": true,
      "process_fastqs": false, "depends_on": ["analysis_3"],
      "inputs": {
        "input_bam": {"$dnanexus_link": {"analysis": "analysis_1", "field": "output_bam"}},
        "sample_id": "INPUT-SAMPLE-NAME",
        "cnvkit_cn_reference": {"$dnanexus_link": {"analysis": "analysis_3", "field": "cn_reference"}},
        "ploidy_cap_purity_threshold": 0.35
      },
      "output_dirs": {"workflow-ATLAS": "/OUT-FOLDER/STAGE-NAME"}
    }
  }
}
```

`analysis_0` (an upstream alignment step) is shown for the raw-BAM source; wire it to match
the host assay, or use `INPUT-*` tokens. The atlas step's `input_bam` reuses the
`analysis_1` chr-prefixed BAM — stage 0 detects it is already prefixed and passes it
through. **`inputs_filter` is an include filter** (eggd_conductor keeps only samples whose
name matches): the pattern `^[0-9]+S[0-9]+` matches real specimen IDs (e.g. `25330S0047`)
so controls (`NA12878`, `Q`-numbers) are not gathered into the PoN — adjust the pattern to
the host assay's sample-ID convention. `hold: true` on the PoN step is required so conductor
waits and gathers the full coverage array before building.

## 5. Environment variables (build-time)

| Var | Required | Default | Description |
|---|---|---|---|
| DX_PROJECT | yes | — | Build project context (`dx select`) |
| (resource_ids.env vars) | yes | — | All app + file IDs from §2 |

No run-time environment variables — all inputs are DNAnexus job inputs.

## 6. External dependencies

| Dependency | Purpose | Version |
|---|---|---|
| dxpy / dx CLI | build & run apps/workflow | current |
| Java (JRE) | AMBER/COBALT/SAGE/PURPLE | openjdk-21 (execDepends) |
| samtools, tabix | BAM/VCF handling | apt (execDepends) |
| R + ggplot2/tidyr/dplyr/cowplot, circos | PURPLE charts | apt (execDepends) |
| CNVkit (Docker image) | coverage/reference/batch | `cgp-cnvkit:1.0.0` |
| igv.js | plotter HTML tracks | 3.8.3 (CDN, view-time) |
| pandas (bundled wheel) | plotter data conversion | offline wheel |
| pytest | helper unit tests | ≥8 |

## 7. Useful links

- PURPLE args: https://github.com/hartwigmedical/hmftools/tree/master/purple (Optional Fitting Arguments)
- eggd_chr_prefix: https://github.com/eastgenomics/eggd_chr_prefix
- eggd_purple_plotter / eggd_cgp-cnvkit-* : github.com/eastgenomics
- eggd_conductor: https://github.com/eastgenomics/eggd_conductor (assay config `executables`, `per_sample`, `depends_on`, `hold`, `inputs_filter`)
- eggd_conductor_configs: https://github.com/eastgenomics/eggd_conductor_configs (Dias CEN/HRD real examples)
- DNAnexus app spec: `dx build --help`, dxapp.json reference

**Scalar outputs / stage linking (DESIGN §2.4):**
- Job Input and Output (JBORs, output hash): https://documentation.dnanexus.com/developer/api/running-analyses/job-input-and-output
- I/O and Run Specifications (output classes, `$dnanexus_link`): https://documentation.dnanexus.com/developer/api/running-analyses/io-and-run-specifications
- Workflows and Analyses (stage input binding): https://documentation.dnanexus.com/developer/api/running-analyses/workflows-and-analyses
- Introduction to Building Workflows (`workflowInputField`, stage links): https://documentation.dnanexus.com/developer/workflows/intro-to-building-workflows
- Path resolution — JBORs: https://documentation.dnanexus.com/user/projects/path-resolution
- Bash apps / `dx-jobutil-add-output`: https://documentation.dnanexus.com/developer/apps/bash

## 8. Glossary

| Term | Meaning |
|---|---|
| PoN | Panel of Normals — CNVkit reference `.cnn` of per-interval median/spread |
| BAF | B-allele frequency (AMBER, per germline het site) |
| WGD | Whole-genome duplication — PURPLE flag; unreliable below ~35% purity |
| ploidy cap | Upper bound on PURPLE's fitted average ploidy (`-max_ploidy`) |
| per-run step | An `eggd_conductor` executable with `per_sample: false` — runs once per sequencing run; its output is reused across all per-sample jobs |
| scalar output | A typed (float/int/string) DNAnexus job output, linkable between stages |
| amplicon mode | CNVkit run with no antitarget bins (targeted panel) |
| chr-prefix | UCSC chromosome naming (chr1…chrM) vs Ensembl (1…MT) |
| nochr | Ensembl-named copy of a CNV file (chr prefix stripped, `chrM`→`MT`), produced by `eggd_cnv_chr_strip` for downstream apps |
| NO_TUMOR | PURPLE status when no tumour signal fits; best-fit recovered from range TSV |
