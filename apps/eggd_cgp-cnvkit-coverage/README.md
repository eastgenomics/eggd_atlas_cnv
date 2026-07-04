<!-- dx-header -->
# eggd_cgp-cnvkit-coverage (DNAnexus Platform App)

Computes CNVkit on-target read coverage for a single tumour-only CGP panel sample.
Runs as **Step 1** of the three-app CGP CNVkit pipeline:

```
eggd_cgp-cnvkit-coverage  (×N, parallel) → eggd_cgp-cnvkit-pon (×1) → eggd_cgp-cnvkit-batch (×N, parallel)
```

## What does this app do?

Runs `cnvkit.py coverage` in **amplicon mode** against a gene panel BED file, producing
a `.targetcoverage.cnn` file containing per-interval read depth for one sample.

Amplicon mode is used (no antitarget bins) because this panel is co-captured with a
genome-wide CNV backbone (~57,000 SNP probes at 50 kb spacing). Off-target reads from the
gene panel fall partly on backbone probe positions and would produce inflated, unreliable
antitarget depth. The backbone genome-wide signal is handled by a separate tool (PURPLE).

CNVkit is installed from PyPI at job start via a virtual environment to avoid conflicts
with Ubuntu 24.04 system-managed packages.

## What are the typical use cases for this app?

- Building a panel-of-normals (PoN) reference for a new assay cohort: run this app
  on all normal or tumour samples, then pass the `.cnn` files to `eggd_cgp-cnvkit-pon`
- Per-sample analysis: run alongside `eggd_cgp-cnvkit-batch`, which also runs its own
  coverage step internally — this app is only needed separately when coverage computation
  and analysis are decoupled (e.g. for PoN construction)

**Assay context:** Designed for the Twist Oncology DNA CGP panel (729 genes, ~1.93 Mb,
20,912 exonic intervals) combined with the East Genomics CNV backbone. The BED file must
be deduplicated on coordinates before use — the Twist CGP BED contains 7 pairs of
intervals with identical chr:start:end but different gene names (MSI microsatellite
markers sharing coordinates with primary gene entries). Use `drop_duplicates(subset=[chr,start,end])`
to reduce 20,912 → 20,905 intervals before uploading to DNAnexus.

## What are the inputs?

| Input | Class | Required | Description |
|---|---|---|---|
| `tumour_bam` | file | ✅ | Tumour BAM (chr-prefix chromosome names) |
| `tumour_bai` | file | ✅ | BAM index (.bai or .bam.bai) |
| `sample_id` | string | ✅ | Output file stem: `<sample_id>.targetcoverage.cnn` |
| `baits` | file | ✅ | Gene panel BED (chr-prefix, deduplicated; 20,905 intervals for Twist CGP) |
| `target_avg_size` | int | ➖ | Target bin size in bp (default: 250). Must match across coverage, PoN, and batch runs. |

## What are the outputs?

| Output | Class | Description |
|---|---|---|
| `coverage_cnn` | file | Per-interval read depth in CNVkit format (`<sample_id>.targetcoverage.cnn`) |

The `.targetcoverage.cnn` file contains one row per BED interval with columns:
`chromosome`, `start`, `end`, `gene`, `depth`, `log2`.

## How to run this app from the command line?

```bash
dx run eggd_cgp-cnvkit-coverage \
  -itumour_bam=file-xxxx \
  -itumour_bai=file-xxxx \
  -isample_id="25330S0047" \
  -ibaits=file-xxxx \
  -itarget_avg_size=250 \
  --destination "project-xxxx:/cnvkit/coverage/25330S0047/" \
  --instance-type mem1_ssd1_v2_x4 \
  --priority high \
  -y
```

Typical runtime: **10–15 min** on `mem1_ssd1_v2_x4`.
All samples in a cohort can be submitted simultaneously (one job per sample).

## Dependencies

CNVkit is installed from PyPI at job start:
```
python3 -m venv /tmp/cnvkit-env
/tmp/cnvkit-env/bin/pip install cnvkit
```
CNVkit version in use: **0.9.13** (latest stable at time of deployment).
No R is required for this step.

## Notes

- **Chromosome naming:** BAM and BED must use consistent chromosome naming (both
  chr-prefix or both no-chr). This app is designed for chr-prefix BAMs.
- **Antitarget bins:** Deliberately omitted. See design rationale above.
- **GC correction:** Not applied at this step. The PoN construction step
  (`eggd_cgp-cnvkit-pon`) accepts an optional `--fasta` for GC/repeat annotation;
  omit it if your reference FASTA uses different chromosome naming from the BAM/BED.
