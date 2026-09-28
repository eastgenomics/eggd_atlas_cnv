# IMPLEMENTATION — eggd_atlas_cnv (TDD build plan)

> **Historical record.** This milestone plan (M1–M11) describes how the repo was
> originally built, when app sources were copied into a local `apps/` directory and
> built from there. That directory no longer exists in this repo — each app now builds
> and publishes itself from its own per-app GitHub repo (see `scripts/app_ids.json` for
> the name→app-ID mapping), and `tests/test_workflow_json.py` validates `dxworkflow.json`
> against the live, deployed app specs via `dx describe`, not a local copy. Any `apps/...`
> path below is historical narrative, not a path that exists today. See
> `specification/README.md` §"Project layout" for the current structure.

## 0. Prerequisites

- Python 3.12, `python3-venv`, `git`.
- DNAnexus CLI (`dx`) logged in to `org-emee_1`; a build project in `aws:eu-central-1`.
- The legacy applet sources for reference (copy into `apps/`, do not modify in place):
  `/home/wook/Documents/cnv-backbone-purple-atlas/applets/{cgp-amber,cgp-cobalt,cgp-sage,cgp-purple,cgp-qc-flags}`.
- The three existing CNVkit apps and the plotter:
  `/home/wook/Documents/eggd_cgp-cnvkit-{coverage,pon,batch}`, `/home/wook/Documents/eggd_purple-plotter`.
- All reference/JAR file IDs from REFERENCE §2 (fill `scripts/resource_ids.env`).

Milestones **M1, M2, M9 and M10** are pure-Python / JSON and run locally with pytest
(M9 also has a DNAnexus `dx build --workflow` step once the apps exist). Milestones
**M3–M8 and M11** require DNAnexus and real reference data — build/run steps are hand-offs
to a human operator (marked ⚠️ DNAnexus).

## 1. Project scaffold

```
eggd_atlas_cnv/
├── pyproject.toml
├── atlas_helpers/{__init__.py, ploidy_gate.py, purity.py, chr_strip.py}
├── apps/…                       (populated M3–M8; chr-prefix runs separately)
├── dxworkflow.json              (M9)
├── conductor/atlas_cnv_pon_{provided,built}.example.json   (M10)
├── scripts/{resource_ids.env.template, build_all.sh, run_e2e.sh}
└── tests/{test_ploidy_gate.py, test_purity.py, test_chr_strip.py, test_workflow_json.py,
           test_conductor_config.py, fixtures/}
```

`pyproject.toml`:

```toml
[project]
name = "atlas_helpers"
version = "0.1.0"
requires-python = ">=3.12"
dependencies = []

[project.optional-dependencies]
dev = ["pytest>=8"]

[tool.setuptools]
packages = ["atlas_helpers"]
```

## 2. Milestone plan

| M | Module(s) | Red tests written | Green when |
|---|---|---|---|
| M1 | scaffold + `atlas_helpers/__init__.py` | — | `pytest` collects; package imports |
| M2 | `ploidy_gate.py`, `purity.py`, `chr_strip.py` | test_ploidy_gate, test_purity, test_chr_strip | all helper tests green |
| M3 ⚠️ | `eggd_cgp-amber` app | (build/smoke) | `dx build --app` OK; smoke run produces `amber_tar` |
| M4 ⚠️ | `eggd_cgp-cobalt` app | (build/smoke) | builds; smoke produces `cobalt_tar` |
| M5 ⚠️ | `eggd_cgp-sage` app | (build/smoke) | builds; smoke produces `somatic_vcf` |
| M6 ⚠️ | `eggd_cgp-purple` app (+ bundled `ploidy_gate.py`) | test_ploidy_gate (M2) + smoke | builds; conditional re-run works; emits scalar `purity`/`ploidy`/`sample_sex` + `cnv_somatic_tsv`/`cnv_gene_tsv` |
| M7 ⚠️ | `eggd_cgp-qc-flags` app | (build/smoke) | builds; produces `qc_report` |
| M8 ⚠️ | build `eggd_cnv_chr_strip` (+ bundled `chr_strip.py`); rebuild coverage/pon/batch/plotter | test_chr_strip (M2) + smoke | strip writes `*.nochr.*` retaining originals; five apps build; plotter accepts `.tsv` genemetrics |
| M9 | `dxworkflow.json` (8 stages) | test_workflow_json | JSON valid; links resolve; chr-prefixed BAM/BAI are workflow inputs; CNV `*_nochr` outputs from cnv_chr_strip; scalar links + `outputs` block present |
| M10 | `conductor/atlas_cnv_pon_{provided,built}.example.json` | test_conductor_config | both topologies valid; sample_id present; include-filter; links resolve |
| M11 ⚠️ | end-to-end | `scripts/run_e2e.sh` | one sample → non-empty `igv_html` |

## 3. Milestone 1 — scaffold (TDD)

### Red

```python
# tests/test_import.py
def test_package_imports():
    import atlas_helpers  # noqa: F401
```

`.venv/bin/pytest tests/test_import.py` → fails (no package).

### Green

Create `atlas_helpers/__init__.py` (empty) and `pyproject.toml` above; `pip install -e ".[dev]"`.

**Verification:** `.venv/bin/pytest tests/ -v` — collection succeeds, import test green.

## 4. Milestone 2 — decision helpers (TDD)
### Red: `tests/test_ploidy_gate.py`

```python
import pytest
from atlas_helpers.ploidy_gate import decide, PloidyMode, PloidyConfigError

def test_none_mode():
    d = decide(max_ploidy=None, threshold=None, cap_value=2)
    assert d.mode is PloidyMode.NONE
    assert d.first_pass_args == []
    assert d.conditional is False

def test_static_mode():
    d = decide(max_ploidy=3, threshold=None, cap_value=2)
    assert d.mode is PloidyMode.STATIC
    assert d.first_pass_args == ["-max_ploidy", "3"]
    assert d.conditional is False

def test_conditional_mode_and_rerun_boundary():
    d = decide(max_ploidy=None, threshold=0.35, cap_value=2)
    assert d.mode is PloidyMode.CONDITIONAL
    assert d.first_pass_args == []
    assert d.conditional is True
    assert d.needs_rerun(0.34) is True     # below threshold -> cap
    assert d.needs_rerun(0.35) is False    # at threshold -> keep
    assert d.needs_rerun(0.36) is False
    assert d.rerun_args() == ["-max_ploidy", "2"]

def test_mutually_exclusive():
    with pytest.raises(PloidyConfigError):
        decide(max_ploidy=3, threshold=0.35, cap_value=2)

@pytest.mark.parametrize("bad", [0, -1])
def test_invalid_values(bad):
    with pytest.raises(PloidyConfigError):
        decide(max_ploidy=bad, threshold=None, cap_value=2)
    with pytest.raises(PloidyConfigError):
        decide(max_ploidy=None, threshold=0.35, cap_value=bad)

@pytest.mark.parametrize("bad", [-0.1, 0.0, 1.5])
def test_invalid_threshold_range(bad):
    with pytest.raises(PloidyConfigError):
        decide(max_ploidy=None, threshold=bad, cap_value=2)
```

### Red: `tests/test_purity.py`

```python
from pathlib import Path
import pytest
from atlas_helpers.purity import read_purity_ploidy, PurityParseError

FIX = Path(__file__).parent / "fixtures"

def test_reads_conditional_fixture():
    fit = read_purity_ploidy(FIX / "purity_conditional.tsv")
    assert fit.purity == pytest.approx(0.31)
    assert fit.ploidy == pytest.approx(3.9)
    assert fit.status == "NORMAL"
    assert fit.sample_sex == "male"       # 'MALE' -> lowercased
    assert fit.ploidy_int() == 4          # round(3.9)

def test_ploidy_int_floor():
    fit = read_purity_ploidy(FIX / "purity_low_ploidy.tsv")  # ploidy 0.4
    assert fit.ploidy_int() == 1          # never below 1

def test_bad_file_raises(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("purity\tploidy\nNA\tNA\n")
    with pytest.raises(PurityParseError):
        read_purity_ploidy(p)
```

Fixtures — `tests/fixtures/purity_conditional.tsv` (from DESIGN §12) and
`purity_low_ploidy.tsv` (same header, row `0.5 0.9 1.0 0.9 0.4 FEMALE NORMAL … false …`).

### Green: implement the two modules

```python
# atlas_helpers/ploidy_gate.py
from dataclasses import dataclass
from enum import Enum

class PloidyConfigError(ValueError): ...

class PloidyMode(str, Enum):
    NONE = "none"; STATIC = "static"; CONDITIONAL = "conditional"

@dataclass(frozen=True)
class PloidyDecision:
    mode: PloidyMode
    first_pass_args: list[str]
    conditional: bool
    cap_value: int
    threshold: float | None
    def needs_rerun(self, fitted_purity: float) -> bool:
        return self.conditional and self.threshold is not None and fitted_purity < self.threshold
    def rerun_args(self) -> list[str]:
        return ["-max_ploidy", str(self.cap_value)]

def decide(max_ploidy, threshold, cap_value=2) -> PloidyDecision:
    if max_ploidy is not None and threshold is not None:
        raise PloidyConfigError("max_ploidy and ploidy_cap_purity_threshold are mutually exclusive")
    if max_ploidy is not None:
        if max_ploidy < 1:
            raise PloidyConfigError("max_ploidy must be >= 1")
        return PloidyDecision(PloidyMode.STATIC, ["-max_ploidy", str(max_ploidy)], False, cap_value, None)
    if threshold is not None:
        if cap_value < 1:
            raise PloidyConfigError("ploidy_cap_value must be >= 1")
        if not (0.0 < float(threshold) <= 1.0):
            raise PloidyConfigError("ploidy_cap_purity_threshold must be in (0, 1]")
        return PloidyDecision(PloidyMode.CONDITIONAL, [], True, cap_value, float(threshold))
    return PloidyDecision(PloidyMode.NONE, [], False, cap_value, None)
```

```python
# atlas_helpers/purity.py
import csv
import math
from dataclasses import dataclass
from pathlib import Path

class PurityParseError(ValueError): ...

@dataclass(frozen=True)
class PurityFit:
    purity: float; ploidy: float; status: str; sample_sex: str
    def ploidy_int(self) -> int:
        return max(1, round(self.ploidy))

def read_purity_ploidy(path) -> PurityFit:
    rows = list(csv.DictReader(open(Path(path)), delimiter="\t"))
    if not rows:
        raise PurityParseError(f"no data rows in {path}")
    d = rows[0]
    try:
        purity = float(d["purity"]); ploidy = float(d["ploidy"])
    except (KeyError, ValueError) as e:
        raise PurityParseError(str(e)) from e
    if math.isnan(purity) or math.isnan(ploidy):     # blank/NA/NaN => broken PURPLE run
        raise PurityParseError(f"non-finite purity/ploidy in {path}")
    sex = d.get("gender", "").strip().lower()
    return PurityFit(purity, ploidy, d.get("status", ""), sex)
```

### Red: `tests/test_chr_strip.py`

```python
import pytest
from atlas_helpers.chr_strip import strip_chrom, strip_file, ChrColumnError

def test_strip_chrom():
    assert strip_chrom("chr1") == "1"
    assert strip_chrom("chr22") == "22"
    assert strip_chrom("chrX") == "X"
    assert strip_chrom("chrY") == "Y"
    assert strip_chrom("chrM") == "MT"
    assert strip_chrom("chrMT") == "MT"
    assert strip_chrom("7") == "7"          # already Ensembl -> unchanged (idempotent)

def test_strip_file_chromosome_column(tmp_path):
    p, o = tmp_path/"in.tsv", tmp_path/"out.tsv"
    p.write_text("chromosome\tstart\tend\tcopyNumber\nchr7\t1\t100\t3\nchrMT\t1\t50\t2\n")
    assert strip_file(p, o) == 2
    lines = o.read_text().splitlines()
    assert lines[0] == "chromosome\tstart\tend\tcopyNumber"   # header preserved
    assert lines[1].split("\t")[0] == "7"
    assert lines[2].split("\t")[0] == "MT"                    # chrMT -> MT
    assert lines[1].split("\t")[1:] == ["1", "100", "3"]      # other columns unchanged

def test_only_chromosome_column_changes(tmp_path):
    # a gene column with chr-like text must NOT be altered
    p, o = tmp_path/"g.tsv", tmp_path/"g.out.tsv"
    p.write_text("gene\tchromosome\tlog2\nCHR7orf\tchr7\t0.5\n")
    strip_file(p, o)
    row = o.read_text().splitlines()[1].split("\t")
    assert row == ["CHR7orf", "7", "0.5"]

def test_header_only_passthrough(tmp_path):
    p, o = tmp_path/"h.tsv", tmp_path/"h.out.tsv"
    p.write_text("chromosome\tstart\tend\n")
    assert strip_file(p, o) == 0
    assert o.read_text().splitlines() == ["chromosome\tstart\tend"]

def test_no_chromosome_column_errors(tmp_path):
    p, o = tmp_path/"x.tsv", tmp_path/"y.tsv"
    p.write_text("foo\tbar\n1\t2\n")
    with pytest.raises(ChrColumnError):
        strip_file(p, o)
```

### Green: implement `chr_strip.py`

```python
# atlas_helpers/chr_strip.py
from pathlib import Path

class ChrColumnError(ValueError): ...

_CHROM_HEADERS = ("chromosome", "chrom", "#chrom", "contig")

def strip_chrom(chrom: str) -> str:
    c = chrom.strip()
    if not c.startswith("chr"):
        return c                       # already Ensembl-style / not prefixed
    rest = c[3:]
    return "MT" if rest in ("M", "MT") else rest

def _chrom_index(header: list[str], chrom_column):
    if chrom_column and chrom_column in header:
        return header.index(chrom_column)
    for name in _CHROM_HEADERS:
        if name in header:
            return header.index(name)
    return None

def strip_file(in_path, out_path, chrom_column=None) -> int:
    """Rewrite ONLY the chromosome column of a TSV, preserving every other field and the
    original line terminators exactly. Returns the number of data rows rewritten."""
    lines = Path(in_path).read_text().splitlines(keepends=True)
    if not lines:
        Path(out_path).write_text(""); return 0
    idx = _chrom_index(lines[0].rstrip("\r\n").split("\t"), chrom_column)
    if idx is None:
        raise ChrColumnError(f"no chromosome column in {in_path}: {lines[0]!r}")
    out, n = [lines[0]], 0                       # header passed through verbatim
    for line in lines[1:]:
        body = line.rstrip("\r\n"); term = line[len(body):]   # keep original terminator
        fields = body.split("\t")
        if len(fields) > idx and fields[idx]:
            fields[idx] = strip_chrom(fields[idx]); n += 1
        out.append("\t".join(fields) + term)
    Path(out_path).write_text("".join(out))
    return n
```

**Verification:** `.venv/bin/pytest tests/ -v` — all helper tests green.

## 5. Milestones 3–5 — convert AMBER / COBALT / SAGE to apps ⚠️ DNAnexus

For each: copy the legacy applet into `apps/eggd_cgp-<tool>/`, then edit `dxapp.json`:

```jsonc
// apps/eggd_cgp-amber/dxapp.json  (delta from the applet)
{
  "name": "eggd_cgp-amber",
  "title": "eggd_cgp-amber",
  "version": "1.0.0",
  "dxapi": "1.0.0",
  "developers": ["org-emee_1"],
  "authorizedUsers": ["org-emee_1"],
  "runSpec": {
    "interpreter": "bash", "file": "src/code.sh",
    "distribution": "Ubuntu", "release": "24.04", "version": "0",
    "timeoutPolicy": { "*": { "hours": 6 } },
    "execDepends": [
      { "name": "openjdk-21-jre-headless" },
      { "name": "samtools" },
      { "name": "tabix" }
    ]
  },
  "regionalOptions": { "aws:eu-central-1": { "systemRequirements": { "*": { "instanceType": "mem1_ssd1_v2_x4" } } } }
}
```

Then **delete the inline `apt-get install ...` line** from `src/code.sh` (deps now come
from `execDepends`). Do not change any tool flag or output name.

- **M3 (AMBER):** instance `mem1_ssd1_v2_x4`, timeout 6h.
- **M4 (COBALT):** instance `mem2_ssd1_v2_x4`, timeout 6h.
- **M5 (SAGE):** instance `mem2_ssd1_v2_x16`, timeout 8h.

**Verification (each):**
```bash
dx build --app apps/eggd_cgp-<tool> --overwrite
grep -q "apt-get install" apps/eggd_cgp-<tool>/src/code.sh && echo FAIL || echo "no inline apt OK"
# smoke run against a known BAM; confirm the documented output object appears
```

## 6. Milestone 6 — `eggd_cgp-purple` app with ploidy cap + scalar outputs ⚠️ DNAnexus

Copy `cgp-purple` into `apps/eggd_cgp-purple/`. Bundle the helper so `code.sh` can call it:
`cp atlas_helpers/ploidy_gate.py atlas_helpers/purity.py apps/eggd_cgp-purple/resources/home/dnanexus/atlas/`.

Add to `dxapp.json` inputSpec:

```jsonc
{ "name": "max_ploidy", "class": "int", "optional": true,
  "help": "Hard cap passed as PURPLE -max_ploidy. Mutually exclusive with ploidy_cap_purity_threshold." },
{ "name": "ploidy_cap_purity_threshold", "class": "float", "optional": true,
  "help": "If set, PURPLE runs unbounded first; if fitted purity < this, PURPLE re-runs with -max_ploidy ploidy_cap_value." },
{ "name": "ploidy_cap_value", "class": "int", "default": 2, "optional": true,
  "help": "The -max_ploidy value applied when the purity threshold triggers." }
```

Add to outputSpec:

```jsonc
{ "name": "purity", "class": "float" },
{ "name": "ploidy", "class": "int" },
{ "name": "sample_sex", "class": "string" },
{ "name": "cnv_somatic_tsv", "class": "file", "patterns": ["*.purple.cnv.somatic.tsv"] },
{ "name": "cnv_gene_tsv",    "class": "file", "patterns": ["*.purple.cnv.gene.tsv"] }
```

`src/code.sh` implements the §4.2 algorithm. Sketch (essential parts only). Add `jq` to the
PURPLE app's `runSpec.execDepends` (the sketch uses it to read the first-pass args JSON;
it is not on a clean Ubuntu image):

```bash
ATLAS=/home/dnanexus/atlas

# 0. sanitise sample_id before any path/rm use (reject path/shell metacharacters)
case "${sample_id}" in
  *[!A-Za-z0-9._-]* | "" | .* ) echo "ERROR: unsafe sample_id" >&2; exit 1 ;;
esac
WORK="out_${sample_id}"          # work dir derived from a validated id

# 1-3. resolve first-pass args + conditional flag as SHELL-SAFE, delimiter-free values.
#      first_pass_args is emitted as a JSON array and read into a bash array via jq;
#      never pack a multi-word arg list into a single space-split variable.
EVAL=$(python3 - "$ATLAS" <<'PY'
import sys, os, json; sys.path.insert(0, sys.argv[1])
from ploidy_gate import decide
d = decide(max_ploidy=(int(os.environ["max_ploidy"]) if os.environ.get("max_ploidy") else None),
           threshold=(float(os.environ["ploidy_cap_purity_threshold"]) if os.environ.get("ploidy_cap_purity_threshold") else None),
           cap_value=int(os.environ.get("ploidy_cap_value", "2")))
print("FIRST_ARGS_JSON=" + json.dumps(json.dumps(d.first_pass_args)))  # quoted JSON string
print("CONDITIONAL=" + ("yes" if d.conditional else "no"))
PY
)
eval "$EVAL"                                   # sets FIRST_ARGS_JSON (a JSON array string) + CONDITIONAL
mapfile -t FIRST_ARGS < <(jq -r '.[]' <<<"$FIRST_ARGS_JSON")

run_purple() {  # "$@" = extra PURPLE args (may be empty)
    java -Xmx10G -jar purple.jar -tumor "${sample_id}" -amber "${AMBER_DIR}" \
      -cobalt "${COBALT_DIR}" -target_regions_bed target_regions.bed \
      -gc_profile GC_profile.1000bp.38.cnp -ensembl_data_dir "${ENSEMBL_DIR}" \
      -ref_genome ref.fasta -ref_genome_version 38 ${SOMATIC_ARG} "$@" \
      -output_dir "${WORK}/"
}

# 4. pass 1 (FIRST_ARGS is [] in NONE/CONDITIONAL modes, or -max_ploidy N in STATIC)
rm -rf "${WORK}"; run_purple "${FIRST_ARGS[@]}"
PTSV="${WORK}/${sample_id}.purple.purity.tsv"

# 5. conditional re-run (at most one extra pass)
if [ "$CONDITIONAL" = "yes" ]; then
  PURITY=$(python3 -c "import sys;sys.path.insert(0,'$ATLAS');from purity import read_purity_ploidy as r;print(r('$PTSV').purity)")
  NEED=$(python3 -c "print('yes' if float('$PURITY') < float('${ploidy_cap_purity_threshold}') else 'no')")
  if [ "$NEED" = "yes" ]; then rm -rf "${WORK}"; run_purple -max_ploidy "${ploidy_cap_value:-2}"; fi
fi

# 6. emit scalars from the final purity.tsv
eval "$(python3 - "$ATLAS" "$PTSV" <<'PY'
import sys; sys.path.insert(0, sys.argv[1]); from purity import read_purity_ploidy
f = read_purity_ploidy(sys.argv[2])
print(f'FPUR={f.purity}; FPLO={f.ploidy_int()}; FSEX={f.sample_sex or "unknown"}')
PY
)"
tar -czf "${sample_id}.purple.tar.gz" "${WORK}/"
dx-jobutil-add-output purple_tar "$(dx upload "${sample_id}.purple.tar.gz" --brief)" --class=file
dx-jobutil-add-output purity "${FPUR}" --class=float
dx-jobutil-add-output ploidy "${FPLO}" --class=int
dx-jobutil-add-output sample_sex "${FSEX}" --class=string
# surface the CNV call TSVs as standalone outputs (were previously only inside purple_tar).
# Guarantee they always exist: write a header-only file if PURPLE omitted it (some NO_TUMOR fits).
CNV_SOMATIC="${WORK}/${sample_id}.purple.cnv.somatic.tsv"
CNV_GENE="${WORK}/${sample_id}.purple.cnv.gene.tsv"
[ -f "${CNV_SOMATIC}" ] || printf 'chromosome\tstart\tend\tcopyNumber\n' > "${CNV_SOMATIC}"
[ -f "${CNV_GENE}" ]    || printf 'chromosome\tgene\tminCopyNumber\tmaxCopyNumber\n' > "${CNV_GENE}"
dx-jobutil-add-output cnv_somatic_tsv "$(dx upload "${CNV_SOMATIC}" --brief)" --class=file
dx-jobutil-add-output cnv_gene_tsv    "$(dx upload "${CNV_GENE}"    --brief)" --class=file
# (also upload purity_tsv, purity_range_tsv, plots_tar as in the legacy applet)
```

> The header rows above are placeholders — use PURPLE 4.4's real column headers for
> `*.purple.cnv.somatic.tsv` / `*.purple.cnv.gene.tsv` when implementing, so a header-only
> fallback is schema-compatible with a populated file.

> **Build-critical (learned in the v1.0.0 build):** `eggd_purple_plotter` reads THREE files
> out of `purple_tar` — `*.amber.baf.tsv.gz`, `*target_region_cn.tsv`, `*purple.segment.tsv`.
> PURPLE emits the last two but NOT the AMBER BAF. Because `WORK` (`out_${sample_id}`) is a
> separate dir from the AMBER/COBALT extract dir (`${sample_id}`), `tar "${WORK}/"` will NOT
> contain `amber.baf.tsv.gz` and the plotter stage fails with
> `Could not find *amber.baf.tsv in ...purple.tar.gz`. Fix: before tarring, copy the PURPLE
> outputs into the AMBER extract dir and tar THAT (`cp -a "${WORK}/." "${AMBER_DIR}/"; tar
> -czf "${sample_id}.purple.tar.gz" "${AMBER_DIR}/"`), so the archive bundles the BAF too.
> (Keep `WORK` distinct from `${sample_id}` — the two-pass `rm -rf "${WORK}"` must never
> delete the AMBER/COBALT inputs.)

**Verification:**
```bash
.venv/bin/pytest tests/test_ploidy_gate.py tests/test_purity.py -v   # helper logic
dx build --app apps/eggd_cgp-purple --overwrite
# smoke: run with -iploidy_cap_purity_threshold=0.35 on a low-purity sample;
# confirm two PURPLE invocations in the log and integer 'ploidy' scalar output.
dx describe <job> --json | jq '.output | {purity, ploidy, sample_sex, cnv_somatic_tsv, cnv_gene_tsv}'
```

## 7. Milestone 7 — `eggd_cgp-qc-flags` app ⚠️ DNAnexus

Copy applet; add app metadata + `"timeoutPolicy": {"*": {"hours": 1}}`; instance
`mem1_ssd1_v2_x2`. No logic change. Leave `sigs_allocation`/`cup_summary` optional & unlinked.

**Verification:** `dx build --app apps/eggd_cgp-qc-flags --overwrite`; smoke run against a
PURPLE `purity_tsv`/`purity_range_tsv`; confirm `qc_report` with the 16-column header.

## 8. Milestone 8 — stage-0 chr_prefix + `eggd_cnv_chr_strip` + rebuild reused apps ⚠️ DNAnexus

No new reference app is built — the PoN build reuses existing apps (see M10 conductor config).

**Fold in `eggd_chr_prefix` as stage 0 (adapted):** copy the public app into
`apps/eggd_chr_prefix/` and adapt it for single-file per-sample use (DESIGN §4.0):

- inputSpec: `input_bam` (file, `*.bam`). Drop `input_file_array`; fix `mode` to `add_chr`.
- outputSpec: `output_bam` (file, `*.bam`), `output_bai` (file, `*.bai`) — **single files**,
  not arrays, so a stage output links into a single-file downstream input.
- `code.sh`: implement the §4.0 algorithm — reheader when needed, else **passthrough** the
  input BAM; always `samtools index` the emitted BAM. Never emit nothing.

```bash
# add-chr header transform (SN: sequence names): 1..22 -> chr1.., X/Y -> chrX/chrY, MT -> chrM
ADD_CHR='s/\tSN:\([0-9][0-9]*\)\t/\tSN:chr\1\t/g; s/\tSN:X\t/\tSN:chrX\t/g; s/\tSN:Y\t/\tSN:chrY\t/g; s/\tSN:MT\t/\tSN:chrM\t/g'
# passthrough guard (essential logic)
samtools view -H in.bam > orig.sam
sed "$ADD_CHR" orig.sam > new.sam
if cmp -s orig.sam new.sam; then cp in.bam output.bam;              # already chr-prefixed
else samtools reheader new.sam in.bam > output.bam; fi
samtools index output.bam output.bam.bai
```

The transform matches the public `eggd_chr_prefix` mapping (note `MT`→`chrM`, not `chrMT`).
Smoke checks must cover: `1→chr1`, `X→chrX`, `Y→chrY`, `MT→chrM`, and an already-`chr`
header passing through unchanged (`cmp -s` true → copy path).

**Build `eggd_cnv_chr_strip` (new, final stage):** a small app that writes chr-stripped
copies of the CNV text files. Bundle `atlas_helpers/chr_strip.py` into
`apps/eggd_cnv_chr_strip/resources/home/dnanexus/atlas/`. All six file inputs are
`optional: true`, and **all six `*_nochr` outputs are `optional: true`** (a standalone run
may omit some inputs; in the workflow all six are linked so all six outputs exist).
`code.sh` downloads the supplied inputs (`dx-download-all-inputs` puts each at
`~/in/<name>/<file>`), validates `sample_id`, and produces a `*.nochr.*` copy per input:

```bash
set -euo pipefail
ATLAS=/home/dnanexus/atlas
dx-download-all-inputs --parallel
case "${sample_id}" in *[!A-Za-z0-9._-]* | "" | .* ) echo "unsafe sample_id" >&2; exit 1 ;; esac

strip() {  # $1 = input field name, $2 = output name, $3 = dx output field
  local f; f=$(ls ~/in/"$1"/* 2>/dev/null || true); [ -n "$f" ] || return 0
  python3 -c "import sys;sys.path.insert(0,'$ATLAS');from chr_strip import strip_file;strip_file(sys.argv[1],sys.argv[2])" "$f" "$2"
  [ -s "$2" ] || { echo "strip produced no output for $1" >&2; exit 1; }
  dx-jobutil-add-output "$3" "$(dx upload "$2" --brief)" --class=file    # only after success
}
strip cnvkit_cnr          "${sample_id}.nochr.cnr"                    cnvkit_cnr_nochr
strip cnvkit_cns          "${sample_id}.nochr.cns"                    cnvkit_cns_nochr
strip cnvkit_call_cns     "${sample_id}.nochr.call.cns"               cnvkit_call_cns_nochr
strip cnvkit_genemetrics  "${sample_id}.nochr.genemetrics.tsv"        cnvkit_genemetrics_nochr
strip purple_cnv_somatic  "${sample_id}.purple.cnv.somatic.nochr.tsv" purple_cnv_somatic_nochr
strip purple_cnv_gene     "${sample_id}.purple.cnv.gene.nochr.tsv"    purple_cnv_gene_nochr
```

Inputs never modified in place — CNVkit originals stay in `cnvkit_batch`'s outputs, PURPLE
originals in `purple_tar`. Instance `mem1_ssd1_v2_x2`, timeout 1h; no execDepends beyond
system `python3`.

**Rebuild reused apps:**
- Copy `eggd_cgp-cnvkit-coverage`, `eggd_cgp-cnvkit-pon`, `eggd_cgp-cnvkit-batch`,
  `eggd_purple_plotter` into `apps/`; rebuild each with `dx build --app --overwrite`.
- **genemetrics fix (DESIGN §8.1):** widen the plotter's `cnvkit_genemetrics` pattern to
  include `*.genemetrics.tsv`/`*.tsv`.

**Verification:** `.venv/bin/pytest tests/test_chr_strip.py -v`; all six apps build
(`eggd_chr_prefix`, `eggd_cnv_chr_strip`, + four reused); chr_prefix smoke on an
Ensembl-named BAM produces a single chr-prefixed `output_bam` + `output_bai`, and on an
already-chr BAM passes through unchanged; `eggd_cnv_chr_strip` smoke on a `chr`-prefixed
`.call.cns` yields a `*.nochr.call.cns` with `chr` removed and the original untouched;
plotter renders track 7 from a `.genemetrics.tsv`.

## 9. Milestone 9 — assemble `dxworkflow.json` (TDD)

### Red: `tests/test_workflow_json.py`

```python
import json
from pathlib import Path

WF = json.loads(Path("dxworkflow.json").read_text())
STAGE_IDS = {s["id"] for s in WF["stages"]}

def test_nine_stages():
    assert STAGE_IDS == {"chr_prefix","amber","cobalt","sage","purple","qc_flags",
                         "cnvkit_batch","cnv_chr_strip","purple_plotter"}
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
    assert outs["cnvkit_call_cns_nochr"]["stage"]     == "cnv_chr_strip"
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
```

> **Build dependency:** this test loads `apps/*/dxapp.json`, so the M3–M8 apps must be
> scaffolded (at least their `dxapp.json` outputSpecs) before M9 runs — in particular the
> M6 PURPLE app must already declare `cnv_somatic_tsv`/`cnv_gene_tsv`, or
> `test_purple_promoted_cnv_outputs_exist` fails and `dx build --workflow` would later
> reject the `outputs` block. Complete M6 before building the workflow.

### Green: write `dxworkflow.json`

Nine stages (see REFERENCE §3 for the full link table). Key links:

- `chr_prefix.input.input_bam ← workflow input input_bam` (mode fixed to `add_chr`).
- `amber/cobalt/sage/cnvkit_batch.input.tumour_bam ← chr_prefix.output_bam`,
  `.tumour_bai ← chr_prefix.output_bai`.
- `purple.input.amber_tar ← amber.amber_tar`, `.cobalt_tar ← cobalt.cobalt_tar`,
  `.somatic_vcf ← sage.somatic_vcf`; plus `max_ploidy`, `ploidy_cap_purity_threshold`,
  `ploidy_cap_value` from workflow inputs.
- `qc_flags.input.purity_tsv ← purple.purity_tsv`, `.purity_range_tsv ← purple.purity_range_tsv`.
- `cnvkit_batch.input.purity ← purple.purity`, `.ploidy ← purple.ploidy`,
  `.sample_sex ← purple.sample_sex`, `.cn_reference ← workflow input cnvkit_cn_reference`.
  (No `cnvkit_coverage` stage — batch computes its own coverage.)
- `purple_plotter.input.purple_tar ← purple.purple_tar`, `.qc_report ← qc_flags.qc_report`,
  `.cnvkit_cnr ← cnvkit_batch.copy_ratios`, `.cnvkit_call_cns ← cnvkit_batch.call_segments`,
  `.cnvkit_genemetrics ← cnvkit_batch.genemetrics` (chr-prefixed — for igv.js/hg38).
- `cnv_chr_strip.input` ← `cnvkit_batch.copy_ratios/segments/call_segments/genemetrics` +
  `purple.cnv_somatic_tsv/cnv_gene_tsv` + workflow `sample_id`.

Also add a top-level **`outputs`** block (REFERENCE §3.1) promoting the curated
deliverables via `outputSource` links — CNV calls in their **chr-stripped** form:
`igv_html`, `qc_report`, `purple_cnv_somatic_nochr` / `purple_cnv_gene_nochr` /
`cnvkit_call_cns_nochr` / `cnvkit_genemetrics_nochr` / `cnvkit_cnr_nochr`
(all ← `cnv_chr_strip`), `purple_tar`, `purity`, `ploidy`. The chr-prefixed originals stay
accessible as `stage-<id>.<field>`.

Workflow inputs: `input_bam` (**required** raw BAM), `sample_id`, `cnvkit_cn_reference`
(**required** file), `max_ploidy?`, `ploidy_cap_purity_threshold?`, `ploidy_cap_value?`.
No `tumour_bai` input — stage 0 generates the index. No `pon_coverage`/`pon_bams` inputs —
PoN building is a conductor concern (M10).

**Verification:** `.venv/bin/pytest tests/test_workflow_json.py -v`; then
`dx build --workflow . --overwrite` (⚠️ needs the built apps).

## 10. Milestone 10 — conductor executables examples + validation (TDD)

Produce **two** example files (one per PoN topology — kept separate so each is a valid
standalone `executables` block with no duplicate `name`/`analysis` collisions):
`conductor/atlas_cnv_pon_provided.example.json` and `conductor/atlas_cnv_pon_built.example.json`
(DESIGN §13, REFERENCE §4). They are not run against DNAnexus, but their structure and
links are validated.

### Red: `tests/test_conductor_config.py`

```python
import json
from pathlib import Path

PROVIDED = json.loads(Path("conductor/atlas_cnv_pon_provided.example.json").read_text())
BUILT    = json.loads(Path("conductor/atlas_cnv_pon_built.example.json").read_text())

def _by_name(cfg, name):
    return next(v for v in cfg["executables"].values() if v.get("name") == name)

def test_provided_atlas_uses_workflow_input_cn_reference():
    atlas = _by_name(PROVIDED, "eggd_atlas_cnv")
    assert atlas["per_sample"] is True
    assert "cnvkit_cn_reference" in atlas["inputs"]          # workflow-level input name
    assert "stage-cnvkit_batch.cn_reference" not in atlas["inputs"]  # never stage-qualified
    assert "sample_id" in atlas["inputs"]                    # required input present

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
    found = json.dumps(atlas["inputs"])
    assert f'"{pon["analysis"]}"' in found and '"cn_reference"' in found
```

### Green: write the example JSON

Write the two files from REFERENCE §4.1 (provided) and §4.2 (built). Built-PoN chain:
`eggd_chr_prefix` (per_sample) → `eggd_cgp-cnvkit-coverage` (per_sample, `tumour_bam` from
the chr_prefix analysis, `sample_id` set) → `eggd_cgp-cnvkit-pon` (per_run, `depends_on` the
coverage analysis, `hold:true`, gathers `coverage_files` via
`$dnanexus_link {analysis, field:coverage_cnn}` + an **include** `inputs_filter`) →
`eggd_atlas_cnv` (per_sample, `depends_on` the pon analysis, `input_bam` from the chr_prefix
analysis — stage 0 passthrough — `sample_id` set, and `cnvkit_cn_reference` linked from the
pon analysis's `cn_reference`).

**Verification:** `.venv/bin/pytest tests/test_conductor_config.py -v`.

## 11. Milestone 11 — end-to-end ⚠️ DNAnexus

`scripts/run_e2e.sh` runs the workflow on one raw (Ensembl-named) BAM with a supplied PoN
and a conditional ploidy cap — stage 0 adds the `chr` prefix — waits, and asserts a
non-empty `igv_html`.

**Verification:**
```bash
bash scripts/run_e2e.sh 25330S0047
dx describe <analysis>-<plotter-stage> --json | jq '.output.igv_html' # non-null
```

## N. Final checks

```bash
.venv/bin/pytest tests/ -v                       # all local unit + JSON tests green
grep -rL "timeoutPolicy" apps/*/dxapp.json       # expect: none (all have a timeout)
grep -rn "apt-get install" apps/*/src/code.sh    # expect: none (deps in execDepends/asset)
for d in apps/*/; do python3 -c "import json;json.load(open('$d/dxapp.json'));print('$d ok')"; done
python3 -c "import json;json.load(open('dxworkflow.json'));print('workflow ok')"
python3 -c "import json;[json.load(open(f'conductor/{n}')) for n in ('atlas_cnv_pon_provided.example.json','atlas_cnv_pon_built.example.json')];print('conductor ok')"
# confirm no dxpy import in the pure helpers
grep -rn "import dxpy\|from dxpy" atlas_helpers/ && echo FAIL || echo "helpers dxpy-free OK"
```
