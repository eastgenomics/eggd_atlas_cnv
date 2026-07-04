#!/bin/bash
# eggd_cgp-qc-flags v1.0.0 (converted from cgp-qc-flags applet: metadata + timeoutPolicy; logic unchanged)
# Post-processes PURPLE output to emit structured QC flags.
# Sources of flags (from REFERENCE.md §11):
#   PURITY_FLOOR      purity < 0.20 (near noise floor of backbone panel)
#   WGD_SUSPECT       WGD=true AND purity < 0.35 (unreliable below ~35%)
#   FLAT_GENOME       diploidProportion > 0.85 (weak CN signal, purity anchor uncertain)
#   WIDE_CI           maxPurity - minPurity > 0.40 (wide confidence interval)
#   NO_TUMOR_BESTFIT  status=NO_TUMOR → surfaces best-fit purity from range TSV
#   CUPPA_DISCORDANT  cuppa top-1 prediction conflicts with low CUP flag (informational)
#   LOW_SNV_COUNT     total allocated SNVs < 50 (SIGS/CUPpa SNV features unreliable)
set -eo pipefail

main() {
    echo "====================================================="
    echo " eggd_cgp-qc-flags: QC flag generation"
    echo " Sample  : ${sample_id}"
    echo "====================================================="

    echo "[1/4] Downloading inputs..."
    dx download "${purity_tsv}"       -o purity.tsv
    dx download "${purity_range_tsv}" -o purity_range.tsv

    [[ -n "${sigs_allocation:-}" ]]  && dx download "${sigs_allocation}" -o sigs.tsv  || touch sigs.tsv
    [[ -n "${cup_summary:-}" ]]      && dx download "${cup_summary}"     -o cuppa.tsv || touch cuppa.tsv

    echo "[2/4] Computing QC flags..."
    python3 - "${sample_id}" purity.tsv purity_range.tsv sigs.tsv cuppa.tsv << 'PYEOF'
import sys, csv

sample_id    = sys.argv[1]
purity_file  = sys.argv[2]
range_file   = sys.argv[3]
sigs_file    = sys.argv[4]
cuppa_file   = sys.argv[5]

flags = []
info  = {}

# ── 1. Read purity TSV ────────────────────────────────────────────────────
with open(purity_file) as f:
    reader = csv.DictReader(f, delimiter='\t')
    rows = list(reader)

if not rows:
    print("ERROR: purity.tsv has no data row", file=sys.stderr)
    sys.exit(1)

d = rows[0]
purity  = float(d['purity'])
ploidy  = float(d['ploidy'])
status  = d['status']
wgd     = d['wholeGenomeDuplication'].upper() == 'TRUE'
dip_prop = float(d.get('diploidProportion', 0))
min_pur  = float(d.get('minPurity', purity))
max_pur  = float(d.get('maxPurity', purity))

info['purity']   = f"{purity:.3f}"
info['ploidy']   = f"{ploidy:.3f}"
info['status']   = status
info['wgd']      = 'TRUE' if wgd else 'FALSE'
info['dip_prop'] = f"{dip_prop:.3f}"
info['ci_width'] = f"{max_pur - min_pur:.3f}"

if purity < 0.20:
    flags.append('PURITY_FLOOR')
if wgd and purity < 0.35:
    flags.append('WGD_SUSPECT')
if dip_prop > 0.85:
    flags.append('FLAT_GENOME')
if (max_pur - min_pur) > 0.40:
    flags.append('WIDE_CI')

# ── 2. NO_TUMOR best-fit from range TSV ──────────────────────────────────
bestfit_purity = None
bestfit_ploidy = None
if status == 'NO_TUMOR':
    flags.append('NO_TUMOR_BESTFIT')
    try:
        with open(range_file) as f:
            rdr = csv.DictReader(f, delimiter='\t')
            best_row = min(rdr, key=lambda r: float(r.get('score', '9999')))
            bestfit_purity = float(best_row['purity'])
            bestfit_ploidy = float(best_row['ploidy'])
            info['bestfit_purity'] = f"{bestfit_purity:.3f}"
            info['bestfit_ploidy'] = f"{bestfit_ploidy:.3f}"
            info['bestfit_score']  = best_row.get('score', 'N/A')
    except Exception as e:
        info['bestfit_error'] = str(e)

# ── 3. SIGS — check SNV count ─────────────────────────────────────────────
total_snvs = 0
try:
    with open(sigs_file) as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for row in rdr:
            try:
                total_snvs += float(row.get('allocation', 0))
            except:
                pass
    info['total_snvs'] = f"{total_snvs:.0f}"
    if 0 < total_snvs < 50:
        flags.append('LOW_SNV_COUNT')
except:
    info['total_snvs'] = 'N/A'

# ── 4. CUPpa — top prediction ─────────────────────────────────────────────
try:
    with open(cuppa_file) as f:
        rdr = csv.DictReader(f, delimiter='\t')
        rows_c = list(rdr)
    combined = [r for r in rows_c if r.get('clf_name') == 'dna_combined']
    if combined:
        r = combined[0]
        info['cuppa_top1']   = r.get('pred_class_1', 'N/A')
        info['cuppa_prob1']  = f"{float(r.get('pred_prob_1', 0)):.3f}"
        info['cuppa_top2']   = r.get('pred_class_2', 'N/A')
        info['cuppa_prob2']  = f"{float(r.get('pred_prob_2', 0)):.3f}"
except:
    info['cuppa_top1'] = 'N/A'

# ── 5. Write report ───────────────────────────────────────────────────────
flag_str = ';'.join(flags) if flags else 'PASS'
print(f"\n{'='*55}")
print(f" QC Report: {sample_id}")
print(f"{'='*55}")
print(f"  Flags     : {flag_str}")
print(f"  Purity    : {info['purity']} [{min_pur:.2f},{max_pur:.2f}]")
print(f"  Ploidy    : {info['ploidy']}")
print(f"  Status    : {info['status']}")
print(f"  WGD       : {info['wgd']}")
print(f"  DipProp   : {info['dip_prop']}")
if bestfit_purity:
    print(f"  BestFit   : purity={info['bestfit_purity']} ploidy={info['bestfit_ploidy']}")
print(f"  SNVs      : {info.get('total_snvs','N/A')}")
print(f"  CUPpa     : {info.get('cuppa_top1','N/A')} (p={info.get('cuppa_prob1','N/A')})")
print(f"{'='*55}\n")

# Write TSV report
out = f"{sample_id}.qc_report.tsv"
header_fields = ['sample_id','flags','purity','ploidy','status','wgd','diploidProportion',
                 'ci_width','bestfit_purity','bestfit_ploidy','bestfit_score',
                 'total_snvs','cuppa_top1','cuppa_prob1','cuppa_top2','cuppa_prob2']
with open(out, 'w') as f:
    f.write('\t'.join(header_fields) + '\n')
    row = [
        sample_id, flag_str,
        info['purity'], info['ploidy'], info['status'], info['wgd'],
        info['dip_prop'], info['ci_width'],
        info.get('bestfit_purity',''), info.get('bestfit_ploidy',''), info.get('bestfit_score',''),
        info.get('total_snvs',''), info.get('cuppa_top1',''), info.get('cuppa_prob1',''),
        info.get('cuppa_top2',''), info.get('cuppa_prob2','')
    ]
    f.write('\t'.join(row) + '\n')

print(f"Written: {out}")
PYEOF

    echo "[3/4] Uploading QC report..."
    QC_FILE=$(ls "${sample_id}.qc_report.tsv")
    qc_report=$(dx upload "${QC_FILE}" --brief)
    dx-jobutil-add-output qc_report "${qc_report}" --class=file

    echo "====================================================="
    echo " eggd_cgp-qc-flags DONE: ${sample_id}"
    echo "====================================================="
}
