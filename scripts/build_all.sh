#!/usr/bin/env bash
# build_all.sh — dx build --app for every app, in dependency order, then the workflow.
# Requires: dx login + a build project selected (dx select), and scripts/resource_ids.env.
set -euo pipefail

cd "$(dirname "$0")/.."

# Apps in dependency order (converted + new first, then reused, plotter last).
APPS=(
  eggd_cgp-amber
  eggd_cgp-cobalt
  eggd_cgp-sage
  eggd_cgp-purple
  eggd_cgp-qc-flags
  eggd_cgp-cnvkit-coverage
  eggd_cgp-cnvkit-pon
  eggd_cgp-cnvkit-batch
  eggd_cnv_chr_strip
  eggd_purple_plotter
)

for app in "${APPS[@]}"; do
  echo "=== dx build --app apps/${app} ==="
  dx build --app "apps/${app}" --overwrite
done

echo "=== dx build --workflow . ==="
dx build --workflow . --overwrite

echo "All apps + workflow built."
