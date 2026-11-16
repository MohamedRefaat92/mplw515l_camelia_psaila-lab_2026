#!/usr/bin/env bash
# Link the reference datasets the analysis reads (annotation references, gene sets, orthologs,
# cisTarget, PerturbGen) from the data hub into data/, each as data/<name>, and copy the sample sheet
# into metadata/samples.csv if it isn't there yet. A dataset missing from the hub is reported.
#
#   scripts/sh/link-references.sh [hub]     (default: /well/ludwig/users/nrs939/data_hub)
set -euo pipefail
cd "$(dirname "$0")/../.."
hub=${1:-/well/ludwig/users/nrs939/data_hub}

# Read by the pipeline (_targets.R)
pipeline=(
  ensembl mgi marker_list_claude                                 # orthologs, marker list
  ferchen_2026_nat_immunol mca_mouse_cell_atlas hematlas         # Azimuth references (mouse)
  zeng_2025_blood_cancer_discov suo_2022 jardine_2021            # references projected across species
  nestorowa_2017_blood kucinski_2024_cell_stem_cell isobe_2023_cell_genom gavish_2023_nature
  psaila_wang_2020_molcell pan_2025_cell_rep                     # gene sets
  yoshida_2020_cell human_protein_atlas naba_2012_mcp            # ATAC peaks, secretome, matrisome
)
# Read by the cluster jobs: SCENIC (mm10 cisTarget) and PerturbGen
jobs=(scenic perturbgen)
# Read only by notebooks in notebooks/exploration/daily/
exploration=(wang_2021_cell_stem_cell bandyopadhyay_2024_cell sanjuan-pla_2013 chip_atalas STEAM-v1)

to_link=() missing=()
for name in "${pipeline[@]}" "${jobs[@]}" "${exploration[@]}"; do
  if [[ -L data/$name ]]; then
    continue
  elif [[ -e $hub/$name ]]; then
    to_link+=("$hub/$name")
  else
    missing+=("$name")
  fi
done
if [[ ${#to_link[@]} -gt 0 ]]; then
  tools/link-input.sh "${to_link[@]}" || true   # it reports and skips anything it can't link
fi
echo "linked in data/: $(find data -maxdepth 1 -type l | wc -l)"

# The sample sheet (which tube is which sample) is metadata: a committed copy in metadata/
sheet=$hub/mpn_10x_camelia/tube_sample_map.csv
if [[ -f metadata/samples.csv ]]; then
  echo "metadata/samples.csv exists (not overwritten)"
elif [[ -f $sheet ]]; then
  cp "$sheet" metadata/samples.csv
  echo "copied $sheet → metadata/samples.csv (commit it)"
else
  missing+=("mpn_10x_camelia/tube_sample_map.csv (→ metadata/samples.csv)")
fi

if [[ ${#missing[@]} -gt 0 ]]; then
  echo "not in $hub:" >&2
  printf '  %s\n' "${missing[@]}" >&2
  exit 1
fi
