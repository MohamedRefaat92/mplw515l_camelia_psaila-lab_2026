#!/usr/bin/env bash
# Download a 10x Genomics prebuilt reference and unpack it into the data hub's raw references
# (references with added sequences, e.g. reporter genes, go to reference_genomes/modified/ instead).
#
#   scripts/sh/get-10x-reference.sh refdata-gex-GRCm39-2024-A
#
# → /well/ludwig/users/nrs939/data_hub/10x_resources/reference_genomes/raw/refdata-gex-GRCm39-2024-A/
# Skipped if that folder exists; an interrupted download resumes on the next run.
set -euo pipefail
cd "$(dirname "$0")/../.."
name=${1:?usage: scripts/sh/get-10x-reference.sh <reference name, e.g. refdata-gex-GRCm39-2024-A>}
dest=/well/ludwig/users/nrs939/data_hub/10x_resources/reference_genomes/raw

if [[ -d $dest/$name ]]; then
  echo "[SKIP] $dest/$name exists"
  exit 0
fi
mkdir -p "$dest"
curl -fL -C - --retry 5 -o "$dest/$name.tar.gz.part" "https://cf.10xgenomics.com/supp/cell-exp/$name.tar.gz"
tar -xzf "$dest/$name.tar.gz.part" -C "$dest"
rm "$dest/$name.tar.gz.part"   # the unpacked folder is what Cell Ranger reads
ls "$dest/$name"
