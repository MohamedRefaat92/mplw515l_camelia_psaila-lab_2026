# scripts/containers/

Apptainer / Singularity definition files (`.def`) for tools conda doesn't provide (e.g. Cell Ranger).

The `.def` recipe is committed; the built image (`.sif`) is large and machine-specific, so it is
git-ignored. Build with: `apptainer build <name>.sif scripts/containers/<name>.def`
