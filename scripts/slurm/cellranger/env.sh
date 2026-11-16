# shellcheck shell=bash
# Settings for the Cell Ranger step, shared by scripts/slurm/cellranger/cellranger-multi.sh,
# scripts/slurm/velocity/velocyto.sh and scripts/sh/link-cellranger.sh.
# To change one, export it before submitting (e.g. CR_OUT=/elsewhere sbatch ...).

# Cell Ranger version: the release in output/containers/cellranger.sif (the job checks the image is it)
export CR_VERSION=${CR_VERSION:-10.1.0}
# The sequencing provider's delivery, read where it was delivered; each library's folder below it
# (X204SC.../01.RawData/<library>) is in metadata/cellranger/libraries.csv. rio548's tar archives of
# these same folders (data/2026.01.16_Camelia_Magic_Mouse/FastqFiles) are a backup, not an input.
export FASTQ_ROOT=${FASTQ_ROOT:-/well/ludwig/users/vax834/Cbenlabiod}
# The reference, and the two things scripts/slurm/cellranger/mkref.sh builds it from:
# the 10x prebuilt mouse reference (GRCm39 Ensembl 110 + GENCODE vM33; scripts/sh/get-10x-reference.sh)
export CR_REF_SRC=${CR_REF_SRC:-/well/ludwig/users/nrs939/data_hub/10x_resources/reference_genomes/raw/refdata-gex-GRCm39-2024-A}
# and the reporter cassette as one extra contig: $CR_REPORTER.fa with the matching .gtf, which is
# the transduced LeGO vector's transcript (mCherry, DRE and the WPRE-LTR 3' UTR, where the reads are)
export CR_REPORTER=${CR_REPORTER:-/well/ludwig/users/nrs939/data_hub/10x_resources/reference_additions/LeGO_SFFV_mCherry_DRE_transcript}
export CR_REF=${CR_REF:-/well/ludwig/users/nrs939/data_hub/10x_resources/reference_genomes/modified/GRCm39-2024-A_LeGO_SFFV_mCherry_DRE}
# Results: one folder per tube (outside the project; linked into data/ by scripts/sh/link-cellranger.sh)
export CR_OUT=${CR_OUT:-/well/ludwig/users/nrs939/data_hub/mplw515l_camelia_10x/cellranger_multi}

export CR_LIBRARIES=metadata/cellranger/libraries.csv
