# MPLW515L mouse bone marrow scRNA-seq

How does a heterozygous *Mpl* mutation alter the composition and cell identity of the bone marrow
compartments in a mouse model of myeloproliferative neoplasm? Single-cell RNA-seq (10x Chromium,
gene expression plus TotalSeq-B hashtags) of three sorted compartments from WT and Het mice:

| Compartment | Tubes | Conditions |
| --- | --- | --- |
| Stroma | Rx_1, Rx_2 | WT, Het |
| HSC | Rx_3, Rx_4, Rx_5 | WT, Het mCh−, Het mCh+ |
| MNC | Rx_6, Rx_7, Rx_8 | WT, Het mCh−, Het mCh+ |

Within each tube, hashtags (TotalSeq-B 0301–0304) mark individual samples; deMULTIplex2 assigns
them. For Camélia (Psaila lab). The environment, inputs, pipeline and reports follow the project
template: [tools/README.md](tools/README.md).

## Data

All inputs live outside the project and are linked into `data/`; `metadata/catalogue_data.tsv`
records each file (path, size, md5).

- **Raw reads:** per tube one GEX and one CITE (hashtag) library, read from the sequencing delivery
  (X204SC25013626-Z01-F001) in `/well/ludwig/users/vax834/Cbenlabiod/`; each library's folder:
  `metadata/cellranger/libraries.csv`. rio548 also holds tar archives of these same folders in
  `/well/ludwig/users/rio548/data/2026.01.16_Camelia_Magic_Mouse/FastqFiles/`, as a backup copy.
- **Cell Ranger:** `cellranger multi` 10.1.0 per tube (`scripts/slurm/cellranger/`, in a container),
  against GRCm39 + GENCODE vM33 with the mCherry and DRE reporter genes added (built by rio548,
  `/well/ludwig/users/rio548/data/genomes/GRCm39`), hashtags as Antibody Capture
  (`metadata/cellranger/featureref_*.csv`; GFP for WT tubes, MPL for Het tubes). The setup of the
  first processing of these data, which used Cell Ranger 9.0.0. Results in
  `/well/ludwig/users/nrs939/data_hub/mplw515l_camelia_10x/cellranger_multi/`, linked as
  `data/cellranger_Rx_<n>`.
- **Sample sheet:** `metadata/samples.csv`: tube → sample name (short and long).
- **References:** annotation references (Azimuth: Ferchen 2026, Mouse Cell Atlas, HemAtlas; human,
  across species: Zeng 2025, Suo 2022, Jardine 2021), gene sets (Nestorowa, Kucinski, Isobe, Gavish,
  Psaila & Wang, Pan), orthologs (Ensembl, MGI), ImmGen ATAC peaks, the HPA secretome, the Naba
  matrisome, the mm10 cisTarget databases and PerturbGen. Each is a folder in
  `/well/ludwig/users/nrs939/data_hub/`, linked as `data/<name>` by `scripts/sh/link-references.sh`.

## How to run

On the cluster, from the project root. Steps 1–3 are once per machine; the jobs write their logs to
`output/logs/`.

1. **Environment** (the template's three layers; R packages are compiled, so step c takes hours):

   ```bash
   uv lock                                          # a. pyproject.toml lists the analysis's Python packages
   tools/bootstrap.sh                               # b. ./env, .venv, renv, kernels
   sbatch scripts/slurm/setup/install-r-packages.sh # c. the analysis's R packages into renv; writes renv.lock
   tools/verify-env.sh
   ```

   Then record the runtime layer: `conda list -p ./env --explicit --md5 > conda-linux-64.lock`.
2. **Containers** for the tools that don't fit the project's environment (recipes in
   `scripts/containers/`, images built into `output/containers/`):
   `sbatch scripts/slurm/setup/build-container.sh scripts/containers/pyscenic.def`, the same with `velocyto.def`, and for Cell Ranger, after
   downloading its release from 10x,
   `... cellranger --build-arg CELLRANGER_TAR=/path/to/cellranger-10.1.0.tar.gz`.
3. **References:** copy the data hub, then link it: `scripts/sh/link-references.sh` (also copies the
   sample sheet into `metadata/samples.csv`).
4. **Cell Ranger** (~20 h per tube, 4 tubes at a time):
   `sbatch scripts/slurm/cellranger/cellranger-multi.sh`. Check each tube's `web_summary.html`, then
   `scripts/sh/link-cellranger.sh`, `tools/catalogue-data.sh`, review and copy the catalogue to
   `metadata/catalogue_data.tsv`.
5. **Pipeline and reports:** `sbatch scripts/slurm/pipeline/run-pipeline.sh` (`tools/run-pipeline.sh`
   as a cluster job).
6. **Cluster analyses** that start from the pipeline's `output/objects/hsc.h5ad`
   (`scripts/slurm/velocity/`, `scripts/slurm/analysis/`): velocyto → scVelo → CellRank, PAGA,
   pySCENIC, PerturbGen. How to run each: the header of its job file.
7. **Pipeline again:** it renders the reports that read those outputs.

## Results

64 reports, rendered by the pipeline into `output/reports/latest/` with their provenance; dated
snapshots to share are in `output/reports/archive/` (`tools/snapshot-report.sh`). Grouped as in
`_targets.R`:

| Topic | Reports (`notebooks/<name>/`) |
| --- | --- |
| QC and demultiplexing | `first-analysis-milestone-report` (the main report), `qc-exploration`, `demultiplex_before_qc`, `demultiplex_after_qc`, `gene-biotypes` |
| Dimensionality reduction, clustering, cell cycle | `dim-reduction-and-clustering-RNA/-SCT`, `pca-cell-distribution`, `alt-dim-reduction(-hsc)`, `cell-cycle(-hsc)`, `hsc-per-condition-embedding`, `subcluster-exploration`, `supercell-exploration`, `graph-analysis`, `magic-denoising-hsc`, `mpig6b-expression` |
| Annotation and gene sets | `cluster-annotation`, `reference-annotation(-hsc/-mnc/-stroma)`, `geneset-exploration` |
| Differential abundance | `milo-hsc`, `milo-mnc`, `milo-stroma`, `miloR_to_genes` |
| Gene regulation, communication | `hdwgcna-hsc*` (8), `scenic-hsc(-conditions)`, `upstream-regulators-hsc`, `cellchat-wt/-mch-neg/-mch-pos` |
| Potency and trajectories | `differentiation_potential`, `trajectory-slingshot*` (12), `transition-matrices`, `paga-hsc`, `paga-traj-analysis-*` (3), `dpt-traj-analysis`, `cellrank-hsc*` (3) |

Also in `notebooks/`:

- `exploration/` holds notebooks run by hand: the daily notebooks, meeting notes, the Jupyter
  notebooks (DPT, PAGA, velocity, AlphaGenome) and the metacell notebooks (MC2, SEACells).
- `archive/` holds reports that read targets the pipeline no longer has (each listed in its README),
  and the first version of the report as a book.

## People

- Analysis: Mohamed R. Shoeb (mohamed.shoeb@ludwig.ox.ac.uk)
- Collaborator: Camélia (Psaila lab)
- First Cell Ranger processing and the GRCm39 + reporters reference: rio548
