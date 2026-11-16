# metadata/

Metadata you write about the project and its data: the sample sheet (`samples.csv`), annotations, gene lists,
notes on sources. These are real files, not links, and the pipeline can read them like any input
(e.g. `tar_target(samples_file, "metadata/samples.csv", format = "file")`).
Everything here is committed, except files starting with `auto_`, which scripts generate (git-ignored).

## The data catalogue

| File | Made by | In git |
| --- | --- | --- |
| `auto_catalogue_data.tsv` | `tools/catalogue-data.sh`, from `data/` | no |
| `catalogue_data.tsv` | you, by reviewing the auto file and copying it | yes |

```bash
tools/catalogue-data.sh                                          # regenerate; reports new / changed / missing files
cp metadata/auto_catalogue_data.tsv metadata/catalogue_data.tsv    # accept it after review
```

The generated file starts from your reviewed copy: your descriptions (kind, sample, organism, tissue, ...)
and any columns you added are kept; file facts (path, format, size, owner, md5, ...) are refreshed; new
files are added; files no longer in `data/` are kept and marked `missing`. Describe new rows in
`catalogue_data.tsv`, then commit it.

## In this project

| File | What |
| --- | --- |
| `project_config.yaml` | the analysis parameters: QC thresholds per compartment, processing, sub-clusters, cell-type markers and colours; read by `_targets.R` (one target per section) |
| `samples.csv` | the sample sheet: `Tube`, `sample_name_short`, `sample_name_long` (copied from the data hub by `scripts/sh/link-references.sh`) |
| `cellranger/libraries.csv` | per tube: its feature reference, and the GEX and CITE libraries (FASTQ id, and the library's folder below `$FASTQ_ROOT`) |
| `cellranger/featureref_GFP.csv`, `featureref_MPL.csv` | the TotalSeq-B hashtags (0301–0304) for the WT (GFP) and Het (MPL) tubes |
