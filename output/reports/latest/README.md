# output/reports/latest/

The **latest render** of each report, written here by the pipeline (`tar_quarto()` in `_targets.R`).

- One folder per report, mirroring its source: `notebooks/<report>/<report>.qmd` ->
  `latest/notebooks/<report>/` with `<report>.html`, `<report>_files/` and `artifacts/` (figures and
  tables the report saved with `report_dir()`).
- Paths are fixed so targets can tell whether a report is up to date; each render replaces the previous one.
- Git-ignored (except this README). To keep a version, see `../archive/`.
