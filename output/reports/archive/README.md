# output/reports/archive/

**Frozen snapshots** of rendered reports, never modified after they are made:

    <YYYY-MM-DD>/<report>_<YYYY-MM-DD_HH-MM>/

`tools/run-pipeline.sh` snapshots every report it re-rendered; `tools/snapshot-report.sh <report>`
takes one by hand. Each snapshot is a copy of the report's folder in `../latest/` (HTML, `<report>_files/`,
`artifacts/`), and its provenance block records the commit, input data and pipeline it came from.

Share a snapshot by zipping **the whole folder**: the HTML links to the files next to it.
Git-ignored (except this README).
