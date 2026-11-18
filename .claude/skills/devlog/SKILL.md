---
name: devlog
description: Record a dated session development log — appends a narrative entry to DEVLOG.md, updates CLAUDE.md project history, and commits both with [skip ci]. Use after completing a work session.
---

Record a dated session log for the MPLW515L mouse bone marrow scRNA-seq project. This skill writes a verbose narrative entry to `DEVLOG.md` in the project root, appends a compact subsection to the `## Project development history` section of `CLAUDE.md`, then commits both files with `[skip ci]`.

## Steps

### 1 — Gather context

Run these two commands to understand what changed:

```bash
git -C "$(git rev-parse --show-toplevel)" diff HEAD --stat 2>/dev/null || echo "(no git diff available)"
git -C "$(git rev-parse --show-toplevel)" log --oneline -8 2>/dev/null || echo "(no git log available)"
```

Also review the current conversation to identify:
- What the user asked for
- Which files were edited (`_targets.R`, `src/R/scrna-seq_utilities.R`, `notebooks/<report>/<report>.qmd`, `metadata/project_config.yaml`, `scripts/slurm/`, `src/python/analysis/`, `CLAUDE.md`, etc.)
- Which pipeline targets were added or modified
- Errors encountered and how they were resolved
- Architectural decisions made and the rationale behind them

### 2 — Write DEVLOG.md entry

Read `DEVLOG.md` first (it may already exist). Then **append** (do not overwrite) a new entry using this template. Use today's actual date and time. The horizontal rule `---` separates entries.

```
---

## YYYY-MM-DD HH:MM — <one-line session title>

### What changed
- `_targets.R`: [new or modified targets]
- `src/R/scrna-seq_utilities.R`: [new or modified functions]
- `notebooks/<report>/<report>.qmd`: [sections added/edited]
- `metadata/project_config.yaml`: [parameter changes]
- `CLAUDE.md`: [if updated]
- Other files: [list any others]

### Why
[Scientific or technical rationale — what question drove this session,
what hypothesis or finding motivated the changes]

### Bugs / issues
[Errors encountered, root causes identified, and fixes applied.
If none, write "None."]

### Decisions
[Architectural or analytical choices made this session, with the
reasoning behind each. E.g. "Used DietSeurat after processing to avoid
retaining raw counts layer in downstream targets."]

### Compartment status
- **Stroma** ✅/⏳: [brief update]
- **HSC** ✅/⏳: [brief update]
- **MNC** ✅/⏳: [brief update]

### Next steps
- [ ] [actionable item]
- [ ] [actionable item]
```

### 3 — Update CLAUDE.md

Read `CLAUDE.md` and locate the `## Project development history` section. Append a new `### YYYY-MM-DD — <title>` subsection at the **bottom** of that section. Keep it compact: 3–8 bullet points capturing only the key facts. Match the style of existing subsections already present in the file.

### 4 — Commit

```bash
git add DEVLOG.md CLAUDE.md
git commit -m "docs: session devlog $(date '+%Y-%m-%d') [skip ci]"
```

Report the short commit hash and confirm which files were updated.

## Project context (for writing accurate entries)

- **Pipeline:** R `targets` — analysis targets in `_targets.R`, functions in `src/R/scrna-seq_utilities.R`; it runs as a cluster job (`scripts/slurm/pipeline/run-pipeline.sh`) and renders the reports.
- **Compartments:** Stroma, HSC, MNC — processed via `tar_map` (stroma/hsc/mnc); most HSC analyses use the cell-cycle-regressed `mpn_processed_hsc_fullcc`.
- **Reports:** one Quarto report per folder in `notebooks/` (`quarto_report()` in `_targets.R`), rendered into `output/reports/latest/`; hand-run notebooks in `notebooks/exploration/`, retired ones in `notebooks/archive/`.
- **Config:** `metadata/project_config.yaml` (one config target per section); sample sheet `metadata/samples.csv`.
- **Inputs and outputs:** inputs are links in `data/` (Cell Ranger results `data/cellranger_Rx_<n>`, references from the data hub); outputs in `output/` (`tables/`, `objects/` — `hsc.h5ad` and the Python tools' results).
- **Cluster jobs:** Cell Ranger (`scripts/slurm/cellranger/`), velocity (`scripts/slurm/velocity/`), CellRank / PAGA / pySCENIC / PerturbGen (`scripts/slurm/analysis/`).
- **Known recurring issues:** presto must be installed for FindAllMarkers speed; `SeuratObject::JoinLayers()` required after merge; deMULTIplex2 `$<-` dispatch fails with NAs (use direct slot assignment).
