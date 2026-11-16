# notebooks/exploration/

Notebooks run by hand, from the project root (JupyterLab: `env/bin/jupyter lab`; R: an R session in
the project). The pipeline doesn't render them and Quarto ignores this folder (`_quarto.yml`).

| File | What |
| --- | --- |
| `daily/` | the dated lab-notebook entries (2026-05 to 2026-08) and meeting notes |
| `dpt_hsc.ipynb` | diffusion pseudotime on the HSC compartment: writes `output/objects/dpt/` (`hsc_dpt.h5ad`, `dpt_cell_table.csv`), read by the `dpt-traj-analysis` report and by CellRank with `INPUT_TYPE=dpt` |
| `paga_hsc.ipynb`, `paga_hsc_conditions.ipynb` | PAGA, interactively (`src/python/analysis/run_paga.py` is the batch version) |
| `velocity-hsc.ipynb` | scVelo + CellRank, interactively |
| `alphagenome_calr_klf1.ipynb` | AlphaGenome predictions (CALR, KLF1): needs `uv sync --group alphagenome` and an API key (`ALPHAGENOME_API_KEY`) |
| `mc2-hsc.qmd`, `seacells-hsc.qmd` | metacells with MC2 and SEACells (Python; render with the project's Python kernel) |
| `cellchat.qmd`, `nichenet.qmd`, `computational_meeting_17_06_26.qmd`, `mpn_meeting_presentation.qmd` | scratch notebooks |
| `Untitled.ipynb` | scratch: the first pySCENIC attempt |
| `Untitled1.ipynb` | empty |

Several were written for an interactive session that already held objects (e.g.
`mpn_processed_hsc_fullcc`), so they don't run top to bottom on their own; three daily notebooks
(`2026_05_14.R`, `2026_05_19.qmd`, `2026_07_29.qmd`) contain pasted text that isn't R code. The
Jupyter notebooks keep the outputs of their last run, which still show the old paths.
