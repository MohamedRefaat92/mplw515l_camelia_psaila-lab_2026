# Project title

> Made from the [hyprid_proj](https://github.com/MohamedRefaat92/hyprid_proj) template. How the environment,
> inputs, pipeline and reports work: [tools/README.md](tools/README.md). In a new project, replace this file's
> placeholders (and this note) with the project's own description.

One paragraph: the question, the samples, who it is for.

## Data

Where the inputs come from and how they were produced. They are linked into `data/` and described
in `metadata/catalogue_data.tsv`.

## How to run

```bash
tools/bootstrap.sh       # once per machine: build the environment
tools/verify-env.sh      # check it
tools/run-pipeline.sh    # run the analysis and render the reports into output/reports/latest/
```

Project-specific steps (cluster jobs, containers, other scripts) and the order to run them in.

## Results

Which reports answer which question; where shared snapshots are (`output/reports/archive/`).

## People

Who runs the analysis, who the collaborators are, who to ask.
