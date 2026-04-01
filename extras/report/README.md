# One-Off Measurement Plausibility Report Workflow

This folder contains a self-contained reporting workflow built on top of the analysis functions already present in the package `R/` directory. The workflow is intentionally kept outside the package API.

## Contents

- `helpers.R`: small utilities for paths, availability checks, gender mapping, and table I/O
- `sample_andromeda.R`: deterministic synthetic sample data and sample Andromeda creation
- `pipeline_functions.R`: orchestration for preparation, analysis, output writing, and report rendering
- `run_report_pipeline.R`: end-to-end runner
- `quarto/`: Quarto report assets
- `data/`: generated sample data assets
- `output/`: generated tables, figures, metadata, and rendered report

## Requirements

The workflow expects these R packages to be available:

- `ggplot2`
- `isotree`
- `Andromeda` for creating and reading a real Andromeda object
- `quarto` or a working `quarto` CLI installation for rendering the report

If `Andromeda` is unavailable, the sample generator falls back to a deterministic mock object so the pipeline can still be exercised locally.

## Usage

Run from the repository root:

```r
source("extras/report/helpers.R")
source("extras/report/sample_andromeda.R")
source("extras/report/pipeline_functions.R")
source("extras/report/run_report_pipeline.R")

results <- run_report_pipeline(
  baseDir = getwd(),
  createSample = TRUE,
  renderQuarto = FALSE
)
```

To create the sample data only:

```r
source("extras/report/helpers.R")
source("extras/report/sample_andromeda.R")
write_sample_measurement_andromeda(baseDir = getwd())
```

To render the report after outputs have been written:

```r
source("extras/report/helpers.R")
source("extras/report/pipeline_functions.R")
render_report(
  outputDir = file.path(getwd(), "extras", "report", "output"),
  baseDir = getwd()
)
```

## Output Structure

Generated outputs are written under `extras/report/output/`:

- `tables/`: CSV and RDS tables
- `figures/`: PNG figures
- `metadata/`: run metadata and session information
- `report/`: rendered Quarto report

## Notes

- The pipeline assumes one measurement concept and one harmonized unit per run.
- `care_site_id` is retained in the prepared data and row-level outputs but is not used in grouping or contextual modeling.
- Primary stratification is `age_band x gender_concept_id`.
