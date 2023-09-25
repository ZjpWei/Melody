# CRC Real Data Analysis

## Data analysis

Scripts:

1. `Melody_abess.R` performs meta-analysis on CRC original data using the `miMeta` package.

2. `train_model.R` analyzes CRC original data using other comparison methods.

3. `train_model_batch.R` analyzes CRC batch-corrected data using other comparison methods.

Arguments:

1. arg1 (`s`): Data scenario.

2. arg2 (`tag`): Folder name for analysis.

Examples: 

1. To run `Melody_abess.R` for data in `CRC_all_order` folder with `all` data.
```console
Rscript Melody_abess.R all CRC_all_order
```

2. To run `train_model.R` for data in `CRC_loso_order` folder with leaving out study `1`.
```console
Rscript train_model.R 1 CRC_loso_order
```

## Sensitivity: meta-analysis with different reference taxa

1. `sensitivity.R` performs meta-analysis on CRC original data using `miMeta` package with 5 different reference taxa.
