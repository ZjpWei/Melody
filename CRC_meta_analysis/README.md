# CRC Real Data Analysis

## Prepare data

Excute `Prepare_data.R` to generate data.

## Data analysis

Scripts:

1. `Melody.R` performs meta-analysis on CRC original data using the `miMeta` package.

2. `Compared_methods_original.R` analyzes CRC original data using other comparison methods.

3. `Compared_methods_batch_corrected.R` analyzes CRC batch-corrected data using other comparison methods.

Arguments:

1. arg1 (`s`): Data scenario.

2. arg2 (`tag`): Folder name for analysis.

Examples: 

1. To run `Melody.R` for data in `CRC_all_order` folder with `all` data.
```console
Rscript Melody.R all CRC_all_order
```

1. To run `Compared_methods_original.R` for original data in `CRC_loso_order` folder with leaving out study `1`.
```console
Rscript Compared_methods_original.R 1 CRC_loso_order
```

## Sensitivity: meta-analysis with difference reference taxon

1. `sensitivity.R` performs meta-analysis on CRC original data using `miMeta` package with 5 different reference taxa.
