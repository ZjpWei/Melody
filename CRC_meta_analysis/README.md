# CRC Real Data Meta-analysis

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

## Prediction

## Random-split prediction (RS)
Workflow:
1. Perform meta-analysis on original data using `miMeta` packages.
```console
Rscript Melody_RS.R arg1
```

2. Perform other comparison methods using original data.
```console
Rscript Compared_methods_original_RS.R arg1
```

3. Perform other comparison methods using batch-corrected data.
```console
Rscript Compared_methods_batch_corrected_RS.R arg1
```

4. Summerize the top-rank taxa for all methods
```console
Rscript Top_taxa_RS.R
```

5. Random-split prediction for top-rank taxa
```console
Rscript Prediction_RS.R arg1 arg2 arg3
```
Arguments:
1. `arg1` (`s`): replicate number; from `1` to `100`.
2. `arg2` (`methodology`): method; `Aldex2`, `ANCOM`, `BW`, `clrlasso`, `Melody`.
3. `arg3` (`data.type`): data type; `Original`, `Batch-corrected`.

## Leave-one-study-out prediction (LOSO)
Workflow:
1. Perform meta-analysis on original data using `miMeta` packages.
```console
Rscript Melody_LOSO.R arg1 arg2
```

2. Perform other comparison methods using original data.
```console
Rscript Compared_methods_original_LOSO.R arg1 arg2
```

3. Perform other comparison methods using batch-corrected data.
```console
Rscript Compared_methods_batch_corrected_LOSO.R arg1 arg2
```

4. Summerize the top-rank taxa for all methods
```console
Rscript Top_taxa_LOSO.R
```

5. Random-split prediction for top-rank taxa
```console
Rscript Prediction_LOSO.R arg1 arg2 arg3 arg4
```
Arguments:
1. `arg1` (`s`): replicate number; from `1` to `100`.
2. `arg2` (`ss`): study number for leaving out.
3. `arg3` (`methodology`): method; `Aldex2`, `ANCOM`, `BW`, `clrlasso`, `Melody`.
4. `arg4` (`data.type`): data type; `Original`, `Batch-corrected`.
