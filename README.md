# User manual

Author: Da (Kevin) Kuang

Last update on: 2020-08-31

## Introduction

This repo contains code to:
1. Extract missense variants that had been observed in
"clinical testing" (as opposed to "literature only") in ClinVar.

2. Rank genes based on their number of unique variants of uncertain significance (VUS).

## Selection Criteria

The following selection criteria was applied:

1. missense variants,

2. variants classified to have uncertain significance,

3. variants that are collected through clinical testing and not through
literature curation only.

## Getting started

### (Optional) Download ClinVar data set

**Please note: if a `filtered_variants.csv` file exists, the filtering process will be skipped.**

A list of filtered variants (`filtered_variants.csv`) is provided for users to reproduce the ranked list reported in the manuscript.

However, because new variants are added to the ClinVar database regularly, we recommend users to download the most recent ClinVar dataset.

Download the [variant_summary.txt.gz](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz) and [submission_summary.txt.gz](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz) files from the NCBI ClinVar FTP Site.

Please refer to [this file](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/README) for headers of both files.

Unzip downloaded files and leave the text files in the root directory of this repo.

Please make sure the plain text files are named as `variant_summary.txt` and `submission_summary.txt`.

### Install packacges

In order to run the scripts below, you need to install these R packages:

```{r}
data.table
stringr
```

### Rank genes by unique ClinVar variants

**Script: `rankClinVarGenes.R`**

You may run the script by `source()`-ing it to 
an interactive R session or by executing the following command:

```{bash}
Rscript rankClinVarGenes.R
```

Two CSV files may be generated as output:

1. filtered_variants.csv: ClinVar variants that passed all filtering criteria. *This file will only be generated if it doesn't already exist.*

2. genes_ranked_filtered_variants.csv: genes ranked by their unqiue ClinVar variants.

## LICENSE

MIT