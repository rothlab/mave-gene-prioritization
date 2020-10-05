# User manual

Author: Da (Kevin) Kuang

Last update on: 2020-10-05

## Introduction

This repo contains code to:
1. Extract missense variants that had been observed in
"clinical testing" (as opposed to "literature only") in ClinVar.

2. Calculate movability and reappearance parameters using the aggregated and capped Invitae variant counts.

3. Apply the movability and reappearance parameters to ClinVar genes, calculating their movability- and reappearance-weighted impact score (MARWIS), as well as their difficulty-adjusted impact score (DAIS).

4. Rank ClinVar genes based on
   1. Unique number of VUS in ClinVar,
   2. Mobability- and Reappearance-weighted Impact Score (MARWIS), and
   3. Difficulty-ajusted Impact Score (DAIS).

## Selection Criteria

The following selection criteria was applied:

1. Missense variants,

2. Variants classified to have uncertain significance (i.e. VUS), and

3. Variants that are collected through clinical testing and not through
literature curation only.

## Getting started

### Download ClinVar data set

**Please note: if `filtered_variants.csv` and `missense_variants.csv` files exist, the filtering process will be skipped.**

Download the [variant_summary.txt.gz](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2020-08.txt.gz) and [submission_summary.txt.gz](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/submission_summary_2020-08.txt.gz) files from the NCBI ClinVar FTP Site.

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

Two CSV files are needed as input:

1. `invitae_variant_count.csv`: Aggregated per-gene level variant counts where the occurrence of each variant was capped at 7 (the same capping threshold used and documented in the manuscript).

2. `gene_length.csv`: Protein-coding gene's proten length (i.e. number of amino acids) based on the canonical isoform according the the Ensembl database. 

Three CSV files may be generated as output:

1. `missense_variants.csv`: Missense ClinVar variants.
*This file will only be generated if it doesn't already exist.*

2. `filtered_variants.csv`: ClinVar variants that passed all filtering criteria. *This file will only be generated if it doesn't already exist.*

3. `ranked_clinvar_genes.csv`: ranked ClinVar genes.

## LICENSE

### MIT

Copyright 2020 Kevin Kuang

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.