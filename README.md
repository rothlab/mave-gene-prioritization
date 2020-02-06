# User manual

Author: Kevin Kuang

Last update on: 2020-02-06

## Introduction

This repo contains code to extract missense variants that had been observed in
"clinical testing" (as opposed to "literature only") in ClinVar.

## Criteria

The following selection criteria was applied:

1. missense variants,

2. variants classified to have uncertain significance,

3. variants that are collected through clinical testing and not through
literature curation only.

## Getting started

### Download ClinVar data set

Download the `variant_summary.txt.gz` and `submission_summary.txt.gz` files from the
[NCBI ClinVar FTP Site](ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/).

Unzip downloaded files and leave the text files in the root directory of this repo.

### Filter ClinVar variants

In order to complete this part, you need to install the following R packages:

```{r}
data.table
stringr
magrittr
ggplot2
scales
Cairo
```

Three CSV files will be generated as output:

1. filtered_variants.csv: ClinVar variants that passed all filtering criteria.

2. filtered_submissions.csv: Submitters associated with filtered Clivar variants.

3. num_variants_submitted.csv: The number of filtered variants deposited to
ClinVar by each submitter.

Only the `filtered_variants.csv` file is relevant for downstream analysis
described in the manuscript.

### Pearson correlation of the number of variants deposited between submitters

**This analysis is not in the manuscript.**

In order to complete this part, you need to install the following R packages:

```{r}
data.table
stringr
magrittr
ggplot2
scales
Cairo
viridis
gplots
dendsort
```

A heatmap (`pearson_correlation_heatmap_clustered.png`) will be generated to
compare the pearson correlation of the number of variants deposited between submitters.

If two submitters agree on the number of variants for most of genes in ClinVar,
they will have a high pearson correlation.

A stacked barplot (`submitters_enrichment_in_pearson_correlation.png`) will be
generated to compare the composition of submitters in high PCC and low PCC
categories.
