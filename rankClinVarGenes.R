library(data.table)
library(stringr)

# If a filtered list of variants exists, skip the filtering
hasFilteredList = file.exists("filtered_variants.csv")

if (!hasFilteredList) {
  # Load ClinVar variants and submissions
  variants = fread("variant_summary.txt")
  submissions = fread("submission_summary.txt")
  
  # Filter out variants that are
  # 1) missense variants, and
  filteredVariants = variants[str_detect(Name, "\\(p\\..*\\)") & 
                                !str_detect(Name, "\\(p\\..*(=|Ter)\\)")]
  
  # 2) are classified to have uncertain significance
  filteredVariants = filteredVariants[str_detect(ClinicalSignificance, "Uncertain significance")]
  
  # Here is a list of clinical significance categories that were removed after the filter
  # Uncertain significance
  # Uncertain significance, risk factor
  # Uncertain significance, drug response
  # Uncertain significance, other
  
  # Select submissions that are collected through clinical testing and not through curating the literature
  ids = filteredVariants[, unique(VariationID)]
  filteredSubmissions = submissions[`#VariationID` %in% ids & 
                                      str_detect(CollectionMethod, "clinical testing") & 
                                      !str_detect(CollectionMethod, "literature only"), 
                                    .(var_id = `#VariationID`, gene_symbol = SubmittedGeneSymbol, 
                                      submitter = Submitter, scv = SCV)]
  
  # Here is a list of collection methods that were removed after the filter
  # clinical testing
  # clinical testing;provider interpretation
  # clinical testing;curation
  
  # Remove "-" in gene name
  filteredSubmissions = filteredSubmissions[gene_symbol != "-"]
  
  # Remove some columns to save space
  filteredVariants$PhenotypeList = NULL
  filteredVariants$PhenotypeIDS = NULL
  filteredVariants$RCVaccession = NULL
  filteredVariants$OtherIDs = NULL
  
  # Save the filtered list
  fwrite(filteredVariants, "filtered_variants.csv")
} else {
  filteredVariants = fread("filtered_variants.csv")
}

# Rank genes by their number of unique VUS variants
rankedGenes = unique(filteredVariants[HGNC_ID != "-" & !str_detect(GeneSymbol, "subset|cover")])
rankedGenes = rankedGenes[, .(num_variants = .N), by = c("GeneSymbol", "HGNC_ID")]
rankedGenes = rankedGenes[order(num_variants, decreasing = T)]
rankedGenes = rankedGenes[, .(symbol = GeneSymbol, hgnc_id = HGNC_ID, num_variants)]

# Save the ranked
fwrite(rankedGenes, "genes_ranked_filtered_variants.csv")