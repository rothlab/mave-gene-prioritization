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
  missenseVariants = variants[str_detect(Name, "\\(p\\..*\\)") & 
                                !str_detect(Name, "\\(p\\..*(=|Ter)\\)")]
  
  # 2) are classified to have uncertain significance
  filteredVariants = missenseVariants[str_detect(ClinicalSignificance, "Uncertain significance")]
  
  # Here is a list of clinical significance categories that were removed after the filter
  # Uncertain significance
  # Uncertain significance, risk factor
  # Uncertain significance, drug response
  # Uncertain significance, other
  
  # Select variants that are collected through clinical testing and not through curating the literature
  selectClinicalVariants = function(sourceVariants) {
    ids = sourceVariants[, unique(VariationID)]
    filteredSubmissions = submissions[`#VariationID` %in% ids & 
                                        str_detect(CollectionMethod, "clinical testing") & 
                                        !str_detect(CollectionMethod, "literature only"), 
                                      .(var_id = `#VariationID`, gene_symbol = SubmittedGeneSymbol, 
                                        submitter = Submitter, scv = SCV)]
    filteredSubmissions = filteredSubmissions[gene_symbol != "-"]
    ids = unique(filteredSubmissions$var_id)
    
    return(sourceVariants[VariationID %in% ids])
  }
  missenseVariants = selectClinicalVariants(missenseVariants)
  filteredVariants = selectClinicalVariants(filteredVariants)
  
  # Here is a list of collection methods that were removed after the filter
  # clinical testing
  # clinical testing;provider interpretation
  # clinical testing;curation
  
  # Save the filtered list
  fwrite(filteredVariants, "filtered_variants.csv")
} else {
  filteredVariants = fread("filtered_variants.csv")
}

# Count unique variants
countUniqueVariants = function(sourceVariants) {
  counts = unique(sourceVariants[HGNC_ID != "-" & !str_detect(GeneSymbol, "subset|cover"),
                                 .(GeneSymbol, Name, HGNC_ID)])
  counts = counts[, .(num_variants = .N), by = c("GeneSymbol", "HGNC_ID")]
  
  return(counts)
}

# Count unique missense and missense VUS variants
missenseCount = countUniqueVariants(missenseVariants)
vusCount= countUniqueVariants(filteredVariants)
merged = merge(missenseCount, vusCount, by = c("GeneSymbol", "HGNC_ID"), all = T)
for (j in names(merged)) set(merged,which(is.na(merged[[j]])),j,0)
merged = merged[, .(symbol = GeneSymbol, hgnc_id = HGNC_ID, 
                    num_unique_missense = num_variants.x,
                    num_unique_missense_vus = num_variants.y)]

# Rank by the number of unique missnese VUS variants
rankedGenes = merged[order(num_unique_missense_vus, decreasing = T)]

# Save the ranked
fwrite(rankedGenes, "genes_ranked_filtered_variants.csv")
