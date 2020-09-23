library(data.table)
library(stringr)

###
# STEP 1: filtering ClinVar variants
###
# If a filtered list of variants exists, skip the filtering
hasFilteredList = file.exists("filtered_variants.csv") && file.exists("missense_variants.csv")

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
  fwrite(missenseVariants, "missense_variants.csv")
  fwrite(filteredVariants, "filtered_variants.csv")
} else {
  missenseVariants = fread("missense_variants.csv")
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

###
# STEP 2: Calculate movability and reappearance parameters
###
# Load Invitae variant count
invitae = fread("invitae_variant_count.csv")
invitae = invitae[, .(gene, missense_vus_unique = unique_vus, 
                      missense_vus_movable_unique = movable_vus, 
                      missense_vus_observed = capped_occurance_vus,
                      missense_vus_movable_observed = observed_movable_vus)]

# Load ClinVar variant count
clinvar = merged

# Filter ClinVar genes to remove the ones without any missesn VUS variants
# and format the rest
clinvar = clinvar[num_unique_missense_vus > 0, .(gene = symbol, 
                                                 clinvar_missense_vus = num_unique_missense_vus, 
                                                 hgnc_id = as.integer(substr(hgnc_id, 6, nchar(hgnc_id))))]

# Match Invitae and Clinvar variant counts
merged = merge(clinvar, invitae, by = "gene", all.x = T)
merged[, from_invitae := complete.cases(merged)]

# Derive parameters from Invitae variant dataset
# NOTE: See the Materials and Methods section of the manuscript for detail

# Calculate movability fraction (See Manuscript Section 2.4)
merged[, missense_vus_movability_fraction := missense_vus_movable_unique / missense_vus_unique]
avg_missense_vus_movability_fraction = mean(na.omit(merged$missense_vus_movability_fraction))

# Calculate reappearance
merged[, missense_vus_patients_per_variant := missense_vus_observed / missense_vus_unique]
avg_missense_vus_patients_per_variant = mean(na.omit(merged$missense_vus_patients_per_variant))

# Replace NA with 0 for downstream calculations
for (j in names(merged)) set(merged,which(is.na(merged[[j]])),j,0)

###
# STEP 3: Apply the movability and reappearance parameters to ClinVar genes
###
# Calculate Regularized Movability Fraction (See Manuscript Section 2.4)
pseudo = 8
merged[, weighted_movability_fraction := 
         (missense_vus_unique * missense_vus_movability_fraction + 
            pseudo * avg_missense_vus_movability_fraction) / 
         (missense_vus_unique + pseudo)]

# Calcuate Weighted Patients per Variant (WPV)
merged[, weighted_patients_per_variant := 
         (missense_vus_unique * missense_vus_patients_per_variant + 
            pseudo * avg_missense_vus_patients_per_variant) / 
         (missense_vus_unique + pseudo)]

# Calculate Patient Weighted Movable VUS Count (PWMVC)
merged[, patient_weighted_movable_vus_count := clinvar_missense_vus * 
         weighted_movability_fraction * weighted_patients_per_variant]

# Get A.A. length
geneLenth = fread("gene_length.csv")
merged = merge(merged, geneLenth[, .(hgnc_id, aa_length = canonical_aa_length)], by = "hgnc_id", all.x = T)

# Calculate Difficulty-Adjusted Impact Score (DAIS)
aaLengthFixed = 300
merged[, difficulty_adjusted_impact_score := patient_weighted_movable_vus_count / (aa_length + aaLengthFixed)]

# Handle missing values in number-type columns
for (col in which(sapply(merged, class) == "integer")) {
  merged[[col]] = lapply(merged[[col]], function(e) if (is.na(e)) 0 else e)
}

###
# STEP 4: Rank genes
###
# Rank by three strategies (unique ClinVar missense variants, MARWIS, DAIS)
ranked = merged
ranked$clinvar_missense_vus = unlist(ranked$clinvar_missense_vus)
ranked[, rank_by_clinvar_vus := frank(ranked, -clinvar_missense_vus, ties.method = "random")]
ranked[, rank_by_marwis := frank(ranked, -patient_weighted_movable_vus_count, ties.method = "random")]
ranked[, rank_dais := frank(ranked, -difficulty_adjusted_impact_score, ties.method = "random")]

# Save to file
fwrite(ranked[order(rank_dais)], "ranked_clinvar_genes.csv")
