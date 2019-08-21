library(data.table)
library(stringr)
library(magrittr)
library(ggplot2)
library(scales)
library(Cairo)

# Load ClinVar variants and submissions
variants = fread("variant_summary.txt")
submissions = fread("submission_summary.txt")

# Filter out variants that are
# 1) single-nucleotide, missense variants, and
filteredVariants = variants[Type == "single nucleotide variant" & 
                              str_detect(Name, "\\(p\\..*\\)") & 
                              !str_detect(Name, "\\(p\\..*(=|Ter)\\)")]

# 2) are classified to have uncertain significance
filteredVariants = filteredVariants[str_detect(ClinicalSignificance, "Uncertain significance")]

# Here is a list of clinical significance categories that were left after the filter
# Uncertain significance
# Uncertain significance, risk factor
# Uncertain significance, drug response
# Uncertain significance, other

# Filter out submissions that are collected through clinical testing and not through literature
ids = filteredVariants[, unique(VariationID)]
filteredSubmissions = submissions[`#VariationID` %in% ids & 
                                    str_detect(CollectionMethod, "clinical testing") & 
                                    !str_detect(CollectionMethod, "literature only"), 
                                  .(var_id = `#VariationID`, gene_symbol = SubmittedGeneSymbol, 
                                    submitter = Submitter, scv = SCV)]

# Here is a list of collection methods that were left after the filter
# clinical testing
# clinical testing;provider interpretation
# clinical testing;curation

# Remove some columns to save space
filteredVariants$PhenotypeList = NULL
filteredVariants$PhenotypeIDS = NULL
filteredVariants$RCVaccession = NULL
filteredVariants$OtherIDs = NULL

# Save the filtered variants and submissions
fwrite(filteredVariants, "filtered_variants.csv")
fwrite(filteredSubmissions, "filtered_submissions.csv")

# Plot the distribution of # variants deposited by each submitter
sumVariantsPerSubmitter = filteredSubmissions[, .(num_variants = .N), by = "submitter"]
plot = ggplot(sumVariantsPerSubmitter) +
  geom_histogram(aes(x = num_variants), bins = 100) +
  scale_x_log10(labels=comma) + labs(x = "# Variants submitted to ClinVar", y = "# Submitters") +
  theme_minimal(base_size = 20) + theme(plot.margin = margin(0.1, 0.4, 0.1, 0.1, unit = "in"))
ggsave("num_variants_submitted_histogram.png", plot, type = "cairo-png", 
       width = 8, height = 6, dpi = 300, units = "in")
sumVariantsPerSubmitter[, percent_variants := num_variants / sum(num_variants) * 100]
sumVariantsPerSubmitter = sumVariantsPerSubmitter[order(num_variants, decreasing = T)]
