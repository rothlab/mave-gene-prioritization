library(data.table)
library(ggplot2)
library(Cairo)
library(scales)
library(magrittr)
library(viridis)
library(gplots)
library(dendsort)

# Load filtered submissions and stats
submissions = fread("filtered_submissions.csv")
submissionSummary = fread("num_variants_submitted.csv")

# Select submissions by the top 10 submitters
topSubmitters = submissionSummary[order(num_variants, decreasing = T), submitter][1:10]
topSubmissions = submissions[submitter %in% topSubmitters]

# Format submitter names
shortName = c("Illumina Clinical Services Laboratory,Illumina" = "Illumina",
              "Laboratory for Molecular Medicine,Partners HealthCare Personalized Medicine" = "Partners HealthCare",
              "Genetic Services Laboratory, University of Chicago" = "UC Genetic Services",
              "Athena Diagnostics Inc" = "Athena Diagnostics",
              "EGL Genetic Diagnostics,Eurofins Clinical Diagnostics" = "Eurofins")
topSubmitters[topSubmitters %in% names(shortName)] = shortName[topSubmitters[topSubmitters %in% names(shortName)]]
topSubmitters = sort(topSubmitters, decreasing = T)
topSubmissions[submitter %in% names(shortName), submitter := shortName[submitter]]

# Summarize by genes and submitters
variantsOnGene = topSubmissions[, .(num_variants = .N), by = c("gene_symbol", "submitter")]
combs = t(expand.grid(1:length(topSubmitters), 1:length(topSubmitters))) # This creates a square
# combs = do.call("cbind", sapply(10:1, function(e) matrix(c(seq(1, e), rep(e, e)), nrow = 2, byrow = T))) # This creates a diagnoal
corr = apply(combs, 2, function(comb) {
  # Get submitters
  submitterX = topSubmitters[comb[[1]]]
  submitterY = topSubmitters[comb[[2]]]
  
  # Find common genes
  common = merge(variantsOnGene[submitter == submitterX, .(gene_symbol, x = num_variants)], 
                 variantsOnGene[submitter == submitterY, .(gene_symbol, y = num_variants)],
                 by = "gene_symbol")
  
  # Calculate pearson correlation
  corr = cor(common$x, common$y)
  
  return(data.table(x = submitterX, y = submitterY, corr = corr, num_common_genes = nrow(common)))
})
corr = rbindlist(corr)

# Generate heatmap with clustering
corrMatrix = matrix(corr$corr, nrow = length(topSubmitters), ncol = length(topSubmitters), 
                    dimnames = list(topSubmitters, topSubmitters))
png(file = "pearson_correlation_heatmap_clustered.png", type = "cairo", width = 8, height = 6, res = 300, units = "in")
heatmap.2(corrMatrix, margins = c(10, 12), srtCol = 45, dendrogram = "col",
          trace = "none", col = viridis(100, direction = -1),
          hclustfun=function(x) dendsort(hclust(x, method="ward.D2"), isReverse = T, type = "min"))
dev.off()

stop()

# Generate heatmap
# plot = ggplot(corr, aes(x = x, y = y, fill = corr, color = "")) + 
#   geom_tile() + theme_minimal(base_size = 16) + 
#   geom_text(aes(x = x, y = y, label = round(corr, digits = 2)), color = "black", size = 4) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank()) +
#   scale_fill_gradient2(low = viridis(5, direction = -1)[2], high = viridis(5, direction = -1)[4], 
#                        midpoint = 0, limit = c(-1, 1), name="Pearson\nCorrelation", na.value="lightyellow") +
#   labs(x = NULL, y = NULL) + scale_colour_manual(values = NA) + 
#   guides(colour=guide_legend("No data", override.aes=list(fill="lightyellow", color = "lightgrey")))
# ggsave("pearson_correlation_heatmap.png", plot, type = "cairo-png",
#        width = 8, height = 6, dpi = 300, units = "in")

# Only look at correlations represented by at least 20 genes
corr = corr[num_common_genes >= 20]
plot = ggplot(corr) + geom_histogram(aes(x = corr, y =..density..), bins = 40) +
  geom_density(aes(x = corr, y =..density..)) + xlim(-0.5, 1.5) +
  theme_minimal(base_size = 20) + labs(x = "Pearson Correlation", y = "Density")
ggsave("pearson_correlation_density.png", plot, type = "cairo-png",
       width = 8, height = 6, dpi = 300, units = "in")

# Plot distribution of submitters
lowPearson = as.data.table(t(table(corr[corr <= 0.4, c(x, y)]))) %>% 
  .[, .(submitter = V2, count = N, pearson = "<= 0.4")] %>% .[, percent := count / sum(count)]
highPearson = as.data.table(t(table(corr[corr > 0.4, c(x, y)]))) %>% 
  .[, .(submitter = V2, count = N, pearson = "> 0.4")] %>% .[, percent := count / sum(count)]
plot = ggplot(rbind(lowPearson, highPearson), 
              aes(x = pearson, y = percent, fill = submitter, label = percent(percent))) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  geom_text(size = 5, position = position_stack(vjust = 0.5)) +
  theme_minimal(base_size = 20) + scale_y_continuous(labels = percent) +
  scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", 
                               "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")) +
  labs(x = "Pearson Correlation", y = "% of submitters", fill = "Submitter")
ggsave("submitters_enrichment_in_pearson_correlation.png", plot, type = "cairo-png",
       width = 12, height = 12, dpi = 300, units = "in")
