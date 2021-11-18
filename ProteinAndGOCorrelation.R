
set.seed(8)

# filtered <- read.csv("t.prot.vs.lfq.filtered.csv")
# row.names(filtered) <- filtered$X
# filtered <- filtered[,-1]
# colnames(filtered)
# rownames(filtered[1:10,])

normalised <- read.csv("output_tables/log2.lfq.filtered.csv")
colnames(normalised)
row.names(normalised) <- normalised$X
normalised <- normalised[,-1]
colnames(normalised)
rownames(normalised[1:10,])
normalised$BaMP <- rownames(normalised)

# Transform data for GO annotations
annotated <- read.csv("output_tables/Proteomic_data.csv")
colnames(annotated)

merged <- merge(normalised, annotated, by="BaMP")
colnames(merged)
merged <- merged[,-c(32:33,35:77)]
colnames(merged)

# split the Panzzer GO annotation to give a row for each GO.
library(tidyverse)
df = separate_rows(merged,32,sep = ",")
colnames(df)

# Create a table of normalised GO intensity sums per sample
df <- df[,-1]
df <- df[,c(31,1:30)]
colnames(df)
df.long <- gather(df,
                  key="sample",
                  value = "log2.LFQ",
                  2:31)
colnames(df.long)

summed <- aggregate(.~PANNZER_GO+sample, df.long, sum)
summed[1:10,]$PANNZER_GO
summed <- summed[-1,]
summed[1:10,]$PANNZER_GO
colnames(summed)

df.long.go <- df.long[,-2]
summed_total <- aggregate(.~PANNZER_GO, df.long.go, sum)

hist(summed_total$log2.LFQ,
     breaks = 100,
     col = "#d95f0e")

upstream <- read.csv("input_tables/GO_to_upstream.csv")
geneontology <- read.csv("input_tables/GeneOntology.csv")
geneontology <- rename(geneontology, GO = id)

summed <- rename(summed, GO = PANNZER_GO)
summed_total <- rename(summed_total, GO = PANNZER_GO)

sum_to_up <- merge(summed, upstream, by="GO")
sum_to_up_tot <- merge(summed_total, upstream, by="GO")

upstream <- separate_rows(sum_to_up, Upstream, sep = ", ")
upstream <- upstream$Upstream
upstream <- unique(upstream)
most_specific <- summed[!(summed$GO %in% upstream), ]
most_specific_total <- summed_total[!(summed_total$GO %in% upstream), ]

most_specific_total <- merge(most_specific_total, geneontology, by="GO")
write.csv(most_specific_total, "summed_GOs_none_upstream.csv")

hist(most_specific_total$log2.LFQ,
     breaks = 100,
     col = "#d95f0e")

# convert from long to wide format
# wide_GO_sum <- pivot_wider(summed, # input data frame
#                            names_from = PANNZER_GO, # new column names
#                            values_from = log2.LFQ) # values for new columns
# colnames(wide_GO_sum[,1:10])

wide_GO_spec <- pivot_wider(most_specific, # input data frame
                           names_from = GO, # new column names
                           values_from = log2.LFQ) # values for new columns
colnames(wide_GO_spec[,1:10])

library(caret)
wide_GO_spec[is.na(wide_GO_spec)] <- 0

zv <- apply(wide_GO_spec, 2, function(x) length(unique(x)) == 1)

dfr <- wide_GO_spec[, !zv]

colnames(dfr[,1:10])
dfr <- as.data.frame(dfr)
rownames(dfr) <- dfr$sample
dfr <- dfr[,-1]
rownames(dfr)

n=length(colnames(dfr))

correlationMatrix <- cor(dfr[1:n],
                         use="pairwise.complete.obs")

highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=(0.8),verbose = FALSE)

important_var <- colnames(dfr[,-highlyCorrelated])
not_important_var <- colnames(dfr[,highlyCorrelated])

wide_GO_sum_important <- dfr[, colnames(dfr) %in% important_var]
wide_GO_sum_unimportant <- dfr[, colnames(dfr) %in% not_important_var]
rownames(wide_GO_sum_important)

important_go <- data.frame("GO" = important_var)
important_go <- merge(important_go, geneontology, by="GO")
write.csv(important_go, "output_tables/independant_go_terms.csv")

unimportant_go <- data.frame("GO" = not_important_var)
unimportant_go <- merge(unimportant_go, geneontology, by="GO")
write.csv(unimportant_go, "output_tables/dependant_go_terms.csv")

write.csv(wide_GO_spec, 
          "output_tables/sample_vs_GO_summed_log2_intensity_most_specific.csv")
write.csv(wide_GO_sum_important,
          "output_tables/sample_vs_GO_summed_log2_intensity_important.csv")
write.csv(wide_GO_sum_unimportant,
          "output_tables/sample_vs_GO_summed_log2_intensity_unimportant.csv")
