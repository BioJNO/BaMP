scores <- read.csv("log2.prot.scrs.csv")

ontology <- read.csv("GeneOntology.csv", header = TRUE)

colnames(ontology)
colnames(scores)

library(tidyverse)

ontology <- ontology %>% rename(GO = id)

colnames(ontology)
scores$GO <- gsub('\\.', ':', scores$GO)

merge_GOs <- merge(scores, ontology, by="GO")

significant <- subset(scores, p.adj<=0.01)

merge_sig <- merge(significant, ontology, by="GO")

setdiff(significant$GO, merge_sig$GO)
# write.csv(merge_sig, "significant_gos.csv")

library(caret)
colnames(merge_sig)
gos <- read.csv("wide.go.csv")
gos <- gos[,-c(1,2)]
colnames(gos[,4970:4974])

binary_go <- gos
for (x in 1:ncol(binary_go)) {
  binary_go[,x] <- ifelse(binary_go[,x]>0,1,0)
}
library(reshape2)
long.data <- melt(binary_go, id="Category",
                  variable.name = "GO",
                  value.name = "count")

keep <- subset(long.data, count>=3)
library(dplyr)
keep <- keep %>% distinct(BaMP)
keep <- keep$BaMP
t.boolean.lfq.filtered <- subset(t.boolean.lfq,
                                 rownames(t.boolean.lfq) %in% keep)

GO_cor <- cor(gos, method = "spearman")
highCorr <- sum(abs(GO_cor[upper.tri(GO_cor)]) > .75)
highCorr
correlations <- as.data.frame(GO_cor)
write.csv(correlations, "GO_correlations.csv")

GO <- read.csv("wide.go.csv")
merged <- merge(t.log2.lfq, GO, by="sample")


