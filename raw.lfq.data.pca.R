# Having removed the unique proteins there likely remains quantitative bias
# within the batches.
lfq.filtered <- read.csv("prot.vs.lfq.filtered.csv")
rownames(lfq.filtered) <- lfq.filtered$sample

colnames(lfq.filtered[,3120:3123])

set.seed(16)
pca <- prcomp(lfq.filtered[,4:3123])

library(viridis)
library(factoextra)
fviz_pca_ind(pca,
             label = "all",
             habillage=lfq.filtered$anther.len) +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete = TRUE)

ggsave("raw.batch.pca.svg", dpi=600, width=183, units = "mm")
