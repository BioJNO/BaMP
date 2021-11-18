# Ordination of samples based on protein LFQ intensity data 
prot.vs.lfq <- read.csv("protein_vs_lfq.csv")
colnames(prot.vs.lfq[,1:5])
row.names(prot.vs.lfq) <- prot.vs.lfq$sample

set.seed(16)

# Principal component analysis raw lfq ----------------------------------------
pca <- prcomp(prot.vs.lfq[,4:6434])
library(factoextra)
library(viridis)

# plot eigenvalues
fviz_eig(pca)

# Plot PCAs
fviz_pca_ind(pca,
     label = "none",
     habillage=prot.vs.lfq$batch) +
     scale_color_viridis(discrete=TRUE) +
     scale_fill_viridis(discrete = TRUE)

ggsave("raw_intensity_pca_batch.svg", width=183, units="mm")

fviz_pca_ind(pca,
             label = "all",
             habillage=prot.vs.lfq$anther.len) +
             scale_color_viridis(discrete = TRUE) +
             scale_fill_viridis(discrete = TRUE)

ggsave("raw_intensity_pca_length.svg", width=183, units="mm")


# Principal component analysis filtered lfq ----------------------------------------
filtered.prot.vs.lfq <- read.csv("prot.vs.lfq.filtered.csv")
colnames(filtered.prot.vs.lfq[,1:10])
row.names(filtered.prot.vs.lfq) <- filtered.prot.vs.lfq$sample
pca <- prcomp(filtered.prot.vs.lfq[,4:3123])
library(factoextra)
library(viridis)

# plot eigenvalues
fviz_eig(pca)

# Plot PCAs
fviz_pca_ind(pca,
             label = "all",
             habillage=filtered.prot.vs.lfq$batch) +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete = TRUE)

ggsave("filtered_intensity_pca_batch.svg", width=183, units="mm")

fviz_pca_ind(pca,
             label = "all",
             habillage=filtered.prot.vs.lfq$anther.len) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE)

ggsave("filtered_intensity_pca_length.svg", width=183, units="mm")

# Raw intensity NMDS ----------------------------------------------------------
library(vegan)
library(ggplot2)
nmds <- metaMDS(prot.vs.lfq[,4:6434], trymax = 1000, k=2, na.rm=TRUE)
nmds
# stress=0.08997857
scrs <- as.data.frame(scores(nmds, display = "sites"))

scrs <- cbind(scrs, anther.len = prot.vs.lfq$anther.len)
scrs <- cbind(scrs, batch = prot.vs.lfq$batch)
scrs <- cbind(scrs, sample = prot.vs.lfq$sample)

scrs$anther.len <- as.factor(scrs$anther.len)

p <- ggplot(scrs) + 
  geom_point(mapping = aes(x=NMDS1,
                           y=NMDS2,
                           color = batch,
                           shape = batch,
                           size=3)) +
  scale_color_viridis(discrete = TRUE, name = "Batch") + 
  scale_shape_discrete(name = "Batch")

p

ggsave("raw.intensity.batch.nmds.svg", dpi=600, width=183, units = "mm")

# Filtered intensity NMDS -----------------------------------------------------
library(vegan)
library(ggplot2)
nmds <- metaMDS(filtered.prot.vs.lfq[,4:3123], trymax = 1000, k=2, na.rm=TRUE)
nmds
# stress=0.1064837
scrs <- as.data.frame(scores(nmds, display = "sites"))

scrs <- cbind(scrs, anther.len = prot.vs.lfq$anther.len)
scrs <- cbind(scrs, batch = prot.vs.lfq$batch)
scrs <- cbind(scrs, sample = prot.vs.lfq$sample)

scrs$anther.len <- as.factor(scrs$anther.len)

p <- ggplot(scrs) + 
  geom_point(mapping = aes(x=NMDS1,
                           y=NMDS2,
                           color = batch,
                           shape = batch,
                           size=3)) +
  scale_color_viridis(discrete = TRUE, name = "Batch") + 
  scale_shape_discrete(name = "Batch")

p

ggsave("filtered.intensity.batch.nmds.svg", dpi=600, width=183, units = "mm")




