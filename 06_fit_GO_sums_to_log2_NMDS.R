library(reshape2)
# Having removed the unique proteins there likely remains quantitative bias
# within the batches.
lfq.filtered <- read.csv("prot.vs.lfq.filtered.csv", na.strings = 0)
rownames(lfq.filtered) <- lfq.filtered$sample

# Plot the distribution.
long.data <- melt(lfq.filtered, id=c("sample", "batch", "anther.len"),
                  variable.name = "BaMP", value.name = "LFQ")
hist(long.data$LFQ, breaks = 100, col = "#d95f0e")
# Still a very heavy left skew. A lot of NAs, how many?
zero <- long.data[is.na(long.data$LFQ), ]
nrow(zero)/nrow(long.data)
# 23058; 24.6%

# It's obvious there was some batch bias on the runs. Lets try a simple log2
# normalisation
t.lfq.filtered <- t(lfq.filtered[,4:3123])
log2.lfq <- log2(t.lfq.filtered)
log2.lfq <- as.data.frame(log2.lfq)
rownames(log2.lfq) <- rownames(t.lfq.filtered)
colnames(log2.lfq) <- colnames(t.lfq.filtered)
# Save the normalised output.
write.csv(log2.lfq, "log2.lfq.csv")
# Flip it again 
t.log2.lfq <- t(log2.lfq)
t.log2.lfq <- as.data.frame(t.log2.lfq)
t.log2.lfq$sample <- rownames(t.log2.lfq)
# Plot the normalised distribution.
long.data.log2 <- melt(t.log2.lfq, id="sample",
                       variable.name = "BaMP",
                       value.name = "log2.LFQ")
hist(long.data.log2$log2.LFQ, breaks = 100, col = "#d95f0e")

# Principal component analysis raw lfq ----------------------------------------
set.seed(16)
omitted.missing <- na.omit(log2.lfq)
t.log2.lfq <- t(omitted.missing)
t.log2.lfq <- as.data.frame(t.log2.lfq)
t.log2.lfq$sample <- rownames(t.log2.lfq)
colnames(t.log2.lfq)
pca <- prcomp(t.log2.lfq[,1:963])
# Add batch and anther length columns back in
library(stringr)
t.log2.lfq$batch <- str_split_fixed(t.log2.lfq$sample, "_", 2)[,1]
t.log2.lfq$anther.len <- str_split_fixed(t.log2.lfq$sample, "_", 2)[,2]

# Plot eigenvalues
library(factoextra)
fviz_eig(pca)

# Plot PCAs
library(viridis)
fviz_pca_ind(pca,
             label = "all",
             habillage=t.log2.lfq$batch) +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete = TRUE)

ggsave("log2.batch.pca.svg", width=183, units = "mm")

fviz_pca_ind(pca,
             label = "all",
             habillage=t.log2.lfq$anther.len) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE)

ggsave("log2.pca.svg", width=183, units = "mm")

# NMDS ------------------------------------------------------------------------
library(vegan)
library(ggplot2)
nmds <- metaMDS(t.log2.lfq[,1:963], trymax = 100, k=2, na.rm=TRUE)
nmds
# stress = 0.1372168
scrs <- as.data.frame(scores(nmds, display = "sites"))
scrs <- cbind(scrs, anther.len = t.log2.lfq$anther.len)
scrs <- cbind(scrs, batch = t.log2.lfq$batch)
scrs <- cbind(scrs, sample = t.log2.lfq$sample)
scrs$anther.len <- as.factor(scrs$anther.len)

p <- ggplot(scrs) + 
  geom_point(mapping = aes(x=NMDS1,
                           y=NMDS2,
                           color = anther.len,
                           shape = anther.len,
                           size=3)) +
  scale_color_viridis(discrete = TRUE, name = "Anther length (mm)") + 
  scale_shape_discrete(name = "Anther length (mm)")
p

ggsave("Figures/log2.nmds.length.svg", width=183, units = "mm")

p <- ggplot(scrs) + 
  geom_point(mapping = aes(x=NMDS1,
                           y=NMDS2,
                           color = batch,
                           shape = batch,
                           size=3)) +
  scale_color_viridis(discrete = TRUE, name = "Batch") + 
  scale_shape_discrete(name = "Batch")
p

ggsave("Figures/log2.nmds.batch.svg", width=183, units = "mm")

# Fit GO terms to NMDS
GO <- read.csv("wide.go.csv")
merged <- merge(t.log2.lfq, GO, by="sample")
colnames(merged[,960:970])

prot.fit <- envfit(nmds, merged[,968:5941], perm=999)
prot.fit

prot.scrs <- data.frame((prot.fit$vectors)$arrows*sqrt(prot.fit$vectors$r),
                        (prot.fit$vectors)$r,
                        (prot.fit$vectors)$pvals)
prot.scrs <- cbind(prot.scrs, species = rownames(prot.scrs))
library(plyr)
prot.scrs <- rename(prot.scrs,
                    c("X.prot.fit.vectors..r" = "r",
                      "X.prot.fit.vectors..pvals" = "p",
                      "species" = "GO"))
prot.scrs$p.adj <- p.adjust(prot.scrs$p,
                            method = "BH",
                            n = length(prot.scrs$p))

write.csv(prot.scrs, "log2.prot.scrs.csv")

adjsig <- subset(prot.scrs, p.adj<0.01)

write.csv(adjsig, "log2.go.scrs.pless0.01.rgreater0.75.csv")

GeneOntology <- read.csv("GeneOntology.csv")
colnames(GeneOntology)
colnames(adjsig)
GeneOntology <- rename(GeneOntology,
                       c("id" = "GO"))

adjsig$GO <- gsub('\\.', ':', adjsig$GO)

merge_GOs <- merge(adjsig, GeneOntology, by="GO")

write.csv(merge_GOs, "significant_biological_processes.csv")

missing <- setdiff(adjsig$GO, merge_GOs$GO)
# only biological processes

p <- ggplot(scrs) +
     expand_limits(x = 0.015) +
     geom_segment(data = adjsig,
                         aes(x = 0, xend=(NMDS1*0.01),
                             y = 0, yend=(NMDS2*0.01)),
                             arrow = arrow(length = unit(0.25, "cm")),
                             colour = "grey") +
     geom_text(data = adjsig,
               aes(x=(NMDS1*0.01),
                   y = (NMDS2*0.01),
                   label = BaMP,
                   angle = (180/pi) * atan((NMDS2*0.01)/(NMDS1*0.01))),
               size = 3) +
  geom_point(mapping = aes(x=NMDS1,
                           y=NMDS2,
                           color = anther.len,
                           shape = anther.len,
                           size=3)) +
  scale_color_viridis(discrete = TRUE, name ="Anther length (mm)") + 
  scale_shape_discrete(name = "Anther length (mm)")
p

# ggsave("log2.nmds.GO.fit.svg", dpi=600, width=183, units = "mm")


go.nmds <- metaMDS(merged[,968:5941], trymax = 100, k=2)
scrs <- as.data.frame(scores(go.nmds, display = "sites"))

scrs <- cbind(scrs, anther.len = t.log2.lfq$anther.len)
scrs <- cbind(scrs, batch = t.log2.lfq$batch)
scrs <- cbind(scrs, sample = t.log2.lfq$sample)

p <- ggplot(scrs) + 
  geom_point(mapping = aes(x=NMDS1,
                           y=NMDS2,
                           color = anther.len,
                           shape = anther.len,
                           size=3)) +
  scale_color_viridis(discrete = TRUE, name = "Anther length (mm)") + 
  scale_shape_discrete(name = "Anther length (mm)")
  
p


ggsave("nmds.GO.LFQ.sum.svg", dpi=600, width=183, units = "mm")
