
library(tidyverse)

# Having removed the unique proteins there likely remains quantitative bias
# within the batches.
lfq.filtered <- read.csv("output_tables/prot.vs.lfq.filtered.csv",
                         na.strings = 0)
rownames(lfq.filtered) <- lfq.filtered$sample

levels(lfq.filtered$sample)

# Plot the distribution.
long.data <- gather(lfq.filtered, # input data frame to reshape
                    key="BaMP", # name for the new key column (name of gathered columns)
                    value = "intensity", # name for the new value column (protein LFQ intensity)
                    4:3642) # which columns to gather (all protein LFQ intensity columns)
long.data[1:10,]

hist(long.data$intensity,
     breaks = 100,
     col = "#d95f0e")
long.data$logintens <- log10(long.data$intensity)

hist(long.data$logintens,
     breaks = 50,
     col = "#d95f0e")


long.data$sample <- gsub('B5', 'B4', long.data$sample)
long.data$sample <- gsub('B6', 'B5', long.data$sample)
long.data$batch <- gsub('B5', 'B4', long.data$batch)
long.data$batch <- gsub('B6', 'B5', long.data$batch)

long.data <- long.data[!grepl("B5_06", long.data$sample), ]


(maxIn <- aggregate(intensity ~ anther.len+BaMP, data=long.data, FUN=max))


annotations <- read.csv("input_tables/GO_annotations.csv")

maxIn_annot <- merge(maxIn, annotations, by="BaMP", all.x = TRUE)

maxIn_annot <- maxIn_annot[order(maxIn_annot$anther.len, maxIn_annot$intensity, decreasing = TRUE),]

write.csv(maxIn_annot, "Len_x_BaMP_summary.csv", row.names = FALSE)

# top 20 from each 
six <- subset(maxIn_annot, anther.len==6)
seven <- subset(maxIn_annot, anther.len==7)
eight <- subset(maxIn_annot, anther.len==8)
nine <- subset(maxIn_annot, anther.len==9)
ten <- subset(maxIn_annot, anther.len==10)
eleven <- subset(maxIn_annot, anther.len==11)

write.csv(six, "point_six_mm_summary.csv", row.names = FALSE)
write.csv(seven, "point_seven_mm_summary.csv", row.names = FALSE)
write.csv(eight, "point_eight_mm_summary.csv", row.names = FALSE)
write.csv(nine, "point_nine_mm_summary.csv", row.names = FALSE)
write.csv(ten, "one_mm_summary.csv", row.names = FALSE)
write.csv(eleven, "one_point_one_mm_summary.csv", row.names = FALSE)

# Still a very heavy left skew. A lot of NAs, how many?
zero <- long.data[is.na(long.data$intensity), ]
nrow(zero)/nrow(long.data)
# 35843; 32.8%

# It's obvious there was some batch bias on the runs. Lets try a simple log2
# normalisation
t.lfq.filtered <- t(lfq.filtered[,4:3642])
# write.csv(t.lfq.filtered, "t.prot.vs.lfq.filtered.csv")
log2.lfq <- log10(t.lfq.filtered)
log2.lfq <- as.data.frame(log2.lfq)
rownames(log2.lfq) <- rownames(t.lfq.filtered)
colnames(log2.lfq) <- colnames(t.lfq.filtered)
# Save the normalised output.
write.csv(log2.lfq,
          "output_tables/log2.lfq.filtered.csv")

# Flip it again 
t.log2.lfq <- t(log2.lfq)
t.log2.lfq <- as.data.frame(t.log2.lfq)
t.log2.lfq$sample <- rownames(t.log2.lfq)
colnames(t.log2.lfq[,1:10])
# Plot the normalised distribution.
long.data.log2 <- gather(t.log2.lfq, # input data frame to reshape
                    key="BaMP", # name for the new key column (name of gathered columns)
                    value = "log2.LFQ", # name for the new value column (protein LFQ intensity)
                    1:3639) # which columns to gather (all protein LFQ intensity columns)

hist(long.data.log2$log2.LFQ,
     breaks = 100,
     col = "#d95f0e")


# Principal component analysis normalised lfq ----------------------------------------
set.seed(16)
omitted.missing <- na.omit(log2.lfq)
t.log2.lfq <- t(omitted.missing)
t.log2.lfq <- as.data.frame(t.log2.lfq)
t.log2.lfq$sample <- rownames(t.log2.lfq)
pca <- prcomp(t.log2.lfq[,1:963])
# Add batch and anther length columns back in
t.log2.lfq <- separate(t.log2.lfq, # data frame
                   sample, # column to split
                   into = c("batch", "anther.len"), # names to give new columns
                   sep="_", # string that separates values to be split
                   remove = FALSE) # don't remove the column to split

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

ggsave("Figures/log2.batch.pca.svg", dpi=600, width=183, units = "mm")

fviz_pca_ind(pca,
             label = "all",
             habillage=t.log2.lfq$anther.len) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE)

ggsave("Figures/log2.pca.svg", dpi=600, width=183, units = "mm")

# NMDS ------------------------------------------------------------------------
library(vegan)
library(ggplot2)
nmds <- metaMDS(t.log2.lfq[,1:963], trymax = 100, k=2, na.rm=TRUE)
nmds
# Sress = 0.1372168
scrs <- as.data.frame(scores(nmds, display = "sites"))

scrs <- cbind(scrs, anther.len = t.log2.lfq$anther.len)
scrs <- cbind(scrs, batch = t.log2.lfq$batch)
scrs <- cbind(scrs, sample = t.log2.lfq$sample)

#scrs$anther.len <- as.character(scrs$anther.len)
#scrs$anther.len <- as.numeric(scrs$anther.len)

p <- ggplot(scrs) + 
  geom_point(mapping = aes(x=NMDS1,
                           y=NMDS2,
                           color = anther.len,
                           shape = anther.len,
                           size=3)) +
  scale_color_viridis(discrete = TRUE, name = "Anther length (mm)") + 
  scale_shape_discrete(name = "Anther length (mm)")

p

ggsave("Figures/log2.nmds.svg", dpi=600, width=183, units = "mm")

# Fit GO terms to NMDS
GO <- read.csv("output_tables/sample_vs_GO_summed_log2_intensity_most_specific.csv")
colnames(GO[1:10])
GO <- GO[,-1]

# Some of the GO terms will not be truly independent variables i.e. they will
# have a very high correlation coefficiant because they often co-occur.
# Before Fitting GO terms to the ordination we need to remove redundant 
# GO terms (feature selection). 
library(mlbench)
library(caret)

# determine the pairwise correlation coefficient of GOs terms
correlationMatrix <- cor(GO[,2:2341])
write.csv(correlationMatrix, "output_tables/most_specific_GO_cor_matrix.csv")

# Identify GO terms that are highly correlated (R > 0.75) and return the GO in
# each pairwise correlation above 0.75 with the largest mean absolute 
# correlation i.e. which of the two correlated variables is most correlated with 
# all other GO terms. 
highlyCorrelated <- findCorrelation(correlationMatrix, # correlation matrix to search
                                    cutoff=0.75, # cutoff correlation coefficient
                                    verbose = TRUE, # talk to me about what you're doing
                                    exact = TRUE, # re-evaluate correlations at each step
                                    names = TRUE) # return names not column numbers

corrframe <- as.data.frame(correlationMatrix)
highcorrframe <- corrframe[highlyCorrelated, highlyCorrelated]
highcorrframe$GO <- rownames(highcorrframe)
highcorrframe$GO <- gsub('\\.', ':', highcorrframe$GO)
highcorrframe <- merge(highcorrframe, GeneOntology, by="GO")
write.csv(highcorrframe, "output_tables/highly_correlated_GO_terms.csv")

# Remove GO terms identified as redundant by caret from the data frame. 
GO_select <- GO[, !(colnames(GO) %in% highlyCorrelated)]

# Merge GO terms with log2 LFQ values
merged <- merge(t.log2.lfq, GO_select, by="sample")
colnames(merged[,960:970])

prot.fit <- envfit(nmds, merged[,968:1716], perm=999)
prot.fit

prot.scrs <- data.frame((prot.fit$vectors)$arrows*sqrt(prot.fit$vectors$r),
                       (prot.fit$vectors)$r,
                       (prot.fit$vectors)$pvals)
prot.scrs <- cbind(prot.scrs, Species = rownames(prot.scrs))
colnames(prot.scrs)
prot.scrs <- rename(prot.scrs,
                   c(r = X.prot.fit.vectors..r,
                     p = X.prot.fit.vectors..pvals,
                     PANNZER_GO = Species))

prot.scrs$p.adj <- p.adjust(prot.scrs$p,
                           method = "BH",
                           n = length(prot.scrs$p))

# GO <- read.csv("input_tables/GO_to_upstream.csv")
GeneOntology <- read.csv("input_tables/GeneOntology.csv")
# colnames(GO)
colnames(GeneOntology)
GeneOntology <- rename(GeneOntology, GO = id)

prot.scrs$PANNZER_GO <- gsub('\\.', ':', prot.scrs$PANNZER_GO)
prot.scrs <- rename(prot.scrs, GO = PANNZER_GO)

prot.scrs.merge <- merge(prot.scrs, GeneOntology, by="GO")
colnames(prot.scrs.merge)

write.csv(prot.scrs.merge,
          "output_tables/log2.prot.scrs.csv")

adjsig <- subset(prot.scrs.merge, p.adj<0.05)
colnames(adjsig)
# correlated_terms <- read.csv("output_tables/dependant_go_terms.csv")
# uncorrelated_terms <- read.csv("output_tables/independant_go_terms.csv")

# upstream <- separate_rows(sum_to_up, Upstream, sep = ", ")
# parents <- 
# summed <- aggregate(.~PANNZER_GO+sample, df.long, sum)

# adjsig_correlated <- merge(adjsig, correlated_terms, by="GO")
# adjsig_uncorrelated <- merge(adjsig, uncorrelated_terms, by="GO")
# write.csv(adjsig_correlated, "adjsig_correlated.csv")
# write.csv(adjsig_uncorrelated, "adjsig_uncorrelated.csv")

write.csv(adjsig, "output_tables/GO_fits.csv")
# adjeff <- subset(adjsig_uncorrelated, r > 0.4)


# sig_upstream <- separate_rows(adjsig, Upstream, sep=", ")
# sig_upstream <- sig_upstream$Upstream
# 
# none_above <- adjsig[!(adjsig$GO %in% sig_upstream), ]
print(levels(adjsig_uncorrelated$category_name.x))
biol <- subset(adjeff,
               category_name.x == "biological_process")
cell <- subset(adjeff,
               category_name.x == "cellular_component")
mol <- subset(adjeff,
              category_name.x == "molecular_function")

bio_plot <- p + geom_text(data = biol,
            aes(x=(NMDS1*0.01), y = (NMDS2*0.01), label = GO,
                angle = (180/pi) * atan(NMDS2/NMDS1)),
            size = 3) +
    geom_segment(data = biol,
             aes(x = 0, xend = (NMDS1*0.01),
                 y = 0, yend = (NMDS2*0.01)),
             arrow = arrow(length = unit(0.25, "cm")),
             colour = "grey") +
    ylim(-0.01, 0.01) +
    xlim(-0.011, 0.01)
bio_plot

cell_plot <- p + geom_text(data = cell,
                           aes(x=(NMDS1*0.01), y = (NMDS2*0.01), label = GO,
                               angle = (180/pi) * atan(NMDS2/NMDS1)),
                           size = 3) +
            geom_segment(data = cell,
               aes(x = 0, xend = (NMDS1*0.01),
                   y = 0, yend = (NMDS2*0.01)),
               arrow = arrow(length = unit(0.25, "cm")),
               colour = "grey") +
           ylim(-0.01, 0.01) +
           xlim(-0.011, 0.01)
cell_plot

mol_plot <- p + geom_text(data = mol,
                           aes(x=(NMDS1*0.01), y = (NMDS2*0.01), label = GO,
                               angle = (180/pi) * atan(NMDS2/NMDS1)),
                           size = 3) +
  geom_segment(data = mol,
               aes(x = 0, xend = (NMDS1*0.01),
                   y = 0, yend = (NMDS2*0.01)),
               arrow = arrow(length = unit(0.25, "cm")),
               colour = "grey") +
  ylim(-0.01, 0.01) +
  xlim(-0.011, 0.01)
mol_plot

library(cowplot)

plot_grid(bio_plot,
          cell_plot,
          mol_plot,
          ncol=2,
          nrow=2,
          labels="auto",
          align = "h")

ggsave("Figures/GO_NMDS_fits.svg",
       width = 400,
       height = 250,
       units = "mm")
