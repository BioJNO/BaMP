library(datasets)

data(iris)

summary(iris)

colnames(iris)

# Principal component analysis ------------------------------------------------
set.seed(16)
pca <- prcomp(iris[,1:4])
library(factoextra)
library(viridis)

# plot eigenvalues
fviz_eig(pca)

# Plot PCAs
fviz_pca_ind(pca,
             label = "none",
             habillage=iris$Species) +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete = TRUE)

ggsave("iris_pca_species.svg")

# NMDS ------------------------------------------------------------------------
library(vegan)
library(ggplot2)
nmds <- metaMDS(iris[,1:4], trymax = 1000, k=2, na.rm=TRUE)
nmds
# stress=0.03775519
scrs <- as.data.frame(scores(nmds, display = "sites"))

scrs <- cbind(scrs, Species = iris$Species)

p <- ggplot(scrs) + 
  geom_point(mapping = aes(x=NMDS1,
                           y=NMDS2,
                           color = Species,
                           shape = Species,
                           size=3)) +
  scale_color_viridis(discrete = TRUE, name = "Species") + 
  scale_shape_discrete(name = "Species")

p

ggsave("iris.species.nmds.svg")


# Fit metadata to NMDS --------------------------------------------------------
meta.fit <- envfit(nmds, iris[,1:4], perm=999)
meta.fit

meta.scrs <- data.frame((meta.fit$vectors)$arrows*sqrt(meta.fit$vectors$r),
                        (meta.fit$vectors)$r,
                        (meta.fit$vectors)$pvals)
meta.scrs <- cbind(meta.scrs, species = rownames(meta.scrs))
library(plyr)
meta.scrs <- rename(meta.scrs,
                    c("X.meta.fit.vectors..r" = "r",
                      "X.meta.fit.vectors..pvals" = "p",
                      "species" = "Variables"))

meta.scrs$p.adj <- p.adjust(meta.scrs$p,
                            method = "BH",
                            n = length(meta.scrs$p))

adjsig <- subset(meta.scrs, p.adj<0.01)
colnames(adjsig)


# Plot envfit -----------------------------------------------------------------
p <- ggplot(scrs) +
  geom_point(mapping = aes(x=NMDS1,
                           y=NMDS2,
                           color = Species,
                           shape = Species,
                           size=3)) +
  scale_color_viridis(discrete = TRUE, name ="Species") + 
  scale_shape_discrete(name = "Species") +
  geom_segment(data = adjsig,
               aes(x = 0, xend=(NMDS1*0.3),
                   y = 0, yend=(NMDS2*0.3)),
               arrow = arrow(length = unit(0.25, "cm")),
               colour = "grey") +
  geom_text(data = adjsig,
            aes(x=(NMDS1*0.3),
                y = (NMDS2*0.3),
                label = Variables,
                angle = (180/pi) * atan((NMDS2)/(NMDS1))),
            size = 3) +
  ylim(-0.25,0.35)
p

ggsave("iris.species.envfit.nmds.svg")
