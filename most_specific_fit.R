wide_GO_spec <- pivot_wider(most_specific, # input data frame
                            names_from = GO, # new column names
                            values_from = log2.LFQ) # values for new columns
colnames(wide_GO_spec[,1:10])

hist(most_specific$log2.LFQ,
     breaks = 100,
     col = "#d95f0e")

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
correlated_var <- colnames(dfr[,highlyCorrelated])

correlated_var <- data.frame("GO" = correlated_var)



correlated_var <- merge(correlated_var, geneontology, by="GO")

write.csv(correlated_var, "correlated_go_terms.csv")

wide_GO_sum_important <- dfr[, colnames(dfr) %in% important_var]
rownames(wide_GO_sum_important)

write.csv(wide_GO_sum_important,
          "output_tables/sample_vs_GO_summed_log2_intensity_most_specific.csv")

# # Impute missing values -------------------------------------------------------
# library(caret)
# preProc_filtered <- preProcess(filtered, method = c("knnImpute"))
# preProc_normalised <- preProcess(normalised, method = c("knnImpute"))
# 
# library(RANN)
# filtered_imp <- predict(preProc_filtered, filtered)
# normal_imp <- predict(preProc_normalised, normalised)

# # Determine correlation between protein intensities ---------------------------
# filter_correlationMatrix <- cor(filtered, use = "pairwise.complete.obs")
# filter_highcor <- findCorrelation(filter_correlationMatrix, cutoff=0.75)
# # 28 highly correlated proteins before normalisation
# 
# normal_correlationMatrix <- cor(normal_imp, use = "pairwise.complete.obs")
# normal_highcor <- findCorrelation(normal_correlationMatrix, cutoff=0.75)
# # 29 highly correlated proteins after normalisation

library(caret)
tnorm <- t(normalised)
tnorm <- as.data.frame(tnorm)

tnorm[is.na(tnorm)] <- 0

zv <- apply(tnorm, 2, function(x) length(unique(x)) == 1)

dfr <- tnorm[, !zv]

n=length(colnames(dfr))

correlationMatrix <- cor(dfr[,1:n],use="pairwise.complete.obs")

print(correlationMatrix)

highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=(0.8),verbose = FALSE)

print(highlyCorrelated)

important_var=colnames(tnorm[,-highlyCorrelated])

tnorm_important <- tnorm[, colnames(tnorm) %in% important_var]

# PCA with imputation ---------------------------------------------------------
pca <- prcomp(tnorm_important)
library(factoextra)
library(viridis)

# plot eigenvalues
fviz_eig(pca)
library(stringr)
tnorm_important$batch <- str_split_fixed(rownames(tnorm_important), "_", 2)[,1]
tnorm_important$anther.len <- str_split_fixed(rownames(tnorm_important), "_", 2)[,2]

# Plot PCAs
fviz_pca_ind(pca,
             label = "none",
             habillage=tnorm_important$batch) +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete = TRUE)

fviz_pca_ind(pca,
             label = "all",
             habillage=tnorm_important$anther.len) +
  scale_color_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE)




