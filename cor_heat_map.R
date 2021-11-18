# Set a value seed for reproucability -----------------------------------------

# This ensures that the correlation calculations give the same result each time
# the script is run by ensuring that the same number is used as the seed for
# Rs pseudo-random number generator.
set.seed(8)

# Load the required R packages ------------------------------------------------
library(tidyverse) # this is required to reshape the data and make the plot
library(caret) # this is required to calculate the correlation matrix
library(RColorBrewer) # this is required to set the colours in the heatmap


# Load and prepare the data ---------------------------------------------------

# Read in normalised intensity values 
normalised <- read.csv("output_tables/log2.lfq.filtered.csv")

# Set the row names of the data frame to equal the "X" column which is the 
# protein IDs.
colnames(normalised)
row.names(normalised) <- normalised$X

# You can now get rid of the "X" column leaving only the samples
normalised <- normalised[,-1]

# Check that this has worked
colnames(normalised)
rownames(normalised[1:10,])

# Calculate correlation matrix between samples --------------------------------

# Store the length of the columns (the number of samples) as n
n=length(colnames(normalised))

# Next we will rename some of the columns using gsub for a neater plot.
# 
# Here we are removing "LFQ.intensity" from the start of each sample.
#
# Then we are replacing "B5" and "B6" with "B4" and "B5" respectively to
# account for a missed number in the bioreps.
# 
# Then we are converting the sample names to mm by replacing 
# "0" with "0." and so on.
#
# Lastly we are replacing the underscore("_") with a space
#
# Look at the column (sample) names
colnames(normalised)

colnames(normalised) <- gsub('LFQ.intensity.', # replace ""
                             '', # with ""
                             colnames(normalised)) # in
colnames(normalised) <- gsub('B5',
                             'B4',
                             colnames(normalised))
colnames(normalised) <- gsub('B6',
                             'B5',
                             colnames(normalised))
colnames(normalised) <- gsub('_0',
                             '_0.',
                             colnames(normalised))
colnames(normalised) <- gsub('_10',
                             '_1.0',
                             colnames(normalised))
colnames(normalised) <- gsub('_11',
                             '_1.1',
                             colnames(normalised))
colnames(normalised) <- gsub('_',
                             ' ',
                             colnames(normalised))

# Look at the changed sample (column) names
colnames(normalised)

# Then we calculate pairwise correlation between all protein lFQ intensities (rows)
# in all samples (columns)
correlationMatrix <- cor(normalised[1:n],
                         use="pairwise.complete.obs")

# We then convert this result into a data frame
cor_df <- as.data.frame(correlationMatrix)

# Then we can add the sample names in a new column
cor_df$Sample <- rownames(cor_df)

# Check the data structure, you can see the "Sample" column has been added
# at the end
colnames(cor_df)

# Convert data from wide format to long format using pivot_longer
long.data <- pivot_longer(cor_df, # input data frame to reshape
                    1:30, # which columns to gather 
                    names_to = "Sample_B", # name for the new sample column 
                    values_to = "correlation") # name for the correlation value column

# Look at the first ten rows
long.data[1:10,]

# Next we'll add the units to the sample names by inserting a new column with the unit
# then joining the name and unit columns
long.data$unit <- "mm"
long.data <- unite(long.data,
                        Sample,
                        unit,
                        col=Anther_A,
                        sep = " ",
                        remove = FALSE)

# Do the same with the second sample column
long.data <- unite(long.data,
                   Sample_B,
                   unit,
                   col=Anther_B,
                   sep = " ",
                   remove = FALSE)

# Plot the heatmap ------------------------------------------------------------
heat <- ggplot(long.data, # use the long.data data frame
               aes(Anther_A, # plot anther A
                   Anther_B)) + # vs anther B
        geom_tile(aes(fill = correlation), # plot tiles whose fill is the pairwise correlation value
                  colour = "white") + 
        scale_fill_distiller(palette = "RdYlBu") + # use the "RdYlBu" pallette for the tile fills
        theme(axis.text.x = element_text(angle = 90, # turn x axis labels 90 degrees
                                         size = 10, # make them 10 pt
                                         vjust = 0.5), # centre them on the tiles
        axis.text.y = element_text(size = 10), # make the y axis labels 10 pt too
        axis.title.x = element_blank(), # no titles in the x
        axis.title.y = element_blank(), # or y axis
        legend.position = "bottom") + # place the legend below the plot
        labs(fill = "Correlation coefficient") # set the label for the legend
heat

long.data$Anther_A <- as.factor(long.data$Anther_A)
levels(long.data$Anther_A)

long.data$Anther_A <- factor(long.data$Anther_A,
                     levels(long.data$Anther_A)[c(1,7,13,19,25,
                                                  2,8,14,20,26,
                                                  3,9,15,21,27,
                                                  4,10,16,22,28,
                                                  5,11,17,23,29,
                                                  6,12,18,24,30)])

long.data$Anther_B <- as.factor(long.data$Anther_B)
levels(long.data$Anther_B)

long.data$Anther_B <- factor(long.data$Anther_B,
                             levels(long.data$Anther_B)[c(1,7,13,19,25,
                                                          2,8,14,20,26,
                                                          3,9,15,21,27,
                                                          4,10,16,22,28,
                                                          5,11,17,23,29,
                                                          6,12,18,24,30)])

heat <- ggplot(long.data, # use the long.data data frame
               aes(Anther_A, # plot anther A
                   Anther_B)) + # vs anther B
  geom_tile(aes(fill = correlation), # plot tiles whose fill is the pairwise correlation value
            colour = "white") + 
  scale_fill_distiller(palette = "RdYlBu") + # use the "RdYlBu" pallette for the tile fills
  theme(axis.text.x = element_text(angle = 90, # turn x axis labels 90 degrees
                                   size = 10, # make them 10 pt
                                   vjust = 0.5), # centre them on the tiles
        axis.text.y = element_text(size = 10), # make the y axis labels 10 pt too
        axis.title.x = element_blank(), # no titles in the x
        axis.title.y = element_blank(), # or y axis
        legend.position = "bottom") + # place the legend below the plot
  labs(fill = "Correlation coefficient") # set the label for the legend
heat


 # Save the plot as an svg at 180 mm wide
 ggsave("Figure2C_3.svg",
        width=180,
        units = "mm")
 
# # Or as a TIFF at 600 dpi and 180 mm wide
# ggsave("Figure2C.tiff",
#        width=180,
#        units = "mm",
#        dpi=600)
