# First lets get a look at how comparable the samples are. Do the same proteins
# show up in every sample? Across biological replicates? 

# Load packages we'll use.
library(tidyr) # for making a long version of the dataframe with the pivot_longer() function.
library(dplyr) # for the distinct() function to remove dulicates
library(UpSetR) # for plotting the intersections between samples.

# Loading in the data table.
prot.vs.lfq <- read.csv("output_tables/protein_vs_lfq.csv")

# Set the sample names as the row names.
row.names(prot.vs.lfq) <- prot.vs.lfq$sample
colnames(prot.vs.lfq[,1:5])

# Counting zeros --------------------------------------------------------------
# What proportion of total measurements are zero? 

# First we'll use tidyr's gather to create a long version of the data with a row 
# for each protein in each sample and it's LFQ intensity measurement.
long.data <- gather(prot.vs.lfq, # input data frame to reshape
                    key="BaMP", # name for the new key column (name of gathered columns)
                    value = "intensity", # name for the new value column (protein LFQ intensity)
                    4:6434) # which columns to gather (all protein LFQ intensity columns)
long.data[1:10,]

# Subset the zero values i.e. every measurement for a protein in a sample that
# is equal to zero using the subset() function.
zeros <- subset(long.data, # data frame to subset
                intensity==0) # column name == value to subset.

# Calculate the proportion of zeros in relation to all measurements--the number
# of zeros observed divided by the number of measurements for each protein
# in each sample.
total.zero.prop <- nrow(zeros)/nrow(long.data)
# Just under 60% of measurements are zero.

# Draw a density plot to visualise the distribution. 
ggplot(long.data,
       aes(x=log10_intensity)) +
  geom_density(alpha=0.5,
               fill = "#d95f0e") +
  theme(text = element_text(size=12)) +
  xlab("Log10 LFQ intensity") +
  ylab("Density")

ggsave("Fig2A.svg",
       width = 90,
       height = 72.355,
       units = "mm")

# Why so many zeros? ----------------------------------------------------------
# To compare samples they must have shared values. How many of the proteins in
# this dataset are in more than one sample? 
# First lets make a boolean presence/abscence table for proteins where
# present==1 and absent==0 
boolean.lfq <- prot.vs.lfq
for (x in 4:ncol(boolean.lfq)) { # for each column in columns 4 to the last column in the table
  boolean.lfq[,x] <- ifelse(boolean.lfq[,x]>0,1,0) # if the value is more than 0 change it to 1;
  # otherwise keep it as 0.
}

boolean.lfq.len <- aggregate(boolean.lfq[,4:6434],
                               by=list(Biorep=boolean.lfq$anther.len),
                               FUN=sum)
rownames(boolean.lfq.len) <- boolean.lfq.len$Biorep

for (x in 2:ncol(boolean.lfq.len)) { # for each column in columns 4 to the last column in the table
  boolean.lfq.len[,x] <- ifelse(boolean.lfq.len[,x]>0,1,0) # if the value is more than 0 change it to 1;
  # otherwise keep it as 0.
}

rownames(boolean.lfq.len)
colnames(boolean.lfq.len[,1:10])

# Use gather to give a count for each protein across all reps for a given anther
# length
long.len <- gather(boolean.lfq.len, # input data frame to reshape
                     key="BaMP", # name for the new key column (name of gathered columns)
                     value = "count", # name for the new value column (protein LFQ intensity)
                     2:6432) # which columns to gather (all protein LFQ intensity columns)
long.len[1:10,]
colnames(boolean.lfq.len[,6420:6432])

##############################################################################

# Transpose this presence/abscence table and generate a count for each protein
# with the number of times a protein is measured in the entire dataset.
t.boolean.lfq.len <- t(boolean.lfq.len[,2:6432])
t.boolean.lfq.len <- as.data.frame(t.boolean.lfq.len)


# Plot the intersects between all 30 samples (5 replicates; 6 anther lengths)
# using UpSetR.
# Each sample will be a set. Store the sample names (names for protein LFQ 
# intensity columns) as a list.
setnames <- colnames(t.boolean.lfq.len[,1:6])

lenupset <- upset(t.boolean.lfq.len, # data frame 
      order.by = "freq", # sort by largest to smallest
      nintersects = 20, # only plot the top 25 intersects
      sets = setnames, # sets to look for intersects
      keep.order = TRUE,
      text.scale = 2.6,
      set_size.show = TRUE,
      set_size.scale_max = 8000,
      mb.ratio = c(0.6, 0.4),
      point.size = 4)# Keep the order of sets in setnames

svg(file="Fig2B.svg",
    width = 8,
    height = 6.4,
    onefile = FALSE)

lenupset

dev.off()

