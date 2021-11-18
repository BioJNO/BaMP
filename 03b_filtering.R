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

# Draw a quick histogram to visualise the distribution. 
hist(long.data$intensity, # list of values (column)
     breaks = 100, # break histogram in to x bars
     col = "#d95f0e") # set a nice orange colour for the bars

long.data$log10_intensity <- log10(long.data$intensity)

# Draw a quick histogram to visualise the distribution. 
hist(long.data$log10_intensity, # list of values (column)
     breaks = 100, # break histogram in to x bars
     col = "#d95f0e", # set a nice orange colour for the bars
     main = NULL,
     xlab = "Log10 LFQ intensity",
     ylab = "Frequency",
     cex.lab=1.5,
     cex.axis=1.5,
     cex.main=1.5,
     cex.sub=1.5) 

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

# Transpose this presence/abscence table and generate a count for each protein
# with the number of times a protein is measured in the entire dataset.
t.boolean.lfq <- t(boolean.lfq[,4:6434])
t.boolean.lfq <- as.data.frame(t.boolean.lfq)
t.boolean.lfq$sample.spread <- rowSums(t.boolean.lfq)

# How many don't turn up at all?
none <- subset(t.boolean.lfq,
               sample.spread==0)
nrow(none)/nrow(t.boolean.lfq)
# 65 proteins; 0.01%.

# How many are present in only one sample?
specific <- subset(t.boolean.lfq,
                   sample.spread==1)
nrow(specific)/nrow(t.boolean.lfq)
# 2010 proteins; 31.25%

# It seems as though a high proportion of proteins are unique or rare. 
# 
# Plot the intersects between all 30 samples (5 replicates; 6 anther lengths)
# using UpSetR.
# Each sample will be a set. Store the sample names (names for protein LFQ 
# intensity columns) as a list.
setnames <- colnames(t.boolean.lfq[,1:30])

upset(t.boolean.lfq, # data frame 
      order.by = "freq", # sort by largest to smallest
      nintersects = 25, # only plot the top 25 intersects
      sets = setnames, # sets to look for intersects
      keep.order = TRUE) # Keep the order of sets in setnames

# Rep B6 at 9, 10, & 11mm has a very high number of completely unique proteins.
# Rep B5 at 6mm is missing a number of proteins present in all other samples. 
# Lets start by removing the completely unique proteins. 
non.unique <- subset(t.boolean.lfq,
                     sample.spread > 1)

# How do the intersects look now?
upset(non.unique,
      order.by = "freq",
      nintersects = 25,
      sets = setnames,
      keep.order = TRUE)

# Still some batch specificity in technical reps. Looks like a particular problem in B6.

# Lets look again for batch specificity.
row.names(boolean.lfq)
colnames(boolean.lfq[,1:5])

# # Miriam filters out B6_06 here.
# boolean.lfq.batch <- boolean.lfq.batch[-25,]

boolean.lfq.batch <- aggregate(boolean.lfq[,4:6434],
                               by=list(Biorep=boolean.lfq$batch),
                               FUN=sum)
rownames(boolean.lfq.batch) <- boolean.lfq.batch$Biorep

for (x in 2:ncol(boolean.lfq.batch)) { # for each column in columns 4 to the last column in the table
  boolean.lfq.batch[,x] <- ifelse(boolean.lfq.batch[,x]>0,1,0) # if the value is more than 0 change it to 1;
  # otherwise keep it as 0.
}

rownames(boolean.lfq.batch)
colnames(boolean.lfq.batch[,1:10])

# Use gather to give a count for each protein across all reps for a given anther
# length
long.batch <- gather(boolean.lfq.batch, # input data frame to reshape
                    key="BaMP", # name for the new key column (name of gathered columns)
                    value = "count", # name for the new value column (protein LFQ intensity)
                    2:6432) # which columns to gather (all protein LFQ intensity columns)
long.batch[1:10,]
colnames(boolean.lfq.batch[,6420:6432])

##############################################################################

# Transpose this presence/abscence table and generate a count for each protein
# with the number of times a protein is measured in the entire dataset.
t.boolean.batch <- t(boolean.lfq.batch[,2:6432])
t.boolean.batch <- as.data.frame(t.boolean.batch)


# Plot the intersects between all 30 samples (5 replicates; 6 anther lengths)
# using UpSetR.
# Each sample will be a set. Store the sample names (names for protein LFQ 
# intensity columns) as a list.
setnames <- colnames(t.boolean.batch[,1:5])

repupset <- upset(t.boolean.batch, # data frame 
      order.by = "freq", # sort by largest to smallest
      nintersects = 25, # only plot the top 25 intersects
      sets = setnames, # sets to look for intersects
      keep.order = TRUE,
      text.scale = 1.2) # Keep the order of sets in setnames

#############################################################################

# Get a list of protein IDs to keep--those present in at least three samples 
# in any biorep.
keep <- subset(long.batch,
               count>=3)
# remove duplicate protein IDs
keep <- keep %>% distinct(BaMP)
keep <- keep$BaMP
# subset only the proteins in our keep list.
t.boolean.lfq.filtered <- subset(t.boolean.lfq,
                                 rownames(t.boolean.lfq) %in% keep)

# Theoretically we've removed batch effects, at least in terms of presence abscence.
# Plot sample intersects of filtered data. 
setnames <- colnames(t.boolean.lfq.filtered[1:30])
upset(t.boolean.lfq.filtered,
      order.by = "freq",
      nintersects = 25,
      sets = setnames,
      keep.order = TRUE)

# Get a list of column numbers which we want to keep
col.num <- which(colnames(prot.vs.lfq) %in% keep)
# Subset from orginal dataset the proteins we want to keep.
filtered.lfq <- prot.vs.lfq[,c(1:3, # keep the sample, batch, and anthre length columns 
                               col.num)]

# Write out the filtered protein data.
write.csv(filtered.lfq,
          "output_tables/prot.vs.lfq.filtered.csv",
          row.names = FALSE)


