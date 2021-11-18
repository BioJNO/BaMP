# Join GO and LFQ data tables -------------------------------------------------

# Load in the data.
GO <- read.csv("input_tables/GO_annotations.csv")
LFQ <- read.csv("input_tables/Protein_LFQs.csv")
colnames(GO)
colnames(LFQ)

# GO annotation data has both Mikado and BaMP protein names. LFQ data has only
# Mikado names. We want BaMP IDs for the 6431 LFQ data proteins instead of 
# Mikado. Both data frames share the Mikado name column so we can merge them
# using the values in this column. 

# Check that there are no LFQ data Mikado IDs which are not present in the GO
# annotation data. 
setdiff(LFQ$Mikado, GO$Mikado)

# There are no Mikados in LFQ that aren't in the GO data. Merge the data frames. 
prot.data <- merge(GO, # data frame 1 
                   LFQ, # data frame 2 
                   by="Mikado") # shared column with shared values

# NB: we did not want the GO annotation data which was not present in the LFQ
# data frame. If we did we would have added "all.x=TRUE" to the merge function
# call. Similarly, if the GO annotation file had been missing some ID's we could
# have kept the LFQ data for these proteins by adding "all.y=TRUE" or all data
# from both data frames by adding "all=TRUE". 

# Write out the merged data frame. 
write.csv(prot.data, # data to write
          "output_tables/Proteomic_data.csv", # where to write/what to call it 
          row.names = FALSE) # don't write out the row names.
