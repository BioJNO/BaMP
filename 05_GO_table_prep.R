# Transform data for GO annotations

annotated <- read.csv("Proteomic_data.csv")
colnames(annotated)

# split the Panzzer GO annotation to give a row for each GO.
library(tidyr)
df = separate_rows(annotated,5,sep = ",")
colnames(df)

# Create a table of GO sums per sample
library(reshape2)
df <- df[,c(5,19:48)]
df.melt <- melt(df,
                id="PANNZER_GO",
                variable.name = "sample",
                value.name = "LFQ")

summed <- aggregate(.~PANNZER_GO+sample, df.melt, sum)

wide.go <- pivot_wider(summed, id_cols = sample, names_from = PANNZER_GO, values_from = LFQ)

# remove "LFQ.Intensity." string from sample name column
wide.go$sample <- gsub("LFQ.intensity.", "", wide.go$sample)

# write.csv(wide.go, "wide.nogg.go.csv", row.names = FALSE)

# Create a table of GOs to proteins
df = separate_rows(annotated,5,sep = ",")
colnames(df)

df <- df[,c(3,5)]

library(dplyr)

GO_count <- count(df,PANNZER_GO)
colnames(GO_count)

library(ggplot2)

ggplot(data=GO_count,
       mapping= aes(x=reorder(GO_count$PANNZER_GO, n), n)) +
  geom_bar(stat="identity", aes(fill="#d95f0e")) +
  coord_flip() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank())

GeneOntology <- read.csv("GeneOntology.csv")
colnames(GeneOntology)

library(plyr)
GeneOntology <- rename(GeneOntology,
                       c("id" = "GO"))

GO_count <- rename(GO_count,
                       c("eggNOG_GO" = "GO"))

merge_GOs <- merge(GO_count, GeneOntology, by="GO")
colnames(merge_GOs)

# write.csv(merge_GOs, "GO_count_nogg.csv")

children = separate_rows(merge_GOs,7,sep = ", ")
colnames(children)
children <- children$children

parents = separate_rows(merge_GOs,14,sep = ", ")
parents <- parents$parent_ids

no_children <- merge_GOs[!(merge_GOs$GO %in% children), ]
remove_parents <- merge_GOs[!(merge_GOs$GO %in% parents), ]

# write.csv(no_children, "GO_count_no_children.csv")
# write.csv(remove_parents, "no_parent_ids.csv")

ggplot(data=no_children,
       mapping= aes(x=reorder(PANNZER_GO, n), n)) +
  geom_bar(stat="identity", aes(fill="#d95f0e")) +
  coord_flip() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank())

ggplot(data=remove_parents, mapping= aes(x=reorder(GO, n), n)) +
  geom_bar(stat="identity") + coord_flip()


cyto <- c("GO:0000153", "GO:0000164", "GO:0000177", "GO:0000229", "GO:0000308", "GO:0000407", "GO:0002096", "GO:0002753", "GO:0005739", "GO:0005773", "GO:0005783", "GO:0005793", "GO:0005794", "GO:0005829", "GO:0005850", "GO:0005851", "GO:0005852", "GO:0005853", "GO:0005854", "GO:0005881", "GO:0005883", "GO:0005938", "GO:0005946", "GO:0005948", "GO:0005950", "GO:0005951", "GO:0005960", "GO:0005964", "GO:0005965", "GO:0005968", "GO:0005969", "GO:0006888", "GO:0009317", "GO:0009320", "GO:0009321", "GO:0009328", "GO:0009331", "GO:0009332", "GO:0009333", "GO:0009336", "GO:0009339", "GO:0009340", "GO:0009345", "GO:0009346", "GO:0009348", "GO:0009357", "GO:0009359", "GO:0009382", "GO:0009504", "GO:0009524", "GO:0009525", "GO:0009536", "GO:0009574", "GO:0010169", "GO:0010938", "GO:0016281", "GO:0016528", "GO:0017102", "GO:0017109", "GO:0018444", "GO:0019037", "GO:0019197", "GO:0019812", "GO:0019813", "GO:0020022", "GO:0030094", "GO:0030096", "GO:0030117", "GO:0030700", "GO:0030877", "GO:0030929", "GO:0031201", "GO:0031209", "GO:0031250", "GO:0031251", "GO:0031410", "GO:0031414", "GO:0031680", "GO:0032019", "GO:0032047", "GO:0032068", "GO:0032449", "GO:0032477", "GO:0032527", "GO:0032921", "GO:0033254", "GO:0034081", "GO:0034270", "GO:0034274", "GO:0034430", "GO:0034499", "GO:0034709", "GO:0034715", "GO:0034719", "GO:0034741", "GO:0034743", "GO:0034744", "GO:0034745", "GO:0034746", "GO:0034748", "GO:0034749", "GO:0034750", "GO:0034973", "GO:0034996", "GO:0035550", "GO:0036007", "GO:0036457", "GO:0036464", "GO:0036501", "GO:0042566", "GO:0042579", "GO:0042587", "GO:0042716", "GO:0042735", "GO:0043265", "GO:0043291", "GO:0043292", "GO:0043614", "GO:0044207", "GO:0044222", "GO:0044223", "GO:0044227", "GO:0044312", "GO:0045169", "GO:0045170", "GO:0045239", "GO:0045251", "GO:0045253", "GO:0045254", "GO:0045495", "GO:0048471", "GO:0048492", "GO:0048500", "GO:0055087", "GO:0055107", "GO:0055108", "GO:0060417", "GO:0061908", "GO:0062073", "GO:0070088", "GO:0070214", "GO:0070422", "GO:0070441", "GO:0070695", "GO:0070826", "GO:0070937", "GO:0070992", "GO:0070993", "GO:0071026", "GO:0071075", "GO:0071142", "GO:0071203", "GO:0071212", "GO:0071254", "GO:0071513", "GO:0071521", "GO:0071563", "GO:0072516", "GO:0090120", "GO:0090651", "GO:0090652", "GO:0097433", "GO:0098545", "GO:0098593", "GO:0099568", "GO:0099568", "GO:0106232", "GO:0110158", "GO:0120098", "GO:0120099", "GO:0120217", "GO:1905720", "GO:1990036", "GO:1990316", "GO:1990565", "GO:1990658", "GO:1990728", "GO:1990730", "GO:1990783", "GO:1990917")

cyto_terms <- merge_GOs[(merge_GOs$GO %in% cyto), ]

setdiff(cyto, children)
