
library(plyr)

prot.scrs <- read.csv("log2.prot.scrs.csv")

adjsig <- subset(prot.scrs, p.adj<0.01)

GeneOntology <- read.csv("GeneOntology.csv")
colnames(GeneOntology)
colnames(adjsig)
GeneOntology <- rename(GeneOntology,
                       c("id" = "GO"))

adjsig$GO <- gsub('\\.', ':', adjsig$GO)

merge_GOs <- merge(adjsig, GeneOntology, by="GO")

#write.csv(merge_GOs, "significant_biological_processes.csv")

colnames(adjsig)

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
                label = GO,
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

#ggsave("log2.nmds.GO.fit.all.svg", width=183, units = "mm")

#write.csv(merge_GOs, "merge_gos.csv")

colnames(merge_GOs)

merge_GOs[1:10,]$GO

GO_to_path <- read.csv("GO_to_upstream.csv")

colnames(GO_to_path)

merge_GOs_full_path <- merge(merge_GOs, GO_to_path, by="GO")

colnames(merge_GOs_full_path)

# children = separate_rows(merge_GOs,12,sep = ", ")
# colnames(children)
# children <- children$children

# parents = separate_rows(merge_GOs,19,sep = ", ")
# parents <- parents$parent_ids

library(tidyr)
path = separate_rows(merge_GOs_full_path,25,sep = ", ")
above <- path$Upstream

none_above <- merge_GOs_full_path[!(merge_GOs$GO %in% above), ]

# no_parent <- merge_GOs[!(merge_GOs$GO %in% parents), ]
# no_children <- merge_GOs[!(merge_GOs$GO %in% children), ]
# 
# ind <- no_children[!(no_children$GO %in% parents), ]

levels(no_children$category_name)

biol_bottom <- subset(none_above, category_name=="biological_process")
cell_bottom <- subset(none_above, category_name=="cellular_component")
molecular_bottom <- subset(none_above, category_name=="molecular_function")

p <- ggplot(scrs) +
  expand_limits(x = 0.015) +
  geom_segment(data = biol_bottom,
               aes(x = 0, xend=(NMDS1*0.01),
                   y = 0, yend=(NMDS2*0.01)),
               arrow = arrow(length = unit(0.25, "cm")),
               colour = "grey") +
  geom_text(data = biol_bottom,
            aes(x=(NMDS1*0.01),
                y = (NMDS2*0.01),
                label = GO,
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


#ggsave("log2.nmds.GO.fit.bottom.level.biological.svg")

