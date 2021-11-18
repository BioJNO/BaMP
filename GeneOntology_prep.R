GeneOntology <- read.csv("GeneOntology.csv")
colnames(GeneOntology)

for GO in GeneOntology$id {
  children = separate_rows(merge_GOs,7,sep = ", ")
  colnames(children)
  children <- children$children
}