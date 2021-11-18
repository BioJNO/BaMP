library(tidyverse)

filenames <- list.files("ProteusComparisonTables",
                        pattern = "*.csv",
                        full.names = TRUE)

ldf <- lapply(filenames,
              read.csv)

ldf <- lapply(ldf,
              function (x) { x[,-c(1, 9:13)]})


filenames <- gsub("ProteusComparisonTables/2020-05-06_DifferentialAbundanceOfProteins",
                  "",
                  filenames)

filenames <- gsub(".csv",
                  "",
                  filenames)

ldf <- ldf[-c(3,5,15)]
filenames <- filenames[-c(3,5,15)]

result <- mapply(cbind,
              ldf,
              "comparison" = filenames,
              SIMPLIFY = FALSE)

merged.data.frame = Reduce(function(...) merge(..., all=TRUE), result)

sup3 <- read.csv("SupplementaryTable3.csv")
prot.data <- read.csv("output_tables/Proteomic_data.csv")

stats <- rename(merged.data.frame, 
                "Mikado" = "protein")

colnames(prot.data)

prot.data <- rename(prot.data,
                    "DAP" = "BaMP")

prot.data$DAP <- gsub("BaMP",
                       "BaAP",
                       prot.data$DAP)

prot.data <- prot.data[,1:3]

merged <- merge(stats,
                prot.data,
                by = "Mikado")

setdiff(stats$Mikado, prot.data$Mikado)
missing <- setdiff(stats$Mikado, prot.data$Mikado)

missing <- stats[stats$Mikado %in% missing, ]

missing <- separate_rows(missing,
                         1,
                         sep = ";",
                         convert = FALSE)

stats <- rbind(stats, missing)

merged <- merge(stats,
                prot.data,
                by = "Mikado")

setdiff(merged.data.frame$protein, merged$Mikado)

merged2 <- merge(sup3,
                 merged,
                 by = "DAP")

missing2 <- setdiff(merged$DAP, sup3$DAP)
missing2 <- merged[merged$DAP %in% missing2,]

setdiff(sup3$DAP, merged$DAP)

colnames(merged2)

merged2 <- merged2[,c(1,10,2,4:9)]

colnames(merged2)

write.csv(merged2,
          "Sup3.csv",
          row.names = FALSE)

