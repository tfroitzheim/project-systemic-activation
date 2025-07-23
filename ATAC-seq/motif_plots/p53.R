library(dplyr)
library(pheatmap)
library(viridis)
load("./Jun102024_6motifs.RData")

P53 <- P53[!duplicated(P53$transcriptId), ]
P53 <- P53[!grepl("^AMEX60", P53[, "transcriptId"]), ]
row.names(P53) <- NULL


# save.image("./Jun102024_6motifs.RData")

P53_filtered <- P53[rowSums(P53[3:18])>100,]

intact_cols <- select(P53_filtered, ends_with('intact'))
contra_cols <- select(P53_filtered, ends_with('contra'))

P53_lfc <- log2((contra_cols+1)/(intact_cols+1))
P53_lfc$transcriptId <- P53_filtered$transcriptId

P53_fibro <- P53_lfc[1:8]
P53_fibro <- abs(P53_fibro)
P53_fibro$max <- apply(P53_fibro, 1, max)
P53_fibro$rank <- rank(P53_fibro$max)
P53_fibro$transcriptId <- P53_filtered$transcriptId




P53_top <- P53_fibro %>%
  arrange(desc(rank)) %>%
  head(100)

P53_counts <- as.data.frame(P53_top$transcriptId)
names(P53_counts) <- "transcriptId"
P53_counts <- left_join(P53_counts, P53[2:18])


pheatmap(SMAD3_fibro[2:17],
         scale="row",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "Genes with P53 Binding Site from 4C single cell data",
         color = plasma(26),
         #color =  colorRampPalette(c("black", "purple", "yellow"))(10000),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", P53_counts$transcriptId)),
         treeheight_row = 100,
         #breaks = seq(0, 6, by = 1),
         #breaks = NA,
         treeheight_col =  5,
)
