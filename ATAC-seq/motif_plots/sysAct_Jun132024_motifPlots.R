seu <- readRDS("~/Projects/sysAct_allFiles/sysAct4C_feb29_2024.rds")
load("~/Projects/sysAct_allFiles/newATAC_PCA_SIDU/Jun102024_6motifs.RData")

SMAD3_filtered <- SMAD3[rowSums(SMAD3[3:18])>100,]

intact_cols <- select(SMAD3_filtered, ends_with('intact'))
contra_cols <- select(SMAD3_filtered, ends_with('contra'))

SMAD3_lfc <- log2((contra_cols+1)/(intact_cols+1))
SMAD3_lfc$transcriptId <- SMAD3_filtered$transcriptId

SMAD3_fibro <- SMAD3_lfc[1:8]
SMAD3_fibro <- abs(SMAD3_fibro)
SMAD3_fibro$max <- apply(SMAD3_fibro, 1, max)
SMAD3_fibro$rank <- rank(SMAD3_fibro$max)
SMAD3_fibro$transcriptId <- SMAD3_filtered$transcriptId

SMAD3_top <- SMAD3_fibro %>%
  arrange(desc(rank)) %>%
  head(100)

SMAD3_counts <- as.data.frame(SMAD3_top$transcriptId)
names(SMAD3_counts) <- "transcriptId"
SMAD3_counts <- left_join(SMAD3_counts, SMAD3[2:18])

pheatmap(SMAD3_counts[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "100 top genes with SMAD3 Binding Site from 4C single cell data",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(400),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", SMAD3_counts$transcriptId)),
         treeheight_row = 100,
         breaks = seq(1, 400, by = 1),
         treeheight_col =  5,
)


SMAD3_top50 <- SMAD3_fibro %>%
  arrange(desc(rank)) %>%
  head(50)

SMAD3_counts50 <- as.data.frame(SMAD3_counts$transcriptId)
names(SMAD3_counts50) <- "transcriptId"
SMAD3_counts50 <- left_join(SMAD3_counts50, SMAD3[2:18])
scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
scl <- scale_values(SMAD3_counts50[2:17])
pheatmap(scale_values(SMAD3_counts50[2:17]),
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "50 top genes with SMAD3 Binding Site from 4C single cell data",
         color = viridis(2000),
         #color =  colorRampPalette(c("purple", "yellow"))(20),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", SMAD3_counts50$transcriptId)),
         treeheight_row = 100,
         #breaks = seq(0, 1, by = .004),
         treeheight_col =  5,
)

SMAD3_top5000 <- SMAD3_fibro %>%
  arrange(desc(rank)) %>%
  head(5000)

SMAD3_counts5000 <- as.data.frame(SMAD3_top5000$transcriptId)
names(SMAD3_counts5000) <- "transcriptId"
SMAD3_counts5000 <- left_join(SMAD3_counts5000, SMAD3[2:18])
# set the top smad3 and p53 top 50 and top 100.  No scale, look for adr

SMAD3_counts5000_pt1 <- head(SMAD3_counts5000,1125)
SMAD3_counts5000_pt2 <- head(SMAD3_counts5000,2250) %>% tail(1125)
SMAD3_counts5000_pt3 <- head(SMAD3_counts5000,3375) %>% tail(1125)
SMAD3_counts5000_pt4 <- tail(SMAD3_counts5000,1131)

pheatmap(SMAD3_counts5000_pt1[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "All genes with SMAD3 Binding Site from 4C single cell data pt1",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(2000),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", SMAD3_counts5000_pt1$transcriptId)),
         treeheight_row = 100,
         breaks = seq(0, 2000, by = 1),
         treeheight_col =  5,
         fontsize_row = 7)

pheatmap(SMAD3_counts5000_pt2[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "All genes with SMAD3 Binding Site from 4C single cell data pt2",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(2000),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", SMAD3_counts5000_pt2$transcriptId)),
         treeheight_row = 100,
         breaks = seq(0, 2000, by = 1),
         treeheight_col =  5,
         fontsize_row = 7)


pheatmap(SMAD3_counts5000_pt3[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "All genes with SMAD3 Binding Site from 4C single cell data pt3",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(2000),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", SMAD3_counts5000_pt3$transcriptId)),
         treeheight_row = 100,
         breaks = seq(0, 2000, by = 1),
         treeheight_col =  5,
         fontsize_row = 7)

pheatmap(SMAD3_counts5000_pt4[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "All genes with SMAD3 Binding Site from 4C single cell data pt4",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(2000),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", SMAD3_counts5000_pt4$transcriptId)),
         treeheight_row = 100,
         breaks = seq(0, 2000, by = 1),
         treeheight_col =  5,
         fontsize_row = 7)





# Start for TWIST2
TWIST2_filtered <- TWIST2[rowSums(TWIST2[3:18])>100,]

intact_cols <- select(TWIST2_filtered, ends_with('intact'))
contra_cols <- select(TWIST2_filtered, ends_with('contra'))

TWIST2_lfc <- log2((contra_cols+1)/(intact_cols+1))
TWIST2_lfc$transcriptId <- TWIST2_filtered$transcriptId

TWIST2_fibro <- TWIST2_lfc[1:8]
TWIST2_fibro <- abs(TWIST2_fibro)
TWIST2_fibro$max <- apply(TWIST2_fibro, 1, max)
TWIST2_fibro$rank <- rank(TWIST2_fibro$max)
TWIST2_fibro$transcriptId <- TWIST2_filtered$transcriptId

TWIST2_top100 <- TWIST2_fibro %>%
  arrange(desc(rank)) %>%
  head(100)

TWIST2_counts <- as.data.frame(TWIST2_top100$transcriptId)
names(TWIST2_counts) <- "transcriptId"
TWIST2_counts <- left_join(TWIST2_counts, TWIST2[2:18])

pheatmap(TWIST2_counts[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "100 top genes with TWIST2 Binding Site from 4C single cell data",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(400),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", TWIST2_counts$transcriptId)),
         treeheight_row = 100,
         breaks = seq(1, 400, by = 1),
         treeheight_col =  5,
)


TWIST2_top50 <- TWIST2_fibro %>%
  arrange(desc(rank)) %>%
  head(50)

TWIST2_counts <- as.data.frame(TWIST2_top50$transcriptId)
names(TWIST2_counts) <- "transcriptId"
TWIST2_counts <- left_join(TWIST2_counts, TWIST2[2:18])
pheatmap(TWIST2_counts[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "50 top genes with TWIST2 Binding Site from 4C single cell data",
         #color = viridis(2000),
         color =  colorRampPalette(c("black", "purple", "yellow"))(400),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", TWIST2_counts$transcriptId)),
         treeheight_row = 100,
         breaks = seq(1, 400, by = 1),
         treeheight_col =  5
)

TWIST2_top5000 <- TWIST2_fibro %>%
  arrange(desc(rank)) %>%
  head(5000)

TWIST2_counts5000 <- as.data.frame(TWIST2_top5000$transcriptId)
names(TWIST2_counts5000) <- "transcriptId"
TWIST2_counts5000 <- left_join(TWIST2_counts5000, TWIST2[2:18])
# set the top TWIST2 and p53 top 50 and top 100.  No scale, look for adr

TWIST2_counts5000_pt1 <- head(TWIST2_counts5000,1000)
TWIST2_counts5000_pt2 <- head(TWIST2_counts5000,2000) %>% tail(1000)
TWIST2_counts5000_pt3 <- head(TWIST2_counts5000,3000) %>% tail(1000)
TWIST2_counts5000_pt4 <- tail(TWIST2_counts5000,973)

pheatmap(TWIST2_counts5000_pt1[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "All genes with TWIST2 Binding Site from 4C single cell data pt1",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(2000),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", TWIST2_counts5000_pt1$transcriptId)),
         treeheight_row = 100,
         breaks = seq(0, 2000, by = 1),
         treeheight_col =  5,
         fontsize_row = 7)

pheatmap(TWIST2_counts5000_pt2[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "All genes with TWIST2 Binding Site from 4C single cell data pt2",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(2000),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", TWIST2_counts5000_pt2$transcriptId)),
         treeheight_row = 100,
         breaks = seq(0, 2000, by = 1),
         treeheight_col =  5,
         fontsize_row = 7)

pheatmap(TWIST2_counts5000_pt3[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "All genes with TWIST2 Binding Site from 4C single cell data pt3",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(2000),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", TWIST2_counts5000_pt3$transcriptId)),
         treeheight_row = 100,
         breaks = seq(0, 2000, by = 1),
         treeheight_col =  5,
         fontsize_row = 7)

pheatmap(TWIST2_counts5000_pt4[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "All genes with TWIST2 Binding Site from 4C single cell data pt4",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(2000),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", TWIST2_counts5000_pt4$transcriptId)),
         treeheight_row = 100,
         breaks = seq(0, 2000, by = 1),
         treeheight_col =  5,
         fontsize_row = 7)







# Start for OCT4
OCT4_filtered <- OCT4[rowSums(OCT4[3:18])>100,]

intact_cols <- select(OCT4_filtered, ends_with('intact'))
contra_cols <- select(OCT4_filtered, ends_with('contra'))

OCT4_lfc <- log2((contra_cols+1)/(intact_cols+1))
OCT4_lfc$transcriptId <- OCT4_filtered$transcriptId

OCT4_fibro <- OCT4_lfc[1:8]
OCT4_fibro <- abs(OCT4_fibro)
OCT4_fibro$max <- apply(OCT4_fibro, 1, max)
OCT4_fibro$rank <- rank(OCT4_fibro$max)
OCT4_fibro$transcriptId <- OCT4_filtered$transcriptId

OCT4_top100 <- OCT4_fibro %>%
  arrange(desc(rank)) %>%
  head(100)

OCT4_counts <- as.data.frame(OCT4_top100$transcriptId)
names(OCT4_counts) <- "transcriptId"
OCT4_counts <- left_join(OCT4_counts, OCT4[2:18])

pheatmap(OCT4_counts[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "100 top genes with OCT4 Binding Site from 4C single cell data",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(1000),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", OCT4_counts$transcriptId)),
         treeheight_row = 100,
         breaks = seq(1, 1000, by = 1),
         treeheight_col =  5,
)


OCT4_top50 <- OCT4_fibro %>%
  arrange(desc(rank)) %>%
  head(50)

OCT4_counts <- as.data.frame(OCT4_top50$transcriptId)
names(OCT4_counts) <- "transcriptId"
OCT4_counts <- left_join(OCT4_counts, OCT4[2:18])
pheatmap(OCT4_counts[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "50 top genes with OCT4 Binding Site from 4C single cell data",
         #color = viridis(2000),
         color =  colorRampPalette(c("black", "purple", "yellow"))(400),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", OCT4_counts$transcriptId)),
         treeheight_row = 100,
         breaks = seq(1, 400, by = 1),
         treeheight_col =  5
)

OCT4_top5000 <- OCT4_fibro %>%
  arrange(desc(rank)) %>%
  head(5000)

OCT4_counts5000 <- as.data.frame(OCT4_top5000$transcriptId)
names(OCT4_counts5000) <- "transcriptId"
OCT4_counts5000 <- left_join(OCT4_counts5000, OCT4[2:18])
# set the top OCT4 and p53 top 50 and top 100.  No scale, look for adr

OCT4_counts5000_pt1 <- head(OCT4_counts5000,1000)

pheatmap(OCT4_counts5000_pt1[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "All genes with OCT4 Binding Site from 4C single cell data pt1",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(2000),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", OCT4_counts5000_pt1$transcriptId)),
         treeheight_row = 100,
         breaks = seq(0, 2000, by = 1),
         treeheight_col =  5,
         fontsize_row = 7)





# Start for OCT6
OCT6_filtered <- OCT6[rowSums(OCT6[3:18])>100,]

intact_cols <- select(OCT6_filtered, ends_with('intact'))
contra_cols <- select(OCT6_filtered, ends_with('contra'))

OCT6_lfc <- log2((contra_cols+1)/(intact_cols+1))
OCT6_lfc$transcriptId <- OCT6_filtered$transcriptId

OCT6_fibro <- OCT6_lfc[1:8]
OCT6_fibro <- abs(OCT6_fibro)
OCT6_fibro$max <- apply(OCT6_fibro, 1, max)
OCT6_fibro$rank <- rank(OCT6_fibro$max)
OCT6_fibro$transcriptId <- OCT6_filtered$transcriptId

OCT6_top100 <- OCT6_fibro %>%
  arrange(desc(rank)) %>%
  head(100)

OCT6_counts <- as.data.frame(OCT6_top100$transcriptId)
names(OCT6_counts) <- "transcriptId"
OCT6_counts <- left_join(OCT6_counts, OCT6[2:18])

pheatmap(OCT6_counts[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "100 top genes with OCT6 Binding Site from 4C single cell data",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(400),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", OCT6_counts$transcriptId)),
         treeheight_row = 100,
         breaks = seq(1, 400, by = 1),
         treeheight_col =  5,
)


OCT6_top50 <- OCT6_fibro %>%
  arrange(desc(rank)) %>%
  head(50)

OCT6_counts <- as.data.frame(OCT6_top50$transcriptId)
names(OCT6_counts) <- "transcriptId"
OCT6_counts <- left_join(OCT6_counts, OCT6[2:18])
pheatmap(OCT6_counts[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "50 top genes with OCT6 Binding Site from 4C single cell data",
         #color = viridis(2000),
         color =  colorRampPalette(c("black", "purple", "yellow"))(400),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", OCT6_counts$transcriptId)),
         treeheight_row = 100,
         breaks = seq(1, 400, by = 1),
         treeheight_col =  5
)

OCT6_top5000 <- OCT6_fibro %>%
  arrange(desc(rank)) %>%
  head(5000)

OCT6_counts5000 <- as.data.frame(OCT6_top5000$transcriptId)
names(OCT6_counts5000) <- "transcriptId"
OCT6_counts5000 <- left_join(OCT6_counts5000, OCT6[2:18])
# set the top OCT6 and p53 top 50 and top 100.  No scale, look for adr

OCT6_counts5000_pt1 <- head(OCT6_counts5000,1000)

pheatmap(OCT6_counts5000_pt1[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "All genes with OCT6 Binding Site from 4C single cell data pt1",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(2000),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", OCT6_counts5000_pt1$transcriptId)),
         treeheight_row = 100,
         breaks = seq(0, 2000, by = 1),
         treeheight_col =  5,
         fontsize_row = 7)




# Start for FOS
FOS_filtered <- FOS[rowSums(FOS[3:18])>100,]

intact_cols <- select(FOS_filtered, ends_with('intact'))
contra_cols <- select(FOS_filtered, ends_with('contra'))

FOS_lfc <- log2((contra_cols+1)/(intact_cols+1))
FOS_lfc$transcriptId <- FOS_filtered$transcriptId

FOS_fibro <- FOS_lfc[1:8]
FOS_fibro <- abs(FOS_fibro)
FOS_fibro$max <- apply(FOS_fibro, 1, max)
FOS_fibro$rank <- rank(FOS_fibro$max)
FOS_fibro$transcriptId <- FOS_filtered$transcriptId

FOS_top100 <- FOS_fibro %>%
  arrange(desc(rank)) %>%
  head(100)

FOS_counts <- as.data.frame(FOS_top100$transcriptId)
names(FOS_counts) <- "transcriptId"
FOS_counts <- left_join(FOS_counts, FOS[2:18])

pheatmap(FOS_counts[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "100 top genes with FOS Binding Site from 4C single cell data",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(400),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", FOS_counts$transcriptId)),
         treeheight_row = 100,
         breaks = seq(1, 400, by = 1),
         treeheight_col =  5,
)


FOS_top50 <- FOS_fibro %>%
  arrange(desc(rank)) %>%
  head(50)

FOS_counts <- as.data.frame(FOS_top50$transcriptId)
names(FOS_counts) <- "transcriptId"
FOS_counts <- left_join(FOS_counts, FOS[2:18])
pheatmap(FOS_counts[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "50 top genes with FOS Binding Site from 4C single cell data",
         #color = viridis(2000),
         color =  colorRampPalette(c("black", "purple", "yellow"))(400),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", FOS_counts$transcriptId)),
         treeheight_row = 100,
         breaks = seq(1, 400, by = 1),
         treeheight_col =  5
)

FOS_top5000 <- FOS_fibro %>%
  arrange(desc(rank)) %>%
  head(5000)

FOS_counts5000 <- as.data.frame(FOS_top5000$transcriptId)
names(FOS_counts5000) <- "transcriptId"
FOS_counts5000 <- left_join(FOS_counts5000, FOS[2:18])
# set the top FOS and p53 top 50 and top 100.  No scale, look for adr

FOS_counts5000_pt1 <- head(FOS_counts5000,1000)

pheatmap(FOS_counts5000_pt1[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "All genes with FOS Binding Site from 4C single cell data pt1",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(2000),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", FOS_counts5000_pt1$transcriptId)),
         treeheight_row = 100,
         breaks = seq(0, 2000, by = 1),
         treeheight_col =  5,
         fontsize_row = 7)


















# Start for PITX1
PITX1_filtered <- PITX1[rowSums(PITX1[3:18])>100,]

intact_cols <- select(PITX1_filtered, ends_with('intact'))
contra_cols <- select(PITX1_filtered, ends_with('contra'))

PITX1_lfc <- log2((contra_cols+1)/(intact_cols+1))
PITX1_lfc$transcriptId <- PITX1_filtered$transcriptId

PITX1_fibro <- PITX1_lfc[1:8]
PITX1_fibro <- abs(PITX1_fibro)
PITX1_fibro$max <- apply(PITX1_fibro, 1, max)
PITX1_fibro$rank <- rank(PITX1_fibro$max)
PITX1_fibro$transcriptId <- PITX1_filtered$transcriptId

PITX1_top100 <- PITX1_fibro %>%
  arrange(desc(rank)) %>%
  head(100)

PITX1_counts <- as.data.frame(PITX1_top100$transcriptId)
names(PITX1_counts) <- "transcriptId"
PITX1_counts <- left_join(PITX1_counts, PITX1[2:18])

pheatmap(PITX1_counts[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "100 top genes with PITX1 Binding Site from 4C single cell data",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(1000),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", PITX1_counts$transcriptId)),
         treeheight_row = 100,
         breaks = seq(1, 1000, by = 1),
         treeheight_col =  5,
)


PITX1_top50 <- PITX1_fibro %>%
  arrange(desc(rank)) %>%
  head(50)

PITX1_counts <- as.data.frame(PITX1_top50$transcriptId)
names(PITX1_counts) <- "transcriptId"
PITX1_counts <- left_join(PITX1_counts, PITX1[2:18])
pheatmap(PITX1_counts[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "50 top genes with PITX1 Binding Site from 4C single cell data",
         #color = viridis(2000),
         color =  colorRampPalette(c("black", "purple", "yellow"))(400),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", PITX1_counts$transcriptId)),
         treeheight_row = 100,
         breaks = seq(1, 400, by = 1),
         treeheight_col =  5
)

PITX1_top5000 <- PITX1_fibro %>%
  arrange(desc(rank)) %>%
  head(5000)

PITX1_counts5000 <- as.data.frame(PITX1_top5000$transcriptId)
names(PITX1_counts5000) <- "transcriptId"
PITX1_counts5000 <- left_join(PITX1_counts5000, PITX1[2:18])
# set the top PITX1 and p53 top 50 and top 100.  No scale, look for adr

PITX1_counts5000_pt1 <- head(PITX1_counts5000,1000)

pheatmap(PITX1_counts5000_pt1[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "All genes with PITX1 Binding Site from 4C single cell data",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(2000),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", PITX1_counts5000_pt1$transcriptId)),
         treeheight_row = 100,
         breaks = seq(0, 2000, by = 1),
         treeheight_col =  5,
         fontsize_row = 7)














# Start for SMAD3
SMAD3_filtered <- SMAD3[rowSums(SMAD3[3:18])>100,]

intact_cols <- select(SMAD3_filtered, ends_with('intact'))
contra_cols <- select(SMAD3_filtered, ends_with('contra'))

SMAD3_lfc <- log2((contra_cols+1)/(intact_cols+1))
SMAD3_lfc$transcriptId <- SMAD3_filtered$transcriptId

SMAD3_fibro <- SMAD3_lfc[1:8]
SMAD3_fibro <- abs(SMAD3_fibro)
SMAD3_fibro$max <- apply(SMAD3_fibro, 1, max)
SMAD3_fibro$rank <- rank(SMAD3_fibro$max)
SMAD3_fibro$transcriptId <- SMAD3_filtered$transcriptId

SMAD3_top100 <- SMAD3_fibro %>%
  arrange(desc(rank)) %>%
  head(100)

SMAD3_counts <- as.data.frame(SMAD3_top100$transcriptId)
names(SMAD3_counts) <- "transcriptId"
SMAD3_counts <- left_join(SMAD3_counts, SMAD3[2:18])

pheatmap(SMAD3_counts[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "100 top genes with SMAD3 Binding Site from 4C single cell data",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(400),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", SMAD3_counts$transcriptId)),
         treeheight_row = 100,
         breaks = seq(1, 400, by = 1),
         treeheight_col =  5,
)


SMAD3_top50 <- SMAD3_fibro %>%
  arrange(desc(rank)) %>%
  head(50)

SMAD3_counts <- as.data.frame(SMAD3_top50$transcriptId)
names(SMAD3_counts) <- "transcriptId"
SMAD3_counts <- left_join(SMAD3_counts, SMAD3[2:18])
pheatmap(SMAD3_counts[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "50 top genes with SMAD3 Binding Site from 4C single cell data",
         #color = viridis(2000),
         color =  colorRampPalette(c("black", "purple", "yellow"))(400),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", SMAD3_counts$transcriptId)),
         treeheight_row = 100,
         breaks = seq(1, 400, by = 1),
         treeheight_col =  5
)

SMAD3_top5000 <- SMAD3_fibro %>%
  arrange(desc(rank)) %>%
  head(5000)

SMAD3_counts5000 <- as.data.frame(SMAD3_top5000$transcriptId)
names(SMAD3_counts5000) <- "transcriptId"
SMAD3_counts5000 <- left_join(SMAD3_counts5000, SMAD3[2:18])
# set the top SMAD3 and p53 top 50 and top 100.  No scale, look for adr

SMAD3_counts5000_pt1 <- head(SMAD3_counts5000,1000)
SMAD3_counts5000_pt2 <- head(SMAD3_counts5000,2000) %>% tail(1000)
SMAD3_counts5000_pt3 <- head(SMAD3_counts5000,3000) %>% tail(1000)
SMAD3_counts5000_pt4 <- head(SMAD3_counts5000,4000) %>% tail(1000)
SMAD3_counts5000_pt5 <- tail(SMAD3_counts5000,840)

pheatmap(SMAD3_counts5000_pt1[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "All genes with SMAD3 Binding Site from 4C single cell data pt1",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(2000),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", SMAD3_counts5000_pt1$transcriptId)),
         treeheight_row = 100,
         breaks = seq(0, 2000, by = 1),
         treeheight_col =  5,
         fontsize_row = 7)

pheatmap(SMAD3_counts5000_pt2[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "All genes with SMAD3 Binding Site from 4C single cell data pt2",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(2000),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", SMAD3_counts5000_pt2$transcriptId)),
         treeheight_row = 100,
         breaks = seq(0, 2000, by = 1),
         treeheight_col =  5,
         fontsize_row = 7)

pheatmap(SMAD3_counts5000_pt3[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "All genes with SMAD3 Binding Site from 4C single cell data pt3",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(2000),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", SMAD3_counts5000_pt3$transcriptId)),
         treeheight_row = 100,
         breaks = seq(0, 2000, by = 1),
         treeheight_col =  5,
         fontsize_row = 7)

pheatmap(SMAD3_counts5000_pt4[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "All genes with SMAD3 Binding Site from 4C single cell data pt4",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(2000),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", SMAD3_counts5000_pt4$transcriptId)),
         treeheight_row = 100,
         breaks = seq(0, 2000, by = 1),
         treeheight_col =  5,
         fontsize_row = 7)

pheatmap(SMAD3_counts5000_pt5[2:17],
         scale="none",
         cluster_rows=T,
         cluster_cols=F,
         show_rownames = T,
         main = "All genes with SMAD3 Binding Site from 4C single cell data pt5",
         #color = plasma(50),
         color =  colorRampPalette(c("black", "purple", "yellow"))(2000),
         border_color = NA,
         labels_row = as.expression(gsub("-(AMEX\\w+)", "", SMAD3_counts5000_pt4$transcriptId)),
         treeheight_row = 100,
         breaks = seq(0, 2000, by = 1),
         treeheight_col =  5,
         fontsize_row = 7)

library(readr)
write_tsv(SMAD3, "/Users/sjblair/Projects/sysAct_allFiles/newATAC_PCA_SIDU/sysAct_bindingMotfs_deliverables/SMAD3_counts.tsv")
write_tsv(FOS, "/Users/sjblair/Projects/sysAct_allFiles/newATAC_PCA_SIDU/sysAct_bindingMotfs_deliverables/FOS_counts.tsv")
write_tsv(OCT4, "/Users/sjblair/Projects/sysAct_allFiles/newATAC_PCA_SIDU/sysAct_bindingMotfs_deliverables/OCT4_counts.tsv")
write_tsv(OCT6, "/Users/sjblair/Projects/sysAct_allFiles/newATAC_PCA_SIDU/sysAct_bindingMotfs_deliverables/OCT6_counts.tsv")
write_tsv(PITX1, "/Users/sjblair/Projects/sysAct_allFiles/newATAC_PCA_SIDU/sysAct_bindingMotfs_deliverables/PITX1_counts.tsv")
write_tsv(TWIST2, "/Users/sjblair/Projects/sysAct_allFiles/newATAC_PCA_SIDU/sysAct_bindingMotfs_deliverables/TWIST2_counts.tsv")
write_tsv(P53, "/Users/sjblair/Projects/sysAct_allFiles/newATAC_PCA_SIDU/sysAct_bindingMotfs_deliverables/P53_counts.tsv")
