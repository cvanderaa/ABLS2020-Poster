
library(SingleCellExperiment)
library(Features)
library(magrittr)
library(ggplot2)
library(export)
library(cowplot)
library(grid)
source("R/utils.R")

####---- Load the data ----####

## The data was downloaded from https://scope2.slavovlab.net/docs/data

## Load metadata
samps <- t(read.csv("data/Cells.csv", row.names = 1))

## Load data
dat <- read.csv(file = "data/Peptides-raw.csv")
rownames(dat) <- dat[, 2]
colnames(dat)[1:2] <- c("ProteinAccession", "PeptideSequence")

assayDat <- as.matrix(dat[, -(1:2)])
rowDat <- dat[, 1:2]
colDat <- samps[colnames(dat[, -(1:2)]), ]

## Create the SingleCellExperiment object
sce_pep <- SingleCellExperiment(assays = list(assayDat), 
                                rowData = rowDat,
                                colData = colDat)


####---- Data manipulation ----####

Features(experiments = list(peptide = sce_pep), 
         colData = colData(sce_pep)) %>%
  scp_normalizeFeatures("peptide", fcol = "PeptideSequence",
                        FUN = "mean", na.rm = TRUE, 
                        name = "peptide_norm")  %>%
  scp_aggregateFeatures("peptide_norm", fcol = "ProteinAccession", 
                        fun = colMeans, na.rm = TRUE,
                        name = "protein")  %>%
  scp_normalizeFeatures("protein", fcol = "ProteinAccession",
                        FUN = "mean", na.rm = TRUE, 
                        name = "protein_fnorm") %>%
  scp_normalizeSamples("protein_fnorm", fcol = "ProteinAccession",
                       FUN = "median", na.rm = TRUE, 
                       name = "protein_snorm") %>%
  scp_impute("protein_fnorm", fcol = "ProteinAccession",
             method = KNNimputation, k = 3,
             name = "protein_imputed") %>%
  scp_batchCorrect("protein_imputed", fcol = "ProteinAccession",
                   batch = "raw.file", target = "celltype",
                   name = "protein_batchCorrected") -> 
  fts
  
  
####---- Plot UMAP of the processed protein data ----####

library(scater)
exp <- "protein_batchCorrected"
fts[[exp]] <- runPCA(fts[[exp]], 
                     exprs_values = 1, name = "PCA",
                     ncomponents = 2, scale = TRUE)
set.seed(1234)
fts[[exp]] <- runUMAP(exprs_values = 1, name = "UMAP",
                      fts[[exp]], dimred = "PCA",
                      ncomponents = 2, scale = TRUE)
pDimred <- ggplot(data.frame(dimred = rbind(reducedDim(fts[[exp]], "PCA"),
                                 reducedDim(fts[[exp]], "UMAP")),
                  dimredType = rep(c("PCA", "UMAP"), each = ncol(fts[[exp]])),
                  CellType = colData(fts[[exp]])$celltype,
                  Batch = colData(fts[[exp]])$batch_chromatography)) + 
  geom_point(aes(x = dimred.PC1, y = dimred.PC2, 
                 col = CellType, shape = Batch)) + 
  scale_color_manual(name = "Cell type", 
                     values = c(sc_m0 = "skyblue3",
                                sc_u = "coral"),
                     labels = c(sc_m0 = "macrophages",
                                sc_u = "monocytes")) + 
  facet_wrap(~ dimredType, scales = "free") + 
  xlab("Dimension 1") + ylab("Dimension 2") + 
  ggtitle("Batch corrected protein expression data") + 
  theme(plot.title = element_text(face = "bold"))

graph2png(pDimred, file = "figs/dimred.png", height = 3.5, 
          width = 7)

####---- Plot leveled expression data ----####

test <- fts[, , names(fts)[1:3]]
sub <- test["P07355", sampleMap(test)$colname[1:200], ]
sub_df <- data.frame(longFormat(sub, colDataCols = c("celltype")))
sub_df$type <- ifelse(sub_df$celltype == "sc_u", 1, 2)
sub_df$UMAP <- reducedDim(fts[[exp]], "UMAP")[sub_df$colname, 2]
sub_df <- sub_df[sub_df$assay %in% c("peptide", "protein"), ]
ggplot(data = sub_df,
       aes(x = reorder(colname, type),
           y = value,
           group = rowname)) +
  xlab("Cells") + ylab("Relative expression value") +
  geom_line(size = 0.5) + 
  geom_point(aes(col = celltype), size = 0.5) + 
  facet_grid(~ assay) + 
  ggtitle("Annexin A2 peptide and protein expression") + 
  scale_color_manual(name = "Cell type", 
                     values = c(sc_m0 = "skyblue3",
                                sc_u = "coral"),
                     labels = c(sc_m0 = "macrophages",
                                sc_u = "monocytes")) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        legend.position = "bottom")
graph2png(file = "figs/annexinA2.png", height = 4, width = 7)
## After discussion with Laurent, this figure is too messy. I will rather use
## the figure from the processing vignette
data(hlpsms)
hl <- readFeatures(hlpsms, ecol = 1:10, name = "psms")
hl <- aggregateFeatures(hl, "psms", "Sequence", name = "peptides", fun = colMeans)
hl <- aggregateFeatures(hl, "peptides", "ProteinGroupAccessions", name = "proteins", fun = colMeans)
sub <- hl["P42227-2", , ]
sub_df <- data.frame(longFormat(sub))
sub_df$assay <- factor(sub_df$assay,
                       labels = c("PSM", "Peptide", "Protein"),
                       levels = c("psms", "peptides", "proteins"))
sub_df$colname <- gsub("X", "", sub_df$colname)
pStat <- ggplot(data = sub_df,
                aes(x = colname,
                    y = value,
                    group = rowname)) +
  geom_line(col = "grey20") + 
  geom_point(col = "#b3ba82") +
  ggtitle("Three levels of Stat3 expression") + 
  facet_grid(~ assay) + 
  xlab("TMT channel") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(face = "bold"))
graph2png(pStat, file = "figs/Stat.png", height = 3, width = 7)


####---- Plot batch effect ----####

test <- fts[["peptide"]]
nbatch <- 3
ncell <- 4
totcell <- nbatch * ncell
test <- test[, 1:(totcell)]
test <- test[do.call(order, 
                     lapply(seq(1, totcell, by = ncell), 
                            function(i) rowSums(is.na(assays(test)[[1]][, i:(i+ncell-1)])))), ]

image(t(assays(test)[[1]][seq(3200, 1, -30), ]), axes = FALSE, 
      xlab = "", ylab = "Features (peptides)", 
      mgp = c(1,0,0), 
      col = colorRampPalette(c(rgb(221/256, 230/256, 163/256), 
                               rgb(122/256, 128/256, 89/256)))(1000),
      main = "Peptide expression data",
      useRaster = TRUE)
marg <- 1 / nbatch / 2
mtext(text = paste0("Batch ", 1:3), side = 1, line = 0, 
      at = marg + 2 * marg * ((1:nbatch) - 1))
graph2png(file = "figs/batchEffect.png", height = 5, 
          width = 5)


####---- Missingness ----####

sc <- sce_pep
# Format the data
df <- do.call(cbind, lapply(unique(colData(sc)$celltype), function(x){
  .sub <- assays(sc)[[1]][, colData(sc)$celltype == x]
  mis <- rowSums(is.na(.sub))/ncol(.sub)*100
  logFC <- apply(.sub, 1, median, na.rm = TRUE)
  out <- data.frame(mis, logFC)
  colnames(out) <- paste0(c("mis", "logFc"), "_", x)
  return(out)
}))
df$relFC <- df$logFc_sc_m0 - df$logFc_sc_u
df$relFC[df$relFC > 2] <- 2
df$relFC[df$relFC < -2] <- -2
# Scatter plot
sp <- ggplot(data = df, aes(x = mis_sc_m0, y = mis_sc_u, col = relFC)) +
  geom_point(size = 0.5) + 
  scale_color_gradient2(name = "log2(rFC)", low = "darkgreen", breaks = c(-2,2), 
                        labels = c("monocyte", "macrophage"),
                        high = "red3", midpoint = 0, mid = "wheat") +
  theme(legend.direction = "horizontal",plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.background = element_blank()) +
  ylab("Missingness (%) in monocytes") + xlab("Missingness (%) in macrophages")
# Macrophage density plot
dp1 <- ggplot(data = df, aes(mis_sc_m0)) + 
  geom_density(fill = "grey") +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 1), "cm"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        panel.background = element_blank()) + 
  ggtitle("Peptide expression data")
# monocyte density plot
dp2 <- ggplot(data = df, aes(mis_sc_u)) + 
  geom_density(fill = "grey") +
  coord_flip() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 1, 0), "cm"),
        panel.background = element_blank())
# Get legend
leg <- get_legend(sp)
sp <- sp + theme(legend.position = "none")
# Empty plot
blank <- ggplot() + geom_blank() + theme(panel.background = element_blank())
# Combine and save plots
pl.mis <- plot_grid(dp1, blank, sp, dp2, align = "none",
                    ncol = 2, nrow = 2,
                    rel_widths=c(10, 1), rel_heights=c(1, 6))
graph2png(x = pl.mis, file = "./figs/missing.png", aspectr = 1,
          width = 4, height = 4)
dev.off()
graph2png(x = grid.draw(leg), file = "./figs/missing-leg.png", 
          width = 4, height = 1)


