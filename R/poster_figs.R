
library(SingleCellExperiment)
library(Features)
library(magrittr)
library(ggplot2)
library(export)
source("R/utils.R")

####---- Load the data ----####

## The data was downloaded from https://scope2.slavovlab.net/docs/data

## Load metadata
samps <- t(read.csv("data/Cells.csv", row.names = 1))

## Load data
dat <- read.csv(file = "data/Peptides-raw.csv")
rownames(dat) <- dat[, 2]
colnames(dat)[1:2] <- c("ProteinAccession", "PeptideSequence")

## Create the SingleCellExperiment object
sce_pep <- SingleCellExperiment(assays = list(as.matrix(dat[, -(1:2)])), 
                                rowData = dat[, 1:2],
                                colData = samps[colnames(dat[, -(1:2)]), ])


####---- Data manipulation ----####

Features(experiments = list(peptide = sce_pep), 
         colData = colData(sce_pep)) %>%
  scp_normalizeFeatures("peptide", fcol = "PeptideSequence",
                        FUN = "mean", na.rm = TRUE, 
                        name = "peptide_norm")  %>%
  scp_aggregateFeatures("peptide_norm", fcol = "ProteinAccession", 
                        name = "protein", fun = colMeans, na.rm = TRUE)  %>%
  scp_normalizeFeatures("protein", fcol = "ProteinAccession",
                        FUN = "mean", na.rm = TRUE, 
                        name = "protein_fnorm") %>%
  scp_normalizeSamples("protein_fnorm", fcol = "ProteinAccession",
                       FUN = "median", na.rm = TRUE, 
                       name = "protein_snorm") %>%
  scp_impute("protein_fnorm", fcol = "ProteinAccession",
             name = "protein_imputed", method = KNNimputation, k = 3) %>%
  scp_batchCorrect("protein_imputed", fcol = "ProteinAccession",
                   name = "protein_batchCorrected",
                   batch = "raw.file", target = "celltype") -> fts
  
  
####---- Plot UMAP of the processed protein data ----####

set.seed(1234)
library(scater)
exp <- "protein_batchCorrected"
fts[[exp]] <- runUMAP(fts[[exp]], 
                      exprs_values = 1, name = "UMAP",
                      ncomponents = 20, scale = TRUE)
ggplot(data.frame(UMAP = reducedDim(fts[[exp]], "UMAP"),
                  CellType = colData(fts[[exp]])$celltype,
                  Batch = colData(fts[[exp]])$batch_chromatography)) + 
  geom_point(aes(x = UMAP.1, y = UMAP.2, 
                 col = CellType, shape = Batch)) + 
  scale_color_manual(name = "Cell type", 
                     values = c(sc_m0 = "skyblue3",
                                sc_u = "coral"),
                     labels = c(sc_m0 = "macrophages",
                                sc_u = "monocytes")) + 
  ggtitle("UMAP on batch corrected protein expression data")



plotUMAP(fts[[exp]], colour_by = "celltype", shape_by = "batch_chromatography")


####---- Plot PCA of the processed protein data ----####

exp <- "protein_batchCorrected"
fts[[exp]] <- runPCA(fts[[exp]], 
                     exprs_values = 1, name = "PCA",
                     ncomponents = 2, scale = TRUE)
plotPCA(fts[[exp]], 
        colour_by = "celltype", shape_by = "batch_chromatography")


####---- Plot expression values for a single protein ----####

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
graph2png(file = "figs/annexinA2.png", height = 4, 
          width = 7)

####---- Plot batch effect ----####





