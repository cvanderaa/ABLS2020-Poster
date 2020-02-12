
####---- Specht et al. 2019 ---####


# Specht, Harrison, Edward Emmott, Toni Koller, and Nikolai Slavov. 2019. 
# “High-Throughput Single-Cell Proteomics Quantifies the Emergence of Macrophage 
# Heterogeneity.” bioRxiv. https://doi.org/10.1101/665307.

library(Features)
library(SingleCellExperiment)
library(magrittr)
setwd("R/")
source("utils.R")

# The data was downloaded from https://scope2.slavovlab.net/docs/data to 
# scpdata/extdata/specht2019


####---- Experiment metadata ----####

# Create the experiment data for the MSnSet object
expInfo <- list(title = "High-throughput single-cell proteomics quantifies the emergence of macrophage heterogeneity",
                abstract = "The fate and physiology of individual cells are controlled by networks of proteins. Yet, our ability to quantitatively analyze protein networks in single cells has remained limited. To overcome this barrier, we developed SCoPE2. It integrates concepts from Single-Cell ProtEomics by Mass Spectrometry (SCoPE-MS) with automated and miniaturized sample preparation, substantially lowering cost and hands-on time. SCoPE2 uses data-driven analytics to optimize instrument parameters for sampling more ion copies per protein, thus supporting quantification with improved count statistics. These advances enabled us to analyze the emergence of cellular heterogeneity as homogeneous monocytes differentiated into macrophage-like cells in the absence of polarizing cytokines. We used SCoPE2 to quantify over 2,000 proteins in 356 single monocytes and macrophages in about 85 hours of instrument time, and the quantified proteins allowed us to discern single cells by cell type. Furthermore, the data uncovered a continuous gradient of proteome states for the macrophage-like cells, suggesting that macrophage heterogeneity may emerge even in the absence of polarizing cytokines. Our methodology lays the foundation for quantitative analysis of protein networks at single-cell resolution.",
                url = "https://doi.org/10.1101/665307",
                dateStamp = "2019-07-09",
                name = "Harrison Specht, Edward Emmott, David H. Perlman, Antonius Koller, Nikolai Slavov",
                lab = "Slavov Lab",
                instrumentModel = "Q-Exactive Orbitrap",
                instrumentManufacturer = "Thermo Scientific",
                softwareName = "MaxQuant",
                softwareVersion = "1.6.2.3",
                switchingCriteria = "After a precursor scan from 450 to 1600 m/z at 70,000 resolving power, the top 5 most intense precursor ions with charges 2 to 4 and above the AGC min threshold of 20,000 were isolated for MS2 analysis via a 0.7 Th isolation window",
                ionSource = "ESI",
                ionSourceDetails = "Electrospray voltage was set to 2,200V, applied at the end of the analytical column. To reduce atmospheric background ions and enhance peptide signal to noise ratio, an Active Background Ion Reduction Device (ABIRD, by ESI Source Solutons, LLC, Woburn MA, USA) was used at the nanospray interface. The temperature of ion transfer tube was 250 degrees Celsius and the S-lens RF level set to 80.",
                analyser = "orbitrap",
                analyserDetails = "Precursor ions were accumulated for at most 300ms. Then they were fragmented via HCD at a and the fragments analyzed at 70,000 resolving power. Dynamic exclusion was used with a duration of 30 seconds with a mass tolerance of 10ppm.",
                collisionEnergy = "33 eV (normalized to m/z 500, z=1)")


####---- Loading the data ----####

## Load MaxQuant output data 
mqFile <- "../data/ev_updated.txt"
hdr <- colnames(read.table(mqFile, header = TRUE, nrows = 1, sep = "\t")) 
ecol <-  grep("y[.]\\d+$", hdr, value = TRUE)
## Create the SingleCellExperiment object
scp <- readSingleCellExperiment(mqFile, 
                                ecol = ecol, 
                                sep = "\t")
## Clean rowData
fcol <- c(Run = "Raw.file", 
          Sequence = "Sequence", 
          ModifiedSequence =  "Modified.sequence", 
          Charge = "Charge",
          Proteins = "Proteins", 
          Reverse = "Reverse", 
          LeadingRazorProtein = "Leading.razor.protein", 
          PotentialContaminant = "Potential.contaminant", 
          GeneNames = "Gene.names",
          ProteinNames = "Protein.names", 
          Mass = "Mass", 
          CalibratedRT = "Calibrated.retention.time", 
          PIF = "PIF", 
          PEP = "dart_PEP",  
          qVal = "dart_qval", 
          Score = "Score")
rowData(scp) <- rowData(scp)[, fcol]
colnames(rowData(scp)) <- names(fcol)
## Split the SingleCellExperiment object in separate object according to the 
## MS run, that is split by batch
scp <- splitSCE(scp, f = "Run")
## Store the data as a Features object and add the experimental information
scp <- Features(experiments = scp, 
                metadata = list(ExperimentalInfo = expInfo))
## Remove batches that are not single-cell and remove the FP95 and FP96
## experiments
scp <- scp[,, grep(pattern = "col19|col2[0-4]|blank|QC|FP9[56]", 
                   x = names(scp), value = TRUE, invert = TRUE)]

## Add sample information
batch <- read.csv("../data/batch_fp60-97.csv", row.names = 1, 
                  check.names = FALSE)
batchVars <- c(LcBatch = "lcbatch", 
               SortBatch ="sortday", 
               DigestBatch ="digest")
annot <- read.csv("../data/annotation_fp60-97.csv", check.names = FALSE)
annot <- tidyr::pivot_longer(annot, cols = -grep("Set", colnames(annot)), 
                             names_to = "Run", values_to = "type")
annot <- data.frame(type = annot$type, row.names = paste0(annot$Set, annot$Run))
## TODO this is a bad idea as sampleMap is overwritten when adding assays
## Check later if can use colData when MatchedAssayExperiment 
## is relaxed in Features see https://github.com/rformassspectrometry/Features/issues/46
colDat <- batch[as.character(sampleMap(scp)$assay), batchVars]
colDat$SampleType <- annot[paste0(sampleMap(scp)$colname, 
                                  sampleMap(scp)$assay), "type"]

tmp1 <- scp

####---- Clean data ----####

## TODO create issue to make it work on multiple assays at once
for (i in 1:length(scp)) {
  scp <- zeroIsNA(scp, i = i)
  print(i)
}

tmp2 <- scp


####---- Filter PSMs ----####

## Filter based on identification
scp %>%
  filterFeatures(~ Reverse != "+") %>%
  filterFeatures(~ PotentialContaminant != "+") %>%
  ## 
  ## filterFeatures(filter = VariableFilter(LeadingRazorProtein, 
  ##                                        value = "REV",
  ##                                        condition = "startsWith")) %>%
  filterFeatures(~ !PIF %in% NA & PIF > 0.8)
  filterFeatures(~ qVal < 0.01) %>%
  filterFeatures(~ PEP < 0.02) -> scp ## TODO maybe adapt to keep all peptide beloging to proteins having at least 1 significant peptide

## Filter based on quantification 
## TODO think how to improve this
for (i in 1:length(scp)) {
  scp[[i]] <- scp_filterSCR(scp[[i]], samples = 4:11, carrier = 1, thresh = 0.1)
  print(i)
}

tmp3 <- scp


####---- Filter batches ----####

## Keep runs with more than 300 PSMs
scp <- scp[, , nrowAssays(scp) > 300]


####---- Aggregate PSMs to Peptides ----####

for (exp in names(scp)) {
  scp <- aggregateFeaturesSCE(scp, 
                              i = exp, 
                              fcol = "Sequence", 
                              name = paste0(exp, "_peptides"), 
                              fun = colMeans)
}


####---- Combine batches ----####

## Create a new features object containing all batches in a single Features
## object, where columns are single cells and rows are features (PSM)
scp2 <- combineFeatures()



