## Set your working directory
setwd("~/eco_genomics_2026/transcriptomics/mydata")

## Import the libraries that we're likely to need in this session

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)  

####################################################

#### Import our data ####

####################################################

# Import the counts matrix
countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample description table
conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ treatment)
# Filter 
dds <- dds[rowSums(counts(dds) >= 15) >= 28,]

# Subset the DESeqDataSet to the specific level of the "generation" factor
dds_F0 <- subset(dds, select = generation == 'F0')

dds_F4 <- subset(dds, select = generation == 'F4')

dds_all <- subset(dds, select = generation == c('F0', 'F4'))

# Perform DESeq2 analysis on the subset
dds_F0 <- DESeq(dds_F0)

dds_F4 <- DESeq(dds_F4)

dds_all <- DESeq(dds_all)
####################################################

#### Check on the DE results from the DESeq ####

####################################################

## F0 -----------------------------------------------------------------
resultsNames(dds_F0)
# [1] "Intercept"           "treatment_OA_vs_AM"  "treatment_OW_vs_AM"  "treatment_OWA_vs_AM"

res_OWAvsAM_F0 <- results(dds_F0, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM_F0 <- res_OWAvsAM_F0[order(res_OWAvsAM_F0$padj),]
head(res_OWAvsAM_F0)  
summary(res_OWAvsAM_F0)


res_OWvsAM_F0 <- results(dds_F0, name="treatment_OW_vs_AM", alpha=0.05)
res_OWvsAM_F0 <- res_OWvsAM_F0[order(res_OWvsAM_F0$padj),]
head(res_OWvsAM_F0) 
summary(res_OWvsAM_F0)

res_OAvsAM_F0 <- results(dds_F0, name="treatment_OA_vs_AM", alpha=0.05)
res_OAvsAM_F0 <- res_OAvsAM_F0[order(res_OAvsAM_F0$padj),]
head(res_OWAvsAM_F0)
summary(res_OAvsAM_F0)

## F4 ---------------------------------------------------------------
resultsNames(dds_F4)
# [1] "Intercept"           "treatment_OA_vs_AM"  "treatment_OW_vs_AM"  "treatment_OWA_vs_AM"

res_OWAvsAM_F4 <- results(dds_F4, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM_F4 <- res_OWAvsAM_F4[order(res_OWAvsAM_F4$padj),]
head(res_OWAvsAM_F4)  
summary(res_OWAvsAM_F4)


res_OWvsAM_F4 <- results(dds_F4, name="treatment_OW_vs_AM", alpha=0.05)
res_OWvsAM_F4 <- res_OWvsAM_F4[order(res_OWvsAM_F4$padj),]
head(res_OWvsAM_F4) 
summary(res_OWvsAM_F4)

res_OAvsAM_F4 <- results(dds_F4, name="treatment_OA_vs_AM", alpha=0.05)
res_OAvsAM_F4 <- res_OAvsAM_F4[order(res_OAvsAM_F4$padj),]
head(res_OWAvsAM_F4)
summary(res_OAvsAM_F4)

#################################################################

#### PLOT OVERLAPPING DEGS IN VENN EULER DIAGRAM for f0 and f4 ####

#################################################################

# For OW vs AM F0
res_OWvsAM_F0 <- results(dds_F0, name="treatment_OW_vs_AM", alpha=0.05) # pull out the results for the contrast of interest
res_OWvsAM_F0 <- res_OWvsAM_F0[order(res_OWvsAM_F0$padj),] # order them by significance
res_OWvsAM_F0 <- res_OWvsAM_F0[!is.na(res_OWvsAM_F0$padj),] # get rid of any NAs
degs_OWvsAM_F0 <- row.names(res_OWvsAM_F0[res_OWvsAM_F0$padj < 0.05,]) # make a list of significant differentially expressed genes for this contrast

# For OA vs AM F0
res_OAvsAM_F0 <- results(dds_F0, name="treatment_OA_vs_AM", alpha=0.05)
res_OAvsAM_F0 <- res_OAvsAM_F0[order(res_OAvsAM_F0$padj),]
res_OAvsAM_F0 <- res_OAvsAM_F0[!is.na(res_OAvsAM_F0$padj),]
degs_OAvsAM_F0 <- row.names(res_OAvsAM_F0[res_OAvsAM_F0$padj < 0.05,])

# For OWA vs AM F0
res_OWAvsAM_F0 <- results(dds_F0, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM_F0 <- res_OWAvsAM_F0[order(res_OWAvsAM_F0$padj),]
res_OWAvsAM_F0 <- res_OWAvsAM_F0[!is.na(res_OWAvsAM_F0$padj),]
degs_OWAvsAM_F0 <- row.names(res_OWAvsAM_F0[res_OWAvsAM_F0$padj < 0.05,])

# For OW vs AM F4
res_OWvsAM_F4 <- results(dds_F4, name="treatment_OW_vs_AM", alpha=0.05) # pull out the results for the contrast of interest
res_OWvsAM_F4 <- res_OWvsAM_F4[order(res_OWvsAM_F4$padj),] # order them by significance
res_OWvsAM_F4 <- res_OWvsAM_F4[!is.na(res_OWvsAM_F4$padj),] # get rid of any NAs
degs_OWvsAM_F4 <- row.names(res_OWvsAM_F4[res_OWvsAM_F4$padj < 0.05,]) # make a list of significant differentially expressed genes for this contrast

# For OA vs AM F4
res_OAvsAM_F4 <- results(dds_F4, name="treatment_OA_vs_AM", alpha=0.05)
res_OAvsAM_F4 <- res_OAvsAM_F4[order(res_OAvsAM_F4$padj),]
res_OAvsAM_F4 <- res_OAvsAM_F4[!is.na(res_OAvsAM_F4$padj),]
degs_OAvsAM_F4 <- row.names(res_OAvsAM_F4[res_OAvsAM_F4$padj < 0.05,])

# For OWA vs AM F4
res_OWAvsAM_F4 <- results(dds_F4, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM_F4 <- res_OWAvsAM_F4[order(res_OWAvsAM_F4$padj),]
res_OWAvsAM_F4 <- res_OWAvsAM_F4[!is.na(res_OWAvsAM_F4$padj),]
degs_OWAvsAM_F4 <- row.names(res_OWAvsAM_F4[res_OWAvsAM_F4$padj < 0.05,])

# For OWA vs AM all
res_OWAvsAM_all <- results(dds_all, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM_all <- res_OWAvsAM_all[order(res_OWAvsAM_all$padj),]
res_OWAvsAM_all <- res_OWAvsAM_all[!is.na(res_OWAvsAM_all$padj),]
degs_OWAvsAM_all <- row.names(res_OWAvsAM_all[res_OWAvsAM_all$padj < 0.05,])

library(eulerr)
library(gridExtra)
## F0 Euler ---------------------------------------------------------------------------

# Total = how many significant genes are differentially expressed in compared treatments; whole circles of venn
length(degs_OAvsAM_F0)  # 602
length(degs_OWvsAM_F0)  # 5517 
length(degs_OWAvsAM_F0)  # 3918

# Intersections = overlaps of venn diagram
length(intersect(degs_OAvsAM_F0,degs_OWvsAM_F0))  # 444
length(intersect(degs_OAvsAM_F0,degs_OWAvsAM_F0))  # 380
length(intersect(degs_OWAvsAM_F0,degs_OWvsAM_F0))  # 2743

# Shared across all = center of venn
intWA <- intersect(degs_OAvsAM_F0,degs_OWvsAM_F0)
length(intersect(degs_OWAvsAM_F0,intWA)) # 338

# Number unique to each treatment = non-overlapping sections of venn

602-444-380+338 # 116 OA
5517-444-2743+338 # 2668 OW 
3918-380-2743+338 # 1133 OWA

# Number shared in pairs of treatments = adjusted intersects accounting for the middle that is shared among all

444-338 # 106 OA & OW
380-338 # 42 OA & OWA
2743-338 # 2405 OWA & OW

# Now assemble the results
# Note that the names are important and have to be specific to line up the diagram
fit1 <- euler(c("OA" = 116, "OW" = 2668, "OWA" = 1133, "OA&OW" = 106, "OA&OWA" = 42, "OW&OWA" = 2405, "OA&OW&OWA" = 338))

# And make the plot!
p1 <- plot(fit1, main = "F0 DEGS for each Treatment", lty = 1:3, quantities = TRUE,
     fills = c("lightpink", "lightblue2", "thistle"))
# lty changes the lines

#cross check with above lengths of DEGS: the four values, unique, shared with one other, shared with the second other, shared across all treatments, should sum to the length of DEGs for each treatment contrast to AM
2668+2405+338+106 # 5517 total OW
1133+2405+338+42  # 3918 total OWA
116+42+106+338    # 602  total OA
## F4 Euler ---------------------------------------------------------------------------

# Total = how many significant genes are differentially expressed in compared treatments; whole circles of venn
length(degs_OAvsAM_F4)  # 157
length(degs_OWvsAM_F4)  # 153 
length(degs_OWAvsAM_F4)  # 241

# Intersections = overlaps of venn diagram
length(intersect(degs_OAvsAM_F4,degs_OWvsAM_F4))  # 58
length(intersect(degs_OAvsAM_F4,degs_OWAvsAM_F4))  # 41
length(intersect(degs_OWAvsAM_F4,degs_OWvsAM_F4))  # 44

# Shared across all = center of venn
intWA <- intersect(degs_OAvsAM_F4,degs_OWvsAM_F4)
length(intersect(degs_OWAvsAM_F4,intWA)) # 18

# Number unique to each treatment = non-overlapping sections of venn

157-58-41+18 # 76 OA
153-58-44+18 # 69 OW 
241-41-44+18 # 174 OWA

# Number shared in pairs of treatments = adjusted intersects accounting for the middle that is shared among all

58-18 # 40 OA & OW
41-18 # 23 OA & OWA
44-18 # 26 OWA & OW

# Now assemble the results
# Note that the names are important and have to be specific to line up the diagram
fit1 <- euler(c("OA" = 76, "OW" = 69, "OWA" = 174, "OA&OW" = 40, "OA&OWA" = 23, "OW&OWA" = 26, "OA&OW&OWA" = 18))

# And make the plot!
p2 <- plot(fit1, main = "F4 DEGS for each Treatment", lty = 1:3, quantities = TRUE,
     fills = c("lightpink3", "lightblue3", "thistle3"))
# lty changes the lines

#cross check with above lengths of DEGS: the four values, unique, shared with one other, shared with the second other, shared across all treatments, should sum to the length of DEGs for each treatment contrast to AM
2668+2405+338+106 # 5517 total OW
1133+2405+338+42  # 3918 total OWA
116+42+106+338    # 602  total OA

# 2. Arrange into 1 row and 2 columns
grid.arrange(p1, p2, ncol = 2)
#######################################

#### PCA ####

######################################

## Now lets transform our data so we can make a PCA plot
# Check the quality of the data by sample clustering and visualization

# The goal of transformation "is to remove the dependence of the variance on the mean, particularly the high variance of the logarithm of count data when the mean is low."

# this gives log2(n + 1)
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))

# Variance stabilizing transformation
vsd <- vst(dds, blind=FALSE)
meanSdPlot(assay(vsd))
## ABOVE is a standard plot to visualize variance, overall more variation in lower expressed genes

### Look for outliers in sample ###

sampleDists <- dist(t(assay(vsd)))

# cluster samples comparing pairwise,dark is comparing self with self
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$line, vsd$generation, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

data <- plotPCA(vsd, intgroup=c("treatment","generation"), returnData=TRUE)
percentVar <- round(100 * attr(data,"percentVar"))

###########  
# Just looking at F0 generation
dataF0 <- subset(data, generation == 'F0')

F0 <- ggplot(dataF0, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill=treatment, shape=treatment)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylim(-10, 25) + xlim(-40, 10)+ # zoom for F0 with new assembly
  scale_shape_manual(values=c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  scale_fill_manual(values=c("lightblue", "orange2","brown2", "purple"), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank())

# just looking at F0 generation, strong clustering by treatment groups, suggests the treatment groups are most similar to one another
F0
# Above graph portrays developmental plasticity, how do they respond to environment

################################ F4

dataF4 <- subset(data, generation == 'F4')

F4 <- ggplot(dataF4, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill=treatment, shape=treatment)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylim(-40, 25) + xlim(-50, 55)+ # limits with filtered assembly
  scale_shape_manual(values=c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  scale_fill_manual(values=c("lightblue", "orange2","brown2", "purple"), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  guides(shape = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  guides(fill = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank())
# Just looking at F4 generation, largely overlapping, synchronized/homeostasis reached
F4

# put it all together into one plot to look at all 4 plots together
blank <- ggplot() + theme_void()

ggarrange(
  ggarrange(F0, F4, ncol = 2, legend = "right", common.legend = TRUE),
  nrow = 1,
  heights = c(1, 2, 1)
)
ggarrange(F0, F4, nrow = 1, ncol = 2, legend = "right", common.legend = TRUE)
#################################################################

#### GO analysis Plots ####

#################################################################

## F0 -------------------------------------------------------------------------------

#### The first step is to create the saved results files with the abbreviated trinity ids ####
# OWA vs AM
res_OWAvsAM.df <- as.data.frame(res_OWAvsAM_F0)
res_OWAvsAM.df$fullID <- rownames(res_OWAvsAM.df)

parts <- strsplit(res_OWAvsAM.df$fullID, "::")

res_OWAvsAM.df$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OWAvsAM.df,
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/F0_OWAvsAM_results.csv",
  row.names = FALSE
)
# OW vs AM
res_OWvsAM.df <- as.data.frame(res_OWvsAM_F0)
res_OWvsAM.df$fullID <- rownames(res_OWvsAM.df)

parts <- strsplit(res_OWvsAM.df$fullID, "::")

res_OWvsAM.df$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OWvsAM.df,
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/F0_OWvsAM_results.csv",
  row.names = FALSE
)
# OA vs AM
res_OAvsAM.df <- as.data.frame(res_OAvsAM_F0)
res_OAvsAM.df$fullID <- rownames(res_OAvsAM.df)

parts <- strsplit(res_OAvsAM.df$fullID, "::")

res_OAvsAM.df$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OAvsAM.df,
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/F0_OAvsAM_results.csv",
  row.names = FALSE
)

# Prep to run the TopGo analysis
library(topGO)

mappingFile <- "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/mydata/trinotate_annotation_GOblastx_forTopGO.txt"

geneID2GO <- readMappings(
  file = mappingFile
)

cat("Genes in GO mapping:",
    length(geneID2GO),
    "\n")
# Genes in GO mapping: 86453 

#### Filter GO terms by size ####
allGO <- table(unlist(geneID2GO))

keepTerms <- names(allGO)[
  allGO >= 5 &
    allGO <= 500
]

geneID2GO.filtered <- lapply(
  geneID2GO,
  function(x) intersect(x, keepTerms)
)

geneID2GO.filtered <- geneID2GO.filtered[
  lengths(geneID2GO.filtered) > 0
]
length(geneID2GO.filtered)
# filtered list = 71320 
# Create a Function to Run a TopGo contrast
run_topGO_contrast <- function(
    infile,
    outfile,
    ontology = "BP",
    padj.cutoff = 0.05){
  
  deseq <- read.csv(
    infile,
    stringsAsFactors = FALSE
  )
  
  deseq <- subset(
    deseq,
    !is.na(padj)
  )
  
  geneList <- factor(
    as.integer(deseq$padj < padj.cutoff)
  )
  
  names(geneList) <- deseq$shortID
  
  cat("\nGenes tested:",
      length(geneList))
  
  cat("\nSignificant genes:",
      sum(geneList == 1),
      "\n")
  
  GOdata <- new(
    "topGOdata",
    ontology = ontology,
    allGenes = geneList,
    geneSelectionFun = function(x) x == 1,
    annot = annFUN.gene2GO,
    gene2GO = geneID2GO.filtered
  )
  
  resultWeight <- runTest(
    GOdata,
    algorithm = "weight01",
    statistic = "fisher"
  )
  
  GOresults <- GenTable(
    GOdata,
    weightFisher = resultWeight,
    orderBy = "weightFisher",
    topNodes = 100
  )
  
  write.csv(
    GOresults,
    outfile,
    row.names = FALSE
  )
  
  return(GOresults)
}

#### Run our TopGO function for the three different contrasts ####
GO_OA_F0 <- run_topGO_contrast(
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/F0_OAvsAM_results.csv",
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/GO_F0_OAvsAM_BP.csv"
)

GO_OW_F0 <- run_topGO_contrast(
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/F0_OWvsAM_results.csv",
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/GO_F0_OWvsAM_BP.csv"
)

GO_OWA_F0 <- run_topGO_contrast(
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/F0_OWAvsAM_results.csv",
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/GO_F0_OWAvsAM_BP.csv"
)

#### Make a bubble plot of the top 10 GO categories (OWA vs AM) ####
# Read TopGO results
go <- read.csv(
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/GO_F0_OWAvsAM_BP.csv",
  stringsAsFactors = FALSE
)

# Convert p-values to numeric
go$weightFisher <- gsub("^<\\s*", "", go$weightFisher)

go$weightFisher <- as.numeric(go$weightFisher)

# Convert counts to numeric
go$Significant <- as.numeric(go$Significant)

# Create -log10(p)
go$minusLogP <- -log10(go$weightFisher)

# Keep top 10 GO terms
go_top10 <- go %>%
  arrange(weightFisher) %>%
  slice(1:10)

# Order terms for plotting
go_top10$Term <- factor(
  go_top10$Term,
  levels = rev(go_top10$Term)
)

# Bubble plot
F0_OWA <- ggplot(
  go_top10,
  aes(
    x = minusLogP,
    y = Term
  )
) +
  geom_point(
    aes(
      size = Significant/Annotated,
      color = minusLogP
    )
  ) +
  scale_color_viridis_c() +
  theme_bw(base_size = 10) +
  labs(
    title = "Top GO Terms: OWA vs AM",
    x = expression(-logp),
    y = "GO Term",
    color = expression(-logp),
    size = "Proportion Significant\nGenes"
  )
F0_OWA
#### Make a bubble plot of the top 10 GO categories (OW vs AM) ####
# Read TopGO results
go <- read.csv(
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/GO_F0_OWvsAM_BP.csv",
  stringsAsFactors = FALSE
)

# Convert p-values to numeric
go$weightFisher <- gsub("^<\\s*", "", go$weightFisher)

go$weightFisher <- as.numeric(go$weightFisher)

# Convert counts to numeric
go$Significant <- as.numeric(go$Significant)

# Create -log10(p)
go$minusLogP <- -log10(go$weightFisher)

# Keep top 10 GO terms
go_top10 <- go %>%
  arrange(weightFisher) %>%
  slice(1:10)

# Order terms for plotting
go_top10$Term <- factor(
  go_top10$Term,
  levels = rev(go_top10$Term)
)

# Bubble plot
F0_OW <- ggplot(
  go_top10,
  aes(
    x = minusLogP,
    y = Term
  )
) +
  geom_point(
    aes(
      size = Significant/Annotated,
      color = minusLogP
    )
  ) +
  scale_color_viridis_c() +
  theme_bw(base_size = 10) +
  labs(
    title = "Top GO Terms: OW vs AM",
    x = expression(-logp),
    y = "GO Term",
    color = expression(-logp),
    size = "Proportion Significant\nGenes"
  )

#### Make a bubble plot of the top 10 GO categories (OWA vs AM) ####
# Read TopGO results
go <- read.csv(
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/GO_F0_OAvsAM_BP.csv",
  stringsAsFactors = FALSE
)

# Convert p-values to numeric
go$weightFisher <- gsub("^<\\s*", "", go$weightFisher)

go$weightFisher <- as.numeric(go$weightFisher)

# Convert counts to numeric
go$Significant <- as.numeric(go$Significant)

# Create -log10(p)
go$minusLogP <- -log10(go$weightFisher)

# Keep top 10 GO terms
go_top10 <- go %>%
  arrange(weightFisher) %>%
  slice(1:10)

# Order terms for plotting
go_top10$Term <- factor(
  go_top10$Term,
  levels = rev(go_top10$Term)
)

# Bubble plot
F0_OA <- ggplot(
  go_top10,
  aes(
    x = minusLogP,
    y = Term
  )
) +
  geom_point(
    aes(
      size = Significant/Annotated,
      color = minusLogP
    )
  ) +
  scale_color_viridis_c() +
  theme_bw(base_size = 14) +
  labs(
    title = "Top GO Terms: OA vs AM",
    x = expression(-logp),
    y = "GO Term",
    color = expression(-logp),
    size = "Proportion Significant\nGenes"
  )
F0_OA
## F4 -------------------------------------------------------------------------
#### The first step is to create the saved results files with the abbreviated trinity ids ####
# OWA vs AM
res_OWAvsAM.df <- as.data.frame(res_OWAvsAM_F4)
res_OWAvsAM.df$fullID <- rownames(res_OWAvsAM.df)

parts <- strsplit(res_OWAvsAM.df$fullID, "::")

res_OWAvsAM.df$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OWAvsAM.df,
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/F4_OWAvsAM_results.csv",
  row.names = FALSE
)
# OW vs AM
res_OWvsAM.df <- as.data.frame(res_OWvsAM_F4)
res_OWvsAM.df$fullID <- rownames(res_OWvsAM.df)

parts <- strsplit(res_OWvsAM.df$fullID, "::")

res_OWvsAM.df$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OWvsAM.df,
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/F4_OWvsAM_results.csv",
  row.names = FALSE
)
# OA vs AM
res_OAvsAM.df <- as.data.frame(res_OAvsAM_F4)
res_OAvsAM.df$fullID <- rownames(res_OAvsAM.df)

parts <- strsplit(res_OAvsAM.df$fullID, "::")

res_OAvsAM.df$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OAvsAM.df,
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/F4_OAvsAM_results.csv",
  row.names = FALSE
)

# Prep to run the TopGo analysis
library(topGO)

mappingFile <- "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/mydata/trinotate_annotation_GOblastx_forTopGO.txt"

geneID2GO <- readMappings(
  file = mappingFile
)

cat("Genes in GO mapping:",
    length(geneID2GO),
    "\n")
# Genes in GO mapping: 86453 

#### Filter GO terms by size ####
allGO <- table(unlist(geneID2GO))

keepTerms <- names(allGO)[
  allGO >= 5 &
    allGO <= 500
]

geneID2GO.filtered <- lapply(
  geneID2GO,
  function(x) intersect(x, keepTerms)
)

geneID2GO.filtered <- geneID2GO.filtered[
  lengths(geneID2GO.filtered) > 0
]
length(geneID2GO.filtered)
# filtered list = 71320 
# Create a Function to Run a TopGo contrast
run_topGO_contrast <- function(
    infile,
    outfile,
    ontology = "BP",
    padj.cutoff = 0.05){
  
  deseq <- read.csv(
    infile,
    stringsAsFactors = FALSE
  )
  
  deseq <- subset(
    deseq,
    !is.na(padj)
  )
  
  geneList <- factor(
    as.integer(deseq$padj < padj.cutoff)
  )
  
  names(geneList) <- deseq$shortID
  
  cat("\nGenes tested:",
      length(geneList))
  
  cat("\nSignificant genes:",
      sum(geneList == 1),
      "\n")
  
  GOdata <- new(
    "topGOdata",
    ontology = ontology,
    allGenes = geneList,
    geneSelectionFun = function(x) x == 1,
    annot = annFUN.gene2GO,
    gene2GO = geneID2GO.filtered
  )
  
  resultWeight <- runTest(
    GOdata,
    algorithm = "weight01",
    statistic = "fisher"
  )
  
  GOresults <- GenTable(
    GOdata,
    weightFisher = resultWeight,
    orderBy = "weightFisher",
    topNodes = 100
  )
  
  write.csv(
    GOresults,
    outfile,
    row.names = FALSE
  )
  
  return(GOresults)
}

#### Run our TopGO function for the three different contrasts ####
GO_OA_F4 <- run_topGO_contrast(
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/F4_OAvsAM_results.csv",
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/GO_F4_OAvsAM_BP.csv"
)

GO_OW_F4 <- run_topGO_contrast(
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/F4_OWvsAM_results.csv",
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/GO_F4_OWvsAM_BP.csv"
)

GO_OWA_F4 <- run_topGO_contrast(
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/F4_OWAvsAM_results.csv",
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/GO_F4_OWAvsAM_BP.csv"
)

#### Make a bubble plot of the top 10 GO categories (OWA vs AM) ####
# Read TopGO results
go <- read.csv(
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/GO_F4_OWAvsAM_BP.csv",
  stringsAsFactors = FALSE
)

# Convert p-values to numeric
go$weightFisher <- gsub("^<\\s*", "", go$weightFisher)

go$weightFisher <- as.numeric(go$weightFisher)

# Convert counts to numeric
go$Significant <- as.numeric(go$Significant)

# Create -log10(p)
go$minusLogP <- -log10(go$weightFisher)

# Keep top 10 GO terms
go_top10 <- go %>%
  arrange(weightFisher) %>%
  slice(1:10)

# Order terms for plotting
go_top10$Term <- factor(
  go_top10$Term,
  levels = rev(go_top10$Term)
)

# Bubble plot
F4_OWA <- ggplot(
  go_top10,
  aes(
    x = minusLogP,
    y = Term
  )
) +
  geom_point(
    aes(
      size = Significant/Annotated,
      color = minusLogP
    )
  ) +
  scale_color_viridis_c() +
  theme_bw(base_size = 10) +
  labs(
    title = "Top GO Terms: OWA vs AM",
    x = expression(-logp),
    y = "GO Term",
    color = expression(-logp),
    size = "Proportion Significant\nGenes"
  )
F4_OWA
#### Make a bubble plot of the top 10 GO categories (OW vs AM) ####
# Read TopGO results
go <- read.csv(
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/GO_F4_OWvsAM_BP.csv",
  stringsAsFactors = FALSE
)

# Convert p-values to numeric
go$weightFisher <- gsub("^<\\s*", "", go$weightFisher)

go$weightFisher <- as.numeric(go$weightFisher)

# Convert counts to numeric
go$Significant <- as.numeric(go$Significant)

# Create -log10(p)
go$minusLogP <- -log10(go$weightFisher)

# Keep top 10 GO terms
go_top10 <- go %>%
  arrange(weightFisher) %>%
  slice(1:10)

# Order terms for plotting
go_top10$Term <- factor(
  go_top10$Term,
  levels = rev(go_top10$Term)
)

# Bubble plot
F4_OW <- ggplot(
  go_top10,
  aes(
    x = minusLogP,
    y = Term
  )
) +
  geom_point(
    aes(
      size = Significant/Annotated,
      color = minusLogP
    )
  ) +
  scale_color_viridis_c() +
  theme_bw(base_size = 10) +
  labs(
    title = "Top GO Terms: OW vs AM",
    x = expression(-logp),
    y = "GO Term",
    color = expression(-logp),
    size = "Proportion Significant\nGenes"
  )
F4_OW
#### Make a bubble plot of the top 10 GO categories (OWA vs AM) ####
# Read TopGO results
go <- read.csv(
  "/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/myresults/GO_F4_OAvsAM_BP.csv",
  stringsAsFactors = FALSE
)

# Convert p-values to numeric
go$weightFisher <- gsub("^<\\s*", "", go$weightFisher)

go$weightFisher <- as.numeric(go$weightFisher)

# Convert counts to numeric
go$Significant <- as.numeric(go$Significant)

# Create -log10(p)
go$minusLogP <- -log10(go$weightFisher)

# Keep top 10 GO terms
go_top10 <- go %>%
  arrange(weightFisher) %>%
  slice(1:10)

# Order terms for plotting
go_top10$Term <- factor(
  go_top10$Term,
  levels = rev(go_top10$Term)
)

# Bubble plot
F4_OA <- ggplot(
  go_top10,
  aes(
    x = minusLogP,
    y = Term
  )
) +
  geom_point(
    aes(
      size = Significant/Annotated,
      color = minusLogP
    )
  ) +
  scale_color_viridis_c() +
  theme_bw(base_size = 14) +
  labs(
    title = "Top GO Terms: OA vs AM",
    x = expression(-logp),
    y = "GO Term",
    color = expression(-logp),
    size = "Proportion Significant\nGenes"
  )
F4_OA

## Combine all 6 into 1 figure? ##

ggarrange(F0_OWA, F0_OW, F0_OA, F4_OWA, F4_OW, F4_OA, nrow = 2, ncol = 3, legend = "right", common.legend = TRUE)


