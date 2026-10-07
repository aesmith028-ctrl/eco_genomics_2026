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

# Perform DESeq2 analysis on the subset
dds_F0 <- DESeq(dds_F0)

dds_F4 <- DESeq(dds_F4)

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

library(eulerr)
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
plot(fit1, main = "F0 DEGS for each Treatment", lty = 1:3, quantities = TRUE,
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
plot(fit1, main = "F4 DEGS for each Treatment", lty = 1:3, quantities = TRUE,
     fills = c("lightpink3", "lightblue3", "thistle3"))
# lty changes the lines

#cross check with above lengths of DEGS: the four values, unique, shared with one other, shared with the second other, shared across all treatments, should sum to the length of DEGs for each treatment contrast to AM
2668+2405+338+106 # 5517 total OW
1133+2405+338+42  # 3918 total OWA
116+42+106+338    # 602  total OA
