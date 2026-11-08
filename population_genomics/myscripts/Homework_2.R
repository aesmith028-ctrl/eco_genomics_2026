library(tidyverse)

setwd("~/eco_genomics_2026/population_genomics/mydata")   # the folder with your plink output

# Now read in the PCA results from plink:

# the .eigenvec file has no column names, so col_names = FALSE (they become X1, X2, ...)
pca <- read_table("Lg_LD_0.5.eigenvec", col_names = FALSE)

# the .eigenval file is just a list of numbers, so scan() reads it into a vector
eigenval <- scan("Lg_LD_0.5.eigenval")

# columns 1 and 2 are the family and individual IDs, which are identical
# thanks to --double-id, so drop column 2 (remember: [rows, columns])
pca <- pca[, -2]

# rename the columns: the first is the sample ID, the rest are PC1, PC2, ...
# paste0() joins text and numbers, e.g. paste0("PC", 1:3) gives "PC1" "PC2" "PC3"
names(pca)[1] <- "ID"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca) - 1))

head(pca)

#### How much variation does each PC explain? ####
# First let’s see how much of the variation each PC explains. This helps us decide how many PCs are worth looking at.

# convert the eigenvalues to a percentage variance explained (pve)
pve <- data.frame(PC = seq_along(eigenval),
                  pve = eigenval / sum(eigenval) * 100)

# make a folder for our figures

pve_plot <- ggplot(pve, aes(PC, pve)) +
  geom_col() +
  labs(x = "Principal component", y = "Percentage variance explained") +
  geom_bar(fill="plum3", stat="identity") +
  theme_bw()
pve_plot

OUTDIR="~/eco_genomics_2026/population_genomics/myresults/Population_Structure"

ggsave(pve_plot, file=paste0(OUTDIR, "/variance_0.5.png"), width = 6, height = 4)

# This is called a scree plot. Look for the point where the bars level off: PCs beyond that point are mostly noise. In genomic data each PC usually explains only a small percentage of the variation, because most variation across the genome is neutral and differs between individuals rather than between populations.

#### PC1 vs PC2 ####
# To colour the PCA by site we need to join it to our sample metadata. Joining matches the rows of the two tables using the ID column, so the IDs in sample_metadata.csv must match the sample names in the VCF exactly.
setwd("~/eco_genomics_2026/population_genomics/mydata/Homework_2")
meta <- read_csv("sample_metadata.csv")
pca <- left_join(pca, meta, by = "ID")

# check every individual matched a row in the metadata: this should print 0
sum(is.na(pca$Site))

# order the sites from north to south, so the legend follows the coast
pca <- pca %>% mutate(Site = fct_reorder(Site, -Latitude))

# axis labels with the % variance explained, e.g. "PC1 (4.2%)"
pc1_lab <- paste0("PC1 (", round(pve$pve[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(pve$pve[2], 1), "%)")

pca_plot_1 <- ggplot(pca, aes(PC1, PC2, colour = Site)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_colour_viridis_d(end = 0.9) +
  labs(x = pc1_lab, y = pc2_lab, colour = "Site (north to south)") +
  theme_bw()
pca_plot_1

OUTDIR="~/eco_genomics_2026/population_genomics/myresults/Population_Structure"

ggsave(pca_plot_1, file = paste0(OUTDIR, "/Lg_0.5_pca.png"), width = 8, height = 6)

#### PC2 vs PC3 ####
# To colour the PCA by site we need to join it to our sample metadata. Joining matches the rows of the two tables using the ID column, so the IDs in sample_metadata.csv must match the sample names in the VCF exactly.
setwd("~/eco_genomics_2026/population_genomics/mydata/Homework_2")
meta <- read_csv("sample_metadata.csv")
pca <- left_join(pca, meta, by = "ID")

# check every individual matched a row in the metadata: this should print 0
sum(is.na(pca$Site))

# order the sites from north to south, so the legend follows the coast
pca <- pca %>% mutate(Site = fct_reorder(Site, -Latitude))

# axis labels with the % variance explained, e.g. "PC1 (4.2%)"
pc2_lab <- paste0("PC2 (", round(pve$pve[2], 1), "%)") #changed number in brackets to reflect pc2 not pc1
pc3_lab <- paste0("PC3 (", round(pve$pve[3], 1), "%)") #changed number in brackets to reflect pc3 not pc2

pca_plot_2 <- ggplot(pca, aes(PC2, PC3, colour = Site)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_colour_viridis_d(end = 0.9) +
  labs(x = pc2_lab, y = pc3_lab, colour = "Site (north to south)") +
  theme_bw()
pca_plot_2

OUTDIR="~/eco_genomics_2026/population_genomics/myresults/Population_Structure"

ggsave(pca_plot_2, file = paste0(OUTDIR, "/Lg_0.5_pca_2.png"), width = 8, height = 6)

#### Is there a pattern with latitude? ####
# To plot PC1, PC2, and PC3 against latitude side by side, we first reshape the data with pivot_longer(). Instead of one column per PC, we get one row per individual per PC, which lets facet_wrap() give each PC its own panel:
lat_plot_HW <- pca %>%
  select(ID, Site, Latitude, PC1, PC2, PC3) %>% ##Select columnns to use
  pivot_longer(c(PC1, PC2, PC3), names_to = "PC", values_to = "score") %>%
  ggplot(aes(Latitude, score, colour = Site)) +
  geom_point(size = 2, alpha = 0.6) +
  facet_wrap(~PC, ncol = 1, scales = "free_y") + ##plots PC1 and PC2 seperately
  scale_colour_viridis_d(end = 0.9) + 
  labs(x = "Latitude (°N)", y = "PC score", colour = "Site (north to south)") +
  theme_bw()+
  coord_flip()
lat_plot_HW 
OUTDIR="~/eco_genomics_2026/population_genomics/myresults/Population_Structure"

ggsave(lat_plot_HW, file = paste0(OUTDIR, "/Lg_0.5_pca_latitude.png"), width = 8, height = 6)

# **********************************************************************************
#### ADMIXTURE ####

#### Choosing K ####

library(tidyverse)

setwd("~/eco_genomics_2026/population_genomics/mydata/Homework_2")   # the folder with your plink output

cv <- read_tsv("cv_errors.tsv", col_names = c("K", "CV_error"))

cv_plot <- ggplot(cv, aes(K, CV_error)) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = cv$K) +
  labs(x = "K (number of ancestral populations)", y = "Cross-validation error") +
  theme_bw()
cv_plot

OUTDIR="~/eco_genomics_2026/population_genomics/myresults/Population_Structure"

ggsave(cv_plot, file = paste0(OUTDIR, "/cv_error.png"), width = 6, height = 4)

#### Plotting Ancestry Proportions ####
# The .Q files have no sample IDs: row 1 is the first individual in the .fam file, row 2 the second, and so on. So we take the IDs from column 2 of the .fam file. We then read the results for K = 1 to 5 and stack them into one table, so that we can plot them together:
fam  <- read_table("Lg_LD_0.5.fam", col_names = FALSE)
meta <- read_csv("sample_metadata.csv")

# a function that reads the .Q file for one value of K and reshapes it
# into one row per individual per cluster
read_q <- function(K) {
  q <- read_table(paste0("Lg_LD_0.5.", K, ".Q"), col_names = FALSE)
  names(q) <- paste0("Cluster", 1:K)
  q %>%
    mutate(ID = fam$X2, K = K) %>%
    pivot_longer(starts_with("Cluster"), names_to = "Cluster", values_to = "Ancestry")
}

# map_dfr() runs read_q() for each K and binds the results together row-wise
# (change 2:5 to look at other values of K)
q_all <- map_dfr(1:5, read_q) %>%
  left_join(meta, by = "ID") %>%
  # order sites north to south, as in the PCA plots
  mutate(Site = fct_reorder(Site, -Latitude),
         K = paste0("K = ", K))

# check every individual matched the metadata: this should print 0
sum(is.na(q_all$Site))

admix_plot <- ggplot(q_all, aes(ID, Ancestry, fill = Cluster)) +
  geom_col(width = 1) +
  facet_grid(K ~ Site, scales = "free_x", space = "free_x") +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1)) +
  labs(x = NULL, y = "Ancestry proportion") +
  theme_minimal() +
  theme(axis.text.x     = element_blank(),         # too many sample names to read
        panel.grid      = element_blank(),
        panel.spacing.x = unit(0.1, "lines"),      # thin gap between sites
        panel.spacing.y = unit(0.8, "lines"),      # gap between values of K
        strip.text.x    = element_text(angle = 90, hjust = 0),
        legend.position = "none")
admix_plot
OUTDIR="~/eco_genomics_2026/population_genomics/myresults/Population_Structure"

ggsave(admix_plot, file = paste0(OUTDIR, "/admixture_K1-5.png"), width = 12, height = 7)

