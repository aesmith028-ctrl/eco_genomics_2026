#### Preparing the Input ####

# STOP!! Do not run this on the full VCF! ADMIXTURE takes a long time to run even with batchjobs, so I have already done this for you! If you would like to follow along you can try this on your randomly subsetted data and only do K = 1-3

#### Choosing K ####

library(tidyverse)

setwd("~/eco_genomics_2026/population_genomics/mydata/Population_Structure")   # the folder with your plink output

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
fam  <- read_table("Lg_LD_0.2.fam", col_names = FALSE)
meta <- read_csv("sample_metadata.csv")

# a function that reads the .Q file for one value of K and reshapes it
# into one row per individual per cluster
read_q <- function(K) {
  q <- read_table(paste0("Lg_LD_0.2.", K, ".Q"), col_names = FALSE)
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

