library(tidyverse)

setwd("/users/a/s/asmit168/eco_genomics_2026/population_genomics/mydata/Genetic_Diversity")
###point to your results folder
OUTDIR="/users/a/s/asmit168/eco_genomics_2026/population_genomics/myresults/Selection"

env <- read_csv("sample_metadata.csv")
pw  <- read_tsv("fst/pairwise_fst_summary.txt")

# order the populations from south to north
env <- env %>%
  distinct(Population, Latitude, Longitude)
pop_order <- env$Population[order(env$Latitude)]

# add each pair in both directions so the heatmap is symmetrical
pw_both <- bind_rows(pw, select(pw, pop1 = pop2, pop2 = pop1, mean_fst, weighted_fst)) %>%
  mutate(pop1 = factor(pop1, levels = pop_order),
         pop2 = factor(pop2, levels = pop_order))

ggplot(pw_both, aes(x = pop1, y = pop2, fill = weighted_fst)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(weighted_fst, 3)), size = 2.5) +
  scale_fill_gradient(low = "white", high = "hotpink4", name = "Weighted\nFST") +
  labs(x = NULL, y = NULL,
       title = "Genome-wide pairwise FST (populations ordered south to north)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(OUTDIR,"/pairwise_fst_heatmap.png"), width = 7, height = 6)

# Does differentiation increase with the geographic distance between populations?

pw_dist <- pw %>%
  left_join(select(env, Population, lat1 = Latitude), by = c("pop1" = "Population")) %>%
  left_join(select(env, Population, lat2 = Latitude), by = c("pop2" = "Population")) %>%
  mutate(lat_diff = abs(lat1 - lat2))

ggplot(pw_dist, aes(x = lat_diff, y = weighted_fst)) +
  geom_point(size = 2, colour = "darkolivegreen3") +
  geom_smooth(method = "lm", se = FALSE, colour = "lightgray") +
  labs(x = "Difference in latitude (degrees)", y = "Weighted FST",
       title = "Isolation by distance?") +
  theme_classic()

ggsave(paste0(OUTDIR,"/fst_vs_distance.png"), width = 5, height = 4)

#### 3.5 Plotting the scan and finding outlier windows ####

##in case you closed R
OUTDIR="/users/a/s/asmit168/eco_genomics_2026/population_genomics/myresults/Selection"

pop1 <- "ASU"   #the same two populations as in the bash code
pop2 <- "KRU"

setwd("/users/a/s/asmit168/eco_genomics_2026/population_genomics/mydata/Genetic_Diversity/fst")
fst_win <- read_tsv(paste0(pop1, "_vs_", pop2, "_FST_10kb.windowed.weir.fst"))

# Helper: turn (scaffold, position) into a position along the whole genome by laying the scaffolds end to end, in the order in which they appear in the data.
add_genome_pos <- function(df) {
  scaffolds <- df %>% 
    group_by(CHROM) %>% ##group by CHROM
    summarise(length = max(POS)) %>% ## length = max position within a chromosome
    arrange(match(CHROM, unique(df$CHROM))) %>% ##make sure each is unique, arranged
    mutate(offset = cumsum(as.numeric(length)) - length, ###cumulatively number (helps with labels)
           shade  = factor(row_number() %% 2))    ##value for alternating colour 
  df %>%
    left_join(scaffolds, by = "CHROM") %>% ##join this data by given data with CHROM
    mutate(genome_pos = offset + POS) ###cumalative genome pos 
}

fst_win <- fst_win %>%
  mutate(POS = (BIN_START + BIN_END) / 2) %>% 
  add_genome_pos() ## this is our new function

# Now we can plot the scan and mark the windows in the top 1% of the FST distribution:

# empirical outlier threshold: the top 1% of windows
fst_cut <- quantile(fst_win$WEIGHTED_FST, probs = 0.99)
fst_cut

fst_win <- mutate(fst_win, outlier = WEIGHTED_FST >= fst_cut)

ggplot(fst_win, aes(x = genome_pos / 1e6, y = WEIGHTED_FST)) +
  geom_point(aes(colour = shade), size = 0.8) +
  geom_point(data = filter(fst_win, outlier), colour = "red", size = 1.2) +
  geom_hline(yintercept = fst_cut, linetype = "dashed") +
  scale_colour_manual(values = c("grey35", "grey70"), guide = "none") +
  labs(x = "Position along the genome (Mb; scaffolds end to end)",
       y = "Weighted FST (10-kb windows)",
       title = paste(pop1, "vs", pop2)) +
  theme_classic()

ggsave(paste0(OUTDIR,"/fst_scan_", pop1, "_vs_", pop2, ".png"), width = 10, height = 4)

# We can then save the outliers to a table.
fst_outliers <- fst_win %>%
  filter(outlier) %>%
  arrange(desc(WEIGHTED_FST)) %>%
  select(CHROM, BIN_START, BIN_END, N_VARIANTS, WEIGHTED_FST)

head(fst_outliers, 10)
write_tsv(fst_outliers, paste0(pop1, "_vs_", pop2, "_FST_10kb_outliers.txt"))

#### 3.6 Zooming in on the top window ####
# How do the individual SNPs behave in and around the most differentiated window? We can use the per-SNP FST file we made earlier:
top <- slice_max(fst_win, WEIGHTED_FST, n = 1, with_ties = FALSE)
top ##top window only

snp_fst <- read_tsv(paste0(pop1, "_vs_", pop2, ".weir.fst"),
                    col_types = "cid", na = c("", "NA", "nan", "-nan"))

zoom <- snp_fst %>%
  filter(CHROM == top$CHROM,
         POS > top$BIN_START - 100000,
         POS < top$BIN_END + 100000,
         !is.na(WEIR_AND_COCKERHAM_FST))   # drop SNPs that vary in neither population

ggplot(zoom, aes(x = POS / 1e6, y = WEIR_AND_COCKERHAM_FST)) +
  annotate("rect", xmin = top$BIN_START / 1e6, xmax = top$BIN_END / 1e6,
           ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.15) +
  geom_point(alpha = 0.6) +
  labs(x = paste(top$CHROM, "position (Mb)"), y = "Per-SNP FST",
       title = "Per-SNP FST around the top window (shaded)") +
  theme_classic()

ggsave(paste0(OUTDIR, "/fst_zoom_top_window.png"), width = 7, height = 4)

#### 4.5: BayPass’s view of population structure ####

library(tidyverse)

OUTDIR="/users/a/s/asmit168/eco_genomics_2026/population_genomics/myresults/Selection"

setwd("/users/a/s/asmit168/eco_genomics_2026/population_genomics/mydata/Genetic_Diversity")

pops <- readLines("Lottia.baypass.pops.txt")
npop <- length(pops)

env <- read_csv("pop_level_biooracle.csv")
pop_order<-env$pop[order(env$lat)]

omega <- as.matrix(read.table("Lottia_bp_mat_omega.out"))
dimnames(omega) <- list(pops, pops)

omega_cor <- cov2cor(omega)
round(omega_cor, 2)

omega_tree <- hclust(as.dist(1 - omega_cor), method = "average")

### we are not using ggplot so it is saved in a different way:
png(paste0(OUTDIR,"/baypass_omega_tree.png"), width = 1800, height = 1200, res = 300)
plot(omega_tree, main = "Populations clustered by BayPass Omega", xlab = "", sub = "")
dev.off()

omega_df <- as.data.frame(as.table(omega_cor)) %>%
  rename(pop1 = Var1, pop2 = Var2, correlation = Freq) %>%
  mutate(pop1 = factor(pop1, levels = pop_order),
         pop2 = factor(pop2, levels = pop_order))

###if you get any errors for visuals, try dev.off()
ggplot(omega_df, aes(x = pop1, y = pop2, fill = correlation)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(low = "cadetblue", mid = "white", high = "red4",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(x = NULL, y = NULL,
       title = "BayPass Omega as correlations (south to north)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(OUTDIR,"/baypass_omega.png"), width = 7, height = 6)

#### 4.6 XtX*: SNPs that are more differentiated than expected ####
snps <- read_tsv("Lottia.baypass.snps.txt")

snps$MRK<-rownames(snps)
snps$MRK<-as.numeric(snps$MRK)

npop <- 19

xtx <- read.table("Lottia_bp_summary_pi_xtx.out", header = TRUE, check.names = FALSE) %>%
  left_join(snps, by = "MRK") %>%
  mutate(lp_high   = pchisq(XtXst, npop, lower.tail = FALSE, log.p = TRUE),   #
         lp_low    = pchisq(XtXst, npop, lower.tail = TRUE,  log.p = TRUE),   
         direction = if_else(lp_high < lp_low, "high", "low"),
         log10_p   = -(pmin(lp_high, lp_low) + log(2)) / log(10),             
         pval      = 10^(-log10_p),
         qval      = p.adjust(pval, method = "BH"))
ggplot(xtx, aes(x = pval)) +
  geom_histogram(breaks = seq(0, 1, by = 0.025), fill = "plum", colour = "white") +
  labs(x = "p-value", y = "no. of SNPs") +
  theme_classic()

ggsave(paste0(OUTDIR,"/baypass_xtx_pvalues.png"), width = 5, height = 4)

##same helper tool as earlier with FST
add_genome_pos <- function(df) {
  scaffolds <- df %>% 
    group_by(CHROM) %>% ##group by CHROM
    summarise(length = max(POS)) %>% ## length = max position within a chromosome
    arrange(match(CHROM, unique(df$CHROM))) %>% ##make sure each is unique, arranged
    mutate(offset = cumsum(as.numeric(length)) - length, ###cumulatively number (helps with labels)
           shade  = factor(row_number() %% 2))    ##value for alternating colour 
  df %>%
    left_join(scaffolds, by = "CHROM") %>% ##join this data by given data with CHROM
    mutate(genome_pos = offset + POS) ###cumalative genome pos 
}


xtx <- xtx %>%
  add_genome_pos()     ##recall above if reloaded R

sum(xtx$outlier)                           

threshold <- quantile(xtx$qval, 0.001)   # top 0.1% (smallest q)

ggplot(xtx, aes(genome_pos / 1e6, -log10(qval))) +
  geom_point(aes(colour = shade), size = 0.4) +
  geom_hline(yintercept = -log10(threshold), linetype = "dashed", colour = "red") +
  scale_colour_manual(values = c("grey20", "grey80"), guide = "none") +
  labs(x = "Genome pos (Mb)", y = "-log10(q)") +
  theme_classic()

ggsave(paste0(OUTDIR,"/baypass_xtx_scan.png"), width = 10, height = 5)

xtx_outliers <- xtx %>%
  filter(qval < threshold) %>%
  arrange(desc(qval)) %>%
  select(MRK, CHROM, POS, M_XtX, XtXst, pval, qval)

head(xtx_outliers, 10)
write_tsv(xtx_outliers, "Lottia_bp_xtx_outliers.txt")

#### 4.7 Genotype–environment association ####
covariates <- readLines("Lottia.baypass.covariates.txt")

bf <- read.table("Lottia_bp_summary_betai_reg.out", header = TRUE, check.names = FALSE) %>%
  rename(BF_dB = `BF(dB)`) %>%
  mutate(covariate = covariates[COVARIABLE]) %>%
  left_join(snps, by = "MRK") %>%
  add_genome_pos()

bf %>%
  group_by(covariate) %>%
  summarise(n_decisive = sum(BF_dB > 20))

ggplot(bf, aes(x = genome_pos / 1e6, y = BF_dB)) +
  geom_point(aes(colour = shade), size = 0.6) +
  geom_point(data = filter(bf, BF_dB > 20), colour = "red", size = 1) +
  geom_hline(yintercept = 20, linetype = "dashed") +
  facet_wrap(~ covariate, ncol = 1) +
  scale_colour_manual(values = c("grey35", "grey70"), guide = "none") +
  labs(x = "Genome pos (Mb)",
       y = "Bayes factor (dB)") +
  theme_classic()

ggsave(paste0(OUTDIR,"/baypass_bf_scan.png"), width = 10, height = 6)

# What does an association actually look like? Let’s plot the allele frequencies of the SNP with the strongest association with sea surface temperature:
geno <- as.matrix(read.table("Lottia.baypass.geno"))

top_snp <- bf %>%
  filter(covariate == "thetao_mean") %>%
  slice_max(BF_dB, n = 1, with_ties = FALSE)
top_snp

counts <- geno[top_snp$MRK, ]
ref <- counts[seq(1, length(counts), by = 2)]
alt <- counts[seq(2, length(counts), by = 2)]

top_freq <- tibble(pop = pops, alt_freq = alt / (ref + alt)) %>%
  left_join(env, by = "pop")

ggplot(top_freq, aes(x = thetao_mean, y = alt_freq, label = pop)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 3) +
  ylim(0, 1) +
  labs(x = "Mean sea surface temperature (C)", y = "Allele Frequency",
       title = paste0("Top SST-associated SNP: ", top_snp$CHROM, ":", top_snp$POS)) +
  theme_classic()

ggsave(paste0(OUTDIR,"/baypass_top_sst_snp.png"), width = 5, height = 4)

##lets repeat with lat

top_snp <- bf %>%
  filter(covariate == "lat") %>%
  slice_max(BF_dB, n = 1, with_ties = FALSE)
top_snp

counts <- geno[top_snp$MRK, ]
ref <- counts[seq(1, length(counts), by = 2)]
alt <- counts[seq(2, length(counts), by = 2)]

top_freq <- tibble(pop = pops, alt_freq = alt / (ref + alt)) %>%
  left_join(env, by = "pop")

ggplot(top_freq, aes(x = lat, y = alt_freq, label = pop)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 3) +
  ylim(0, 1) +
  labs(x = "Latitude", y = "Allele Frequency",
       title = paste0("Top Latitude-associated SNP: ", top_snp$CHROM, ":", top_snp$POS)) +
  theme_classic()

ggsave(paste0(OUTDIR,"/baypass_top_lat_snp.png"), width = 5, height = 4)
