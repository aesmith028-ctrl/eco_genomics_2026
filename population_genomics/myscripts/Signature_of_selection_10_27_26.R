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
  scale_fill_gradient(low = "white", high = "red4", name = "Weighted\nFST") +
  labs(x = NULL, y = NULL,
       title = "Genome-wide pairwise FST (populations ordered south to north)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(OUTDIR,"/pairwise_fst_heatmap.png"), width = 7, height = 6)