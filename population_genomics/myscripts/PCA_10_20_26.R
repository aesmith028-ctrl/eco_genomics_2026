library(tidyverse)

setwd("~/eco_genomics_2026/population_genomics/mydata")   # the folder with your plink output

# Now read in the PCA results from plink:

# the .eigenvec file has no column names, so col_names = FALSE (they become X1, X2, ...)
pca <- read_table("Lg_LD_0.2.eigenvec", col_names = FALSE)

# the .eigenval file is just a list of numbers, so scan() reads it into a vector
eigenval <- scan("Lg_LD_0.2.eigenval")

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
  theme_bw()
pve_plot

OUTDIR="~/eco_genomics_2026/population_genomics/myresults/Population_Structure"

ggsave(pve_plot, file=paste0(OUTDIR, "/variance_0.2.png"), width = 6, height = 4)

    # This is called a scree plot. Look for the point where the bars level off: PCs beyond that point are mostly noise. In genomic data each PC usually explains only a small percentage of the variation, because most variation across the genome is neutral and differs between individuals rather than between populations.
#### The PCA Plot ####
    # To colour the PCA by site we need to join it to our sample metadata. Joining matches the rows of the two tables using the ID column, so the IDs in sample_metadata.csv must match the sample names in the VCF exactly.
setwd("~/eco_genomics_2026/population_genomics/mydata/Population_Structure")
meta <- read_csv("sample_metadata.csv")
pca <- left_join(pca, meta, by = "ID")

# check every individual matched a row in the metadata: this should print 0
sum(is.na(pca$Site))

# order the sites from north to south, so the legend follows the coast
pca <- pca %>% mutate(Site = fct_reorder(Site, -Latitude))

# axis labels with the % variance explained, e.g. "PC1 (4.2%)"
pc1_lab <- paste0("PC1 (", round(pve$pve[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(pve$pve[2], 1), "%)")

pca_plot <- ggplot(pca, aes(PC1, PC2, colour = Site)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_colour_viridis_d(end = 0.9) +
  labs(x = pc1_lab, y = pc2_lab, colour = "Site (north to south)") +
  theme_bw()
pca_plot

OUTDIR="~/eco_genomics_2026/population_genomics/myresults/Population_Structure"

ggsave(pca_plot, file = paste0(OUTDIR, "/Lg_0.2_pca.png"), width = 8, height = 6)

#### Is there a pattern with latitude? ####
    # To plot PC1 and PC2 against latitude side by side, we first reshape the data with pivot_longer(). Instead of one column per PC, we get one row per individual per PC, which lets facet_wrap() give each PC its own panel:
lat_plot <- pca %>%
  select(ID, Site, Latitude, PC1, PC2) %>% ##Select columnns to use
  pivot_longer(c(PC1, PC2), names_to = "PC", values_to = "score") %>%
  ggplot(aes(Latitude, score, colour = Site)) +
  geom_point(size = 2, alpha = 0.6) +
  facet_wrap(~PC, ncol = 1, scales = "free_y") + ##plots PC1 and PC2 seperately
  scale_colour_viridis_d(end = 0.9) + 
  labs(x = "Latitude (°N)", y = "PC score", colour = "Site (north to south)") +
  theme_bw()+
  coord_flip()
lat_plot 
OUTDIR="~/eco_genomics_2026/population_genomics/myresults/Population_Structure"

ggsave(lat_plot, file = paste0(OUTDIR, "/Lg_0.2_pca_latitude.png"), width = 8, height = 6)

#### Bonus: a map of our sites ####
    # We can also map our sites and put the map next to the PCA (might be useful for your project!). map_data() needs the maps package to be installed (we don’t need to load it with library()), and patchwork lets us combine ggplots. If either isn’t installed, run install.packages(c("maps", "patchwork")) first.

library(patchwork) ##combined plots

world <- map_data("world")
sites <- pca %>% distinct(Site, Latitude, Longitude)

map_plot <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group),
               fill = "grey90", colour = "grey60") +
  geom_point(data = sites, aes(Longitude, Latitude, colour = Site), size = 3) +
  scale_colour_viridis_d(end = 0.9) +
  # zoom to our sites, with a 2-degree margin all round
  coord_quickmap(xlim = range(sites$Longitude) + c(-2, 2),
                 ylim = range(sites$Latitude) + c(-2, 2)) +
  labs(x = "Longitude", y = "Latitude", colour = "Site (north to south)") +
  theme_bw()
map_plot

# with patchwork, + puts plots side by side
# we hide the map's legend, as the PCA's legend already shows the sites
map_pca <- map_plot + theme(legend.position = "none") + pca_plot
map_pca

OUTDIR="~/eco_genomics_2026/population_genomics/myresults/Population_Structure"

ggsave(map_pca, file = paste0(OUTDIR, "/Lg_0.2_pca_map.png"), width = 8, height = 6)

