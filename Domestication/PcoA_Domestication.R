# =============================================================================
# Load Required Libraries
# =============================================================================
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

set.seed(123)

save_dir <- "RF"

# =============================================================================
# Load and Clean Phyloseq Object (Bacteria)
# =============================================================================
load("/VU16s_roottraits_updated.Rdata")
VU16s_bacteria <- VU_16s_root_architecture_chem_dry
VU16s_bacteria <- prune_samples(sample_sums(VU16s_bacteria) > 0, VU16s_bacteria)

metadata_bacteria <- data.frame(sample_data(VU16s_bacteria)) %>%
  mutate(
    SRL = Total.Root.Length.mm / g_dry,
    RTD = Volume.mm3 / g_dry,
    Avg_Diameter = Average.Diameter.mm,
    C_N_Ratio = as.numeric(C_root_percent) / as.numeric(N__root_percent),
    N_root_percent = as.numeric(N__root_percent),
    FineRoots = Root.Length.mm0.0.5 + Root.Length.mm0.5.2.0,
    CoarseRoots = Root.Length.mm2.0.max
  ) %>% filter(C_N_Ratio > 0) %>% na.omit()

VU16s_bacteria <- prune_samples(rownames(metadata_bacteria), VU16s_bacteria)

# =============================================================================
# Load and Clean Phyloseq Object (Fungi)
# =============================================================================
load("/VU_ITS_roottraits_updated.Rdata")
VUITS_fungi <- VU_ITS_filtered_joined
VUITS_fungi <- prune_samples(sample_sums(VUITS_fungi) > 0, VUITS_fungi)

metadata_fungi <- data.frame(sample_data(VUITS_fungi)) %>%
  mutate(
    SRL = Total.Root.Length.mm / g_dry,
    RTD = Volume.mm3 / g_dry,
    Avg_Diameter = Average.Diameter.mm,
    C_N_Ratio = as.numeric(C_root_percent) / as.numeric(N__root_percent),
    N_root_percent = as.numeric(N__root_percent),
    FineRoots = Root.Length.mm0.0.5 + Root.Length.mm0.5.2.0,
    CoarseRoots = Root.Length.mm2.0.max
  ) %>% filter(C_N_Ratio > 0) %>% na.omit()

VUITS_fungi <- prune_samples(rownames(metadata_fungi), VUITS_fungi)

# =============================================================================
# Subset for Wild-Domestic Pairs
# =============================================================================
wild_domestic_pairs <- c("PS", "PW", "CR", "CA", "GS", "GM")

bacteria_wd <- subset_samples(VU16s_bacteria, PlantShort %in% wild_domestic_pairs)
bacteria_wd <- prune_taxa(taxa_sums(bacteria_wd) > 0, bacteria_wd)

fungi_wd <- subset_samples(VUITS_fungi, PlantShort %in% wild_domestic_pairs)
fungi_wd <- prune_taxa(taxa_sums(fungi_wd) > 0, fungi_wd)

metadata_bacteria_wd <- metadata_bacteria %>%
  filter(PlantShort %in% wild_domestic_pairs) %>%
  mutate(Pair = case_when(
    PlantShort %in% c("PS", "PW") ~ "Pair1",
    PlantShort %in% c("CR", "CA") ~ "Pair2",
    PlantShort %in% c("GS", "GM") ~ "Pair3"
  ))

metadata_fungi_wd <- metadata_fungi %>%
  filter(PlantShort %in% wild_domestic_pairs) %>%
  mutate(Pair = case_when(
    PlantShort %in% c("PS", "PW") ~ "Pair1",
    PlantShort %in% c("CR", "CA") ~ "Pair2",
    PlantShort %in% c("GS", "GM") ~ "Pair3"
  ))

# =============================================================================
# Remove NAs from model variables
# =============================================================================
metadata_bacteria_wd <- metadata_bacteria_wd %>%
  filter(!is.na(Domestication_status), !is.na(PlantShort), !is.na(Phosphorus), !is.na(Pair))

metadata_fungi_wd <- metadata_fungi_wd %>%
  filter(!is.na(Domestication_status), !is.na(PlantShort), !is.na(Phosphorus), !is.na(Pair))

# =============================================================================
# Match metadata to OTU tables
# =============================================================================
metadata_bacteria_wd <- metadata_bacteria_wd[rownames(metadata_bacteria_wd) %in% sample_names(bacteria_wd), ]
bacteria_wd <- prune_samples(rownames(metadata_bacteria_wd), bacteria_wd)

metadata_fungi_wd <- metadata_fungi_wd[rownames(metadata_fungi_wd) %in% sample_names(fungi_wd), ]
fungi_wd <- prune_samples(rownames(metadata_fungi_wd), fungi_wd)

# =============================================================================
# Calculate Bray-Curtis Distances
# =============================================================================
otu_bacteria_wd <- as.matrix(otu_table(bacteria_wd))
if (taxa_are_rows(bacteria_wd)) otu_bacteria_wd <- t(otu_bacteria_wd)
bray_bacteria_wd <- vegdist(otu_bacteria_wd, method = "bray")

otu_fungi_wd <- as.matrix(otu_table(fungi_wd))
if (taxa_are_rows(fungi_wd)) otu_fungi_wd <- t(otu_fungi_wd)
bray_fungi_wd <- vegdist(otu_fungi_wd, method = "bray")

# =============================================================================
# PERMANOVA with Strata by Pair (Bacteria)
# =============================================================================
permanova_bacteria_wd <- adonis2(
  bray_bacteria_wd ~ Domestication_status + PlantShort + Phosphorus,
  data = metadata_bacteria_wd,
  strata = metadata_bacteria_wd$Pair,
  permutations = 999,
  by = "terms"
)
print(permanova_bacteria_wd)

# =============================================================================
# Beta Dispersion Tests (Bacteria)
# =============================================================================
for (factor in c("Domestication_status", "PlantShort", "Phosphorus")) {
  group <- metadata_bacteria_wd[[factor]]
  bd <- betadisper(bray_bacteria_wd, group = group)
  cat("\nBacteria -", factor, "\n")
  print(anova(bd))
  print(permutest(bd, permutations = 999))
}

# =============================================================================
# PERMANOVA with Strata by Pair (Fungi)
# =============================================================================
permanova_fungi_wd <- adonis2(
  bray_fungi_wd ~ Domestication_status + PlantShort + Phosphorus,
  data = metadata_fungi_wd,
  strata = metadata_fungi_wd$Pair,
  permutations = 999,
  by = "terms"
)
print(permanova_fungi_wd)

# =============================================================================
# Beta Dispersion Tests (Fungi)
# =============================================================================
for (factor in c("Domestication_status", "PlantShort", "Phosphorus")) {
  group <- metadata_fungi_wd[[factor]]
  bd <- betadisper(bray_fungi_wd, group = group)
  cat("\nFungi -", factor, "\n")
  print(anova(bd))
  print(permutest(bd, permutations = 999))
}

# =============================================================================
# Create PERMANOVA R² Summary Data (rounded)
# =============================================================================
permanova_results <- data.frame(
  Variable = rep(c("Domestication_status", "PlantShort", "Phosphorus"), 2),
  R2 = c(0.10, 0.33, 0.02, 0.03, 0.40, 0.01),
  Dataset = rep(c("Bacteria", "Fungi"), each = 3)
)
# =============================================================================
# Create Variance Explained Plot with Dodged Bars (No Pattern)
# =============================================================================
variance_explained_plot <- ggplot(permanova_results, aes(y = Variable, x = R2, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.8) +
  scale_fill_manual(values = c("Bacteria" = "lightgrey", "Fungi" = "darkgrey")) +
  labs(
    title = "Variance Explained by PERMANOVA",
    x = "Variance Explained (R²)",
    y = "Variable",
    fill = "Dataset"
  ) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlim(0, 0.5)

# Save the plot
ggsave(
  filename = file.path(save_dir, "Domestication_PERMANOVA_Variance_Explained_Barplot.pdf"),
  plot = variance_explained_plot,
  width = 8,
  height = 6,
  dpi = 300,
  scale = 0.5
)

variance_explained_plot

# =============================================================================
# Define Custom Colors
# =============================================================================
custom_colors <- c(
  "PS" = "#7A3DA5", "PW" = "#7A3DA5", 
  "CR" = "#2A9D8F", "CA" = "#2A9D8F", 
  "GS" = "#FFC75F", "GM" = "#FFC75F"
)

# =============================================================================
# Compute PCoAs
# =============================================================================

# Bacteria PCoA Calculation (using cmdscale for Bray-Curtis distances)
pcoa_bacteria <- cmdscale(bray_bacteria_wd, eig = TRUE, k = 2)
pcoa_bacteria_scores <- as.data.frame(pcoa_bacteria$points) %>%
  tibble::rownames_to_column("SampleID") %>%
  left_join(metadata_bacteria_wd %>% mutate(SampleID = rownames(metadata_bacteria_wd)), by = "SampleID")

# Fungi PCoA Calculation (using cmdscale for Bray-Curtis distances)
pcoa_fungi <- cmdscale(bray_fungi_wd, eig = TRUE, k = 2)
pcoa_fungi_scores <- as.data.frame(pcoa_fungi$points) %>%
  tibble::rownames_to_column("SampleID") %>%
  left_join(metadata_fungi_wd %>% mutate(SampleID = rownames(metadata_fungi_wd)), by = "SampleID")

# =============================================================================
# Compute Variance Explained for the First Two PCoA Axes
# =============================================================================

# Bacteria - Variance Explained for the First Two Axes
eigenvalues_bacteria <- pcoa_bacteria$eig
variance_explained_bacteria <- eigenvalues_bacteria / sum(eigenvalues_bacteria) * 100
variance_explained_bacteria_first_two <- variance_explained_bacteria[1:2]

# Fungi - Variance Explained for the First Two Axes
eigenvalues_fungi <- pcoa_fungi$eig
variance_explained_fungi <- eigenvalues_fungi / sum(eigenvalues_fungi) * 100
variance_explained_fungi_first_two <- variance_explained_fungi[1:2]

# Print the variance explained for both bacteria and fungi (first two axes)
cat("Variance Explained for Bacteria (First Two Axes): ", variance_explained_bacteria_first_two, "%\n")
cat("Variance Explained for Fungi (First Two Axes): ", variance_explained_fungi_first_two, "%\n")


# =============================================================================
# Extract PCoA Scores Directly with Metadata
# =============================================================================

# Bacteria
pcoa_bacteria_scores <- as.data.frame(pcoa_bacteria$points) %>%
  tibble::rownames_to_column("SampleID") %>%
  left_join(metadata_bacteria_wd %>% mutate(SampleID = rownames(metadata_bacteria_wd)), by = "SampleID")

# Fungi
pcoa_fungi_scores <- as.data.frame(pcoa_fungi$points) %>%
  tibble::rownames_to_column("SampleID") %>%
  left_join(metadata_fungi_wd %>% mutate(SampleID = rownames(metadata_fungi_wd)), by = "SampleID")


# =============================================================================
# PCoA Plots
# =============================================================================
pcoa_plot_bacteria<-ggplot(na.omit(pcoa_bacteria_scores), aes(x = V1, y = V2, color = PlantShort, shape = Domestication_status)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = c("Wild" = 17, "Domestic" = 16)) +
  labs(
    title = "Bacteria PCoA (Bray-Curtis)",
    x = "PCoA1",
    y = "PCoA2"
  ) +
  theme_classic()+ 
  stat_ellipse()

pcoa_plot_fungi <- ggplot(na.omit(pcoa_fungi_scores), aes(x = V1, y = V2, color = PlantShort, shape = Domestication_status)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = c("Wild" = 17, "Domestic" = 16)) +
  labs(
    title = "Bacteria PCoA (Bray-Curtis)",
    x = "PCoA1",
    y = "PCoA2"
  ) +
  theme_classic()+ 
  stat_ellipse()

# =============================================================================
# Combine Variance Explained Plot with PCoAs and Use Common Legend
# =============================================================================
final_plot <- (pcoa_plot_bacteria | pcoa_plot_fungi) + plot_layout(guides = "collect")

# Save the final plot with a common legend
ggsave(
  filename = file.path(save_dir, "Domestication_BactFungi_PCoAs.pdf"),
  plot = final_plot,
  dpi = 300, 
  scale = 0.4
)

final_plot
