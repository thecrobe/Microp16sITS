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
# Plot Variance Explained
# =============================================================================
variance_explained_plot <- ggplot(permanova_results, aes(y = Variable, x = R2)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(
    title = "Variance Explained by PERMANOVA",
    x = "Variable",
    y = "Variance Explained (R²)",
    fill = "Dataset"
  ) +
  theme_classic() + 
  facet_grid(~Dataset) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlim(0,0.5)

ggsave(
  filename = file.path(save_dir, "Domestication_PERMANOVA_Variance_Explained_Barplot.pdf"),
  plot = variance_explained_plot,
  width = 8,
  height = 6,
  dpi = 300
)

variance_explained_plot

