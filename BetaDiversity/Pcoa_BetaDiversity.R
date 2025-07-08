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
    FineRoots = Root.Length.mm0.0.5 + Root.Length.mm0.5.2.0,
    CoarseRoots = Root.Length.mm2.0.max
  ) %>% na.omit()

VU16s_bacteria <- prune_samples(rownames(metadata_bacteria), VU16s_bacteria)

# =============================================================================
# Load and Clean Phyloseq Object (Fungi)
# =============================================================================
load("VU_ITS_roottraits_updated.Rdata")
VUITS_fungi <- VU_ITS_filtered_joined
VUITS_fungi <- prune_samples(sample_sums(VUITS_fungi) > 0, VUITS_fungi)

metadata_fungi <- data.frame(sample_data(VUITS_fungi)) %>%
  mutate(
    SRL = Total.Root.Length.mm / g_dry,
    RTD = Volume.mm3 / g_dry,
    Avg_Diameter = Average.Diameter.mm,
    C_N_Ratio = as.numeric(C_root_percent) / as.numeric(N__root_percent),
    FineRoots = Root.Length.mm0.0.5 + Root.Length.mm0.5.2.0,
    CoarseRoots = Root.Length.mm2.0.max
  ) %>% na.omit()

VUITS_fungi <- prune_samples(rownames(metadata_fungi), VUITS_fungi)

# =============================================================================
# PCoA Ordination (Bacteria)
# =============================================================================
otu_bacteria <- as.matrix(otu_table(VU16s_bacteria))
if (taxa_are_rows(VU16s_bacteria)) otu_bacteria <- t(otu_bacteria)
bray_bacteria <- vegdist(otu_bacteria, method = "bray")
pcoa_bacteria <- cmdscale(bray_bacteria, eig = TRUE, k = 2)

ordination_bacteria <- as.data.frame(pcoa_bacteria$points)
colnames(ordination_bacteria) <- c("PCoA1", "PCoA2")
ordination_bacteria$PlantShort <- metadata_bacteria$PlantShort

var_explained_bacteria <- pcoa_bacteria$eig / sum(pcoa_bacteria$eig) * 100

centroids_bacteria <- ordination_bacteria %>%
  group_by(PlantShort) %>%
  summarise(
    Mean_PCoA1 = mean(PCoA1, na.rm = TRUE),
    Mean_PCoA2 = mean(PCoA2, na.rm = TRUE)
  )

# =============================================================================
# PCoA Ordination (Fungi)
# =============================================================================
otu_fungi <- as.matrix(otu_table(VUITS_fungi))
if (taxa_are_rows(VUITS_fungi)) otu_fungi <- t(otu_fungi)
bray_fungi <- vegdist(otu_fungi, method = "bray")
pcoa_fungi <- cmdscale(bray_fungi, eig = TRUE, k = 2)

ordination_fungi <- as.data.frame(pcoa_fungi$points)
colnames(ordination_fungi) <- c("PCoA1", "PCoA2")
ordination_fungi$PlantShort <- metadata_fungi$PlantShort

var_explained_fungi <- pcoa_fungi$eig / sum(pcoa_fungi$eig) * 100

centroids_fungi <- ordination_fungi %>%
  group_by(PlantShort) %>%
  summarise(
    Mean_PCoA1 = mean(PCoA1, na.rm = TRUE),
    Mean_PCoA2 = mean(PCoA2, na.rm = TRUE)
  )

# =============================================================================
# Define Colors
# =============================================================================
custom_colors <- c(
  "BD" = "#168AAD", "ZM" = "#1FAF9D",  
  "AT" = "#264653", "BO" = "#2E5A72",  
  "CP" = "#2D6A4F", "CS" = "#3E8C6E",  
  "AH" = "#A52A2A", "LA" = "#C44F41",  
  "VU" = "#9A031E", "PV" = "#C2381C", "PC" = "#E86A33",  
  "GM" = "#EEA04B", "GS" = "#FFC75F",  
  "CA" = "#3A86FF", "CR" = "#2A9D8F",  
  "PW" = "#5A189A", "PS" = "#7A3DA5", "LC" = "#A070D6", "VS" = "#C9A2F4",  
  "MT" = "#F8C300", "MS" = "#F7D564", "SL" = "#005F99"
)

# =============================================================================
# Create Bacteria Plot
# =============================================================================
bacteria_pcoa_plot <- ggplot() +
  geom_point(data = ordination_bacteria, 
             aes(x = PCoA1, y = PCoA2, color = PlantShort), 
             size = 2, alpha = 0.4) +
  geom_point(data = centroids_bacteria, 
             aes(x = Mean_PCoA1, y = Mean_PCoA2, color = PlantShort), 
             size = 5, shape = 21, stroke = 1.5) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Bacteria PCoA Ordination",
    x = paste0("PCoA1 (", round(var_explained_bacteria[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_explained_bacteria[2], 1), "%)")
  ) +
  theme_classic()

# =============================================================================
# Create Fungi Plot
# =============================================================================
fungi_pcoa_plot <- ggplot() +
  geom_point(data = ordination_fungi, 
             aes(x = PCoA1, y = PCoA2, color = PlantShort), 
             size = 2, alpha = 0.4) +
  geom_point(data = centroids_fungi, 
             aes(x = Mean_PCoA1, y = Mean_PCoA2, color = PlantShort), 
             size = 5, shape = 21, stroke = 1.5) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Fungi PCoA Ordination",
    x = paste0("PCoA1 (", round(var_explained_fungi[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_explained_fungi[2], 1), "%)")
  ) +
  theme_classic()

# =============================================================================
# Combine and Save Multipanel Plot
# =============================================================================
multi_panel_pcoa <- bacteria_pcoa_plot / fungi_pcoa_plot +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(
  filename = file.path(save_dir, "PCoA_Ordination_Bacteria_Fungi.pdf"),
  plot = multi_panel_pcoa,
  width = 10,
  height = 10,
  dpi = 300,scale = 0.6,
  device = "pdf"
)


# =============================================================================
# PERMANOVA for Bacteria
# =============================================================================
permanova_bacteria <- adonis2(
  bray_bacteria ~ PlantShort + Phosphorus,
  data = metadata_bacteria,
  strata = metadata_bacteria$Greenhouse_Block,
  permutations = 999,
  by = "terms"
)
print(permanova_bacteria)

# =============================================================================
# PERMANOVA for Fungi
# =============================================================================
permanova_fungi <- adonis2(
  bray_fungi ~ PlantShort + Phosphorus,
  data = metadata_fungi,
  strata = metadata_fungi$Greenhouse_Block,
  permutations = 999,
  by = "terms"
)
print(permanova_fungi)


# =============================================================================
# Beta Dispersion Test for Bacteria
# =============================================================================
group_bacteria <- metadata_bacteria$PlantShort
betadisper_bacteria <- betadisper(bray_bacteria, group = group_bacteria)
anova_betadisper_bacteria <- anova(betadisper_bacteria)
permutest_betadisper_bacteria <- permutest(betadisper_bacteria, permutations = 999)
print(anova_betadisper_bacteria)
print(permutest_betadisper_bacteria)

# =============================================================================
# Beta Dispersion Test for Fungi
# =============================================================================
group_fungi <- metadata_fungi$PlantShort
betadisper_fungi <- betadisper(bray_fungi, group = group_fungi)
anova_betadisper_fungi <- anova(betadisper_fungi)
permutest_betadisper_fungi <- permutest(betadisper_fungi, permutations = 999)
print(anova_betadisper_fungi)
print(permutest_betadisper_fungi)
