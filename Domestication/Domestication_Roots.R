# =============================================================================
# Load Required Libraries
# =============================================================================
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(brms)
library(lme4)
library(lmerTest)
library(emmeans)

set.seed(123)

save_dir <- "/Users/justinstewart/Dropbox/Collaborations/FamilyExperiment/Data/Code/New/Feb2024/Figures/"

# =============================================================================
# Load and Clean Phyloseq Object (Bacteria)
# =============================================================================
load("VU16s_roottraits_updated.Rdata")
VU16s_bacteria <- VU_16s_root_architecture_chem_dry
VU16s_bacteria <- prune_samples(sample_sums(VU16s_bacteria) > 0, VU16s_bacteria)

metadata_bacteria <- data.frame(sample_data(VU16s_bacteria)) %>% mutate(
  SRL = (Total.Root.Length.mm)/10 / g_dry,
  RTD = (Volume.mm3)/10 / g_dry,
  Avg_Diameter = Average.Diameter.mm,
  N__root_percent = as.numeric(N__root_percent),
  FineRoots = Root.Length.mm0.0.5/10 + Root.Length.mm0.5.2.0/10,
  CoarseRoots = Root.Length.mm2.0.max/10,
  Phosphorus = ifelse(Phosphorus == "High", 1, 0),
  Greenhouse_Block = as.factor(Greenhouse_Block),
  PlantShort = as.factor(PlantShort),
  Fine_to_Coarse = FineRoots / (CoarseRoots + 1)
) %>%
  na.omit() %>%
  filter(C_root_percent > 23)

# =============================================================================
# Subset for Wild-Domestic Pairs and Apply 95th Percentile Filter
# =============================================================================
metadata_bacteria_wd <- metadata_bacteria %>%
  filter(PlantShort %in% c("PS", "PW", "CR", "CA", "GS", "GM")) %>%
  mutate(Pair = case_when(
    PlantShort %in% c("PS", "PW") ~ "Pair1",
    PlantShort %in% c("CR", "CA") ~ "Pair2",
    PlantShort %in% c("GS", "GM") ~ "Pair3"
  ))

root_traits <- c("SRL", "RTD", "Avg_Diameter", "N__root_percent", "Fine_to_Coarse", "C_root_percent")

for (trait in root_traits) {
  perc_95 <- quantile(metadata_bacteria_wd[[trait]], 0.95, na.rm = TRUE)
  metadata_bacteria_wd <- metadata_bacteria_wd %>%
    filter(.data[[trait]] <= perc_95)
}

# =============================================================================
# Scale and Standardize Root Traits
# =============================================================================
scaled_metadata_bacteria_wd <- metadata_bacteria_wd %>%
  mutate(across(all_of(root_traits), scale))

# =============================================================================
# Barplot with Mean ± SD
# =============================================================================
scaled_root_traits_long <- scaled_metadata_bacteria_wd %>%
  select(PlantShort, Domestication_status, all_of(root_traits)) %>%
  pivot_longer(cols = all_of(root_traits), names_to = "Trait", values_to = "Scaled_Value") %>%
  filter(Scaled_Value > -3)

root_traits_barplot_bacteria <- ggplot(scaled_root_traits_long, aes(x = PlantShort, y = Scaled_Value, fill = Domestication_status)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.9), alpha = 1, color = "black") +
  geom_jitter(aes(color = "black"), 
              position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9), 
              size = 0.5, alpha = 0.6) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", 
               position = position_dodge(width = 0.9), width = 0.2) +
  facet_wrap(~ Trait, scales = "free_y", ncol = 3) +
  ggsci::scale_color_aaas() +
  ggsci::scale_fill_aaas() +
  scale_shape_manual(values = c("Wild" = 17, "Domestic" = 16)) +
  labs(
    title = "Filtered and Scaled Root Traits (Bacteria)",
    x = "Plant Species",
    y = "Scaled Trait Value",
    fill = "Plant Species",
    color = "Plant Species",
    shape = "Domestication Status"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#ggsave(
#  filename = file.path(save_dir, "Filtered_Scaled_Root_Traits_Barplots_Bacteria.pdf"),
  plot = root_traits_barplot_bacteria,
  width = 12,
  height = 8,
  dpi = 300, scale = 0.5
)

# =============================================================================
# Rename Pairs in Metadata
# =============================================================================
scaled_metadata_bacteria_wd <- scaled_metadata_bacteria_wd %>%
  mutate(Pair = recode(Pair,
                       "Pair1" = "Pisum",
                       "Pair2" = "Cicer",
                       "Pair3" = "Glycine"
  ))

# =============================================================================
# Loop Through Traits and Store Pairwise Contrasts
# =============================================================================
traits <- c("SRL", "RTD", "Avg_Diameter", "N__root_percent", "Fine_to_Coarse", "C_root_percent")
contrast_list <- list()

for (trait in traits) {
  model <- lmer(
    as.formula(paste(trait, "~ Domestication_status * Pair + Phosphorus + (1 | Greenhouse_Block)")),
    data = scaled_metadata_bacteria_wd
  )
  
  emm <- emmeans(model, ~ Domestication_status | Pair)
  contrast <- contrast(emm, method = "pairwise")
  
  contrast_df <- as.data.frame(contrast) %>%
    mutate(
      Trait = trait,
      Pair = recode(Pair,  # Ensure pairs in contrasts are renamed too
                    "Pair1" = "Pisum",
                    "Pair2" = "Cicer",
                    "Pair3" = "Glycine"
      )
    )
  
  contrast_list[[trait]] <- contrast_df
}

# =============================================================================
# Apply FDR Correction to P-values
# =============================================================================
all_contrasts <- bind_rows(contrast_list) %>%
  dplyr::mutate(
    Pair = dplyr::recode(Pair,
                         "Pair1" = "Pisum",
                         "Pair2" = "Cicer",
                         "Pair3" = "Glycine"
    ),
    p.value = signif(p.value, 3),
    p.adj = p.adjust(p.value, method = "fdr"),
    p.adj = signif(p.adj, 3),
    color = case_when(
      p.adj < 0.05 ~ Pair,
      TRUE ~ "NS"
    )
  )

# =============================================================================
# Set Color for Each Genus (Pair) and NS
# =============================================================================
genus_colors <- c("Pisum" = "#7A3DA5",  # Purple
                  "Cicer" = "#2A9D8F",  # Teal
                  "Glycine" = "#FFC75F",  # Yellow
                  "NS" = "lightgrey"  # Light grey for non-significant results
)

# =============================================================================
# Plot: Y-axis as Trait, X-axis as Estimate, Dodge by Genus (Pair)
# =============================================================================
ggplot(all_contrasts, aes(x =estimate, y = Trait, color = color)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +  # Points for estimates
  geom_errorbarh(aes(xmin = estimate - SE, xmax = estimate + SE), 
                 position = position_dodge(width = 0.6), height = 0.2) +  # Horizontal error bars
  scale_color_manual(values = genus_colors) +  # Set custom colors for each genus (Pair) and NS
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Dashed line at 0
  labs(
    title = "Pairwise Contrasts (Wild vs. Domestic) per Trait",
    x = "Estimate (Wild - Domestic)",
    y = "Trait",
    color = "Genus"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis text
        legend.position = "right")  # Display legend for color

# Save the plot as a file
ggsave(
  filename = file.path(save_dir, "Domestication_Pairwise_Contrasts_Estimates_Colored_By_Genus_NS_Gray.pdf"),
  plot = last_plot(),
  width = 8,
  height = 8,
  dpi = 300,
  scale = 0.5
) 
