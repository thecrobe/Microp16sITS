# =============================================================================
# Libraries
# =============================================================================
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(lme4)
library(lmerTest)
library(emmeans)
library(tibble)
library(readr)
library(marginaleffects)

set.seed(123)

# =============================================================================
# Paths
# =============================================================================
save_dir <- ""
dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Load phyloseq object
# =============================================================================
load("VU16s_roottraits_updated.Rdata")

VU16s_bacteria <- VU_16s_root_architecture_chem_dry
VU16s_bacteria <- prune_samples(sample_sums(VU16s_bacteria) > 0, VU16s_bacteria)

# =============================================================================
# Build metadata
# =============================================================================
metadata_bacteria <- data.frame(sample_data(VU16s_bacteria)) %>%
  mutate(
    SRL = (Total.Root.Length.mm / 10) / g_dry,
    RTD = (Volume.mm3 / 10) / g_dry,
    Avg_Diameter = Average.Diameter.mm,
    N__root_percent = as.numeric(N__root_percent),
    FineRoots = Root.Length.mm0.0.5 / 10 + Root.Length.mm0.5.2.0 / 10,
    CoarseRoots = Root.Length.mm2.0.max / 10,
    Phosphorus = ifelse(Phosphorus == "High", 1, 0),
    Greenhouse_Block = factor(Greenhouse_Block),
    PlantShort = factor(PlantShort),
    Fine_to_Coarse = FineRoots / (CoarseRoots + 1)
  ) %>%
  na.omit() %>%
  filter(C_root_percent > 23)

# =============================================================================
# Subset wild-domestic pairs
# =============================================================================
metadata_bacteria_wd <- metadata_bacteria %>%
  filter(PlantShort %in% c("PS", "PW", "CR", "CA", "GS", "GM")) %>%
  mutate(
    Pair = case_when(
      PlantShort %in% c("PS", "PW") ~ "Pisum",
      PlantShort %in% c("CR", "CA") ~ "Cicer",
      PlantShort %in% c("GS", "GM") ~ "Glycine"
    ),
    Pair = factor(Pair, levels = c("Pisum", "Cicer", "Glycine"))
  )

# =============================================================================
# Traits
# =============================================================================
root_traits <- c(
  "SRL",
  "RTD",
  "Avg_Diameter",
  "N__root_percent",
  "Fine_to_Coarse",
  "C_root_percent"
)

# =============================================================================
# Apply 95th percentile filter trait-by-trait
# =============================================================================
for (trait in root_traits) {
  cutoff_95 <- quantile(metadata_bacteria_wd[[trait]], 0.95, na.rm = TRUE)
  metadata_bacteria_wd <- metadata_bacteria_wd %>%
    filter(.data[[trait]] <= cutoff_95)
}

# =============================================================================
# Scale root traits
# =============================================================================
scaled_metadata_bacteria_wd <- metadata_bacteria_wd %>%
  mutate(across(all_of(root_traits), ~ as.numeric(scale(.x))))

# =============================================================================
# Optional barplot
# =============================================================================
scaled_root_traits_long <- scaled_metadata_bacteria_wd %>%
  select(PlantShort, Domestication_status, all_of(root_traits)) %>%
  pivot_longer(
    cols = all_of(root_traits),
    names_to = "Trait",
    values_to = "Scaled_Value"
  ) %>%
  filter(Scaled_Value > -3)

root_traits_barplot_bacteria <- ggplot(
  scaled_root_traits_long,
  aes(x = PlantShort, y = Scaled_Value, fill = Domestication_status)
) +
  stat_summary(
    fun = mean,
    geom = "bar",
    position = position_dodge(width = 0.9),
    color = "black"
  ) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9),
    size = 0.5,
    alpha = 0.6
  ) +
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    position = position_dodge(width = 0.9),
    width = 0.2
  ) +
  facet_wrap(~ Trait, scales = "free_y", ncol = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Filtered and Scaled Root Traits",
    x = "Plant species",
    y = "Scaled trait value",
    fill = "Domestication"
  )

ggsave(
  filename = file.path(save_dir, "Filtered_Scaled_Root_Traits_Barplots_Bacteria.pdf"),
  plot = root_traits_barplot_bacteria,
  width = 12,
  height = 8,
  dpi = 300
)

# =============================================================================
# Fit models and extract outputs
# =============================================================================
model_list <- list()
anova_list <- list()
contrast_list <- list()
marginal_list <- list()

for (trait in root_traits) {
  
  formula_i <- as.formula(
    paste(trait, "~ Domestication_status * Pair + Phosphorus + (1 | Greenhouse_Block)")
  )
  
  model_i <- lmer(formula_i, data = scaled_metadata_bacteria_wd)
  model_list[[trait]] <- model_i
  
  # ---------------------------------------------------------------------------
  # Full ANOVA table
  # ---------------------------------------------------------------------------
  anova_i <- anova(model_i) %>%
    data.frame() %>%
    rownames_to_column("Effect") %>%
    mutate(Trait = trait) %>%
    relocate(Trait, .before = Effect)
  
  anova_list[[trait]] <- anova_i
  
  # ---------------------------------------------------------------------------
  # Wild vs domestic contrast within each genus
  # ---------------------------------------------------------------------------
  emm_i <- emmeans(model_i, ~ Domestication_status | Pair)
  
  contrast_i <- contrast(emm_i, method = "pairwise") %>%
    data.frame() %>%
    mutate(Trait = trait) %>%
    relocate(Trait, .before = Pair)
  
  contrast_list[[trait]] <- contrast_i
  
  # ---------------------------------------------------------------------------
  # Marginal domestication effect averaged across genera
  # ---------------------------------------------------------------------------
  avg_i <- avg_predictions(model_i, by = "Domestication_status") %>%
    data.frame() %>%
    mutate(Trait = trait) %>%
    relocate(Trait)
  
  marginal_list[[trait]] <- avg_i
}

# =============================================================================
# Combine outputs
# =============================================================================
anova_all <- bind_rows(anova_list) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 4))
  )

contrast_all <- bind_rows(contrast_list) %>%
  mutate(
    p.value = signif(p.value, 3),
    p.adj = signif(p.adjust(p.value, method = "fdr"), 3)
  )

marginal_all <- bind_rows(marginal_list) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 4))
  )

# =============================================================================
# Add FDR-adjusted p values to ANOVA table
# =============================================================================
anova_all <- anova_all %>%
  group_by(Effect) %>%
  mutate(
    p.adj = p.adjust("Pr..F.", method = "fdr"),
    p.adj = signif(p.adj, 3)
  ) %>%
  ungroup() %>%
  mutate(
    `Pr(>F)` = signif(`Pr(>F)`, 3)
  )

# =============================================================================
# Manuscript-ready ANOVA summary table
# =============================================================================
anova_summary <- anova_all %>%
  select(Trait, Effect, NumDF, DenDF, "F.value",Pr..F.) %>%
  mutate(
    across(c(NumDF, DenDF, F.value), ~ round(.x, 2))
  ) %>%
  arrange(Trait, Effect)

# =============================================================================
# Manuscript-ready marginal domestication table
# =============================================================================
marginal_summary <- marginal_all %>%
  select(Trait, Domestication_status, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 3))
  ) %>%
  arrange(Trait, Domestication_status)

# =============================================================================
# Save tables as CSV
# =============================================================================
write_csv(anova_summary, file.path(save_dir, "Trait_models_ANOVA_summary.csv"))
write_csv(contrast_all, file.path(save_dir, "Trait_models_pairwise_contrasts.csv"))
write_csv(marginal_summary, file.path(save_dir, "Trait_models_marginal_domestication_effects.csv"))

# =============================================================================
# Also print clean tables to console
# =============================================================================
anova_summary
contrast_all
marginal_summary