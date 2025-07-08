library(phyloseq)
library(vegan)
library(dplyr)
library(tidyr)
library(caret)
library(ggplot2)
library(forcats)
library(ranger)
library(fastshap)

set.seed(123)

# =============================================================================
# Load and process data
# =============================================================================
load("/VU16s_roottraits_updated.Rdata")
metadata_bact_pcoa <- data.frame(sample_data(VU_16s_root_architecture_chem_dry))
metadata_bact_pcoa$N__root_percent <- as.numeric(metadata_bact_pcoa$N__root_percent)

metadata_bact_pcoa <- metadata_bact_pcoa %>%
  mutate(
    SRL = (Total.Root.Length.mm)/10 / g_dry,
    RTD = (Volume.mm3)/10 / g_dry,
    Avg_Diameter = Average.Diameter.mm,
    N__root_percent = N__root_percent,
    FineRoots = Root.Length.mm0.0.5/10 + Root.Length.mm0.5.2.0/10,
    CoarseRoots = Root.Length.mm2.0.max/10,
    Phosphorus = ifelse(Phosphorus == "High", 1, 0),
    Greenhouse_Block = as.factor(Greenhouse_Block),
    PlantShort = as.factor(PlantShort),
    Fine_to_Coarse = FineRoots / (CoarseRoots + 1)
  ) %>%
  na.omit() %>%
  filter(C_root_percent > 0) 

physeq_obj <- prune_samples(rownames(metadata_bact_pcoa), VU_16s_root_architecture_chem_dry)
otu_matrix <- as.matrix(otu_table(physeq_obj))
if (taxa_are_rows(physeq_obj)) otu_matrix <- t(otu_matrix)
bray_dist <- vegdist(otu_matrix, method = "bray")
pcoa_res <- cmdscale(bray_dist, eig = TRUE, k = 5)
pcoa_scores <- as.data.frame(pcoa_res$points)
colnames(pcoa_scores) <- paste0("PCoA", 1:5)
metadata_bact_pcoa <- bind_cols(metadata_bact_pcoa, pcoa_scores)

# =============================================================================
# Scale, center, and apply 95th percentile filter
# =============================================================================
scale_center_filter <- function(df, predictors) {
  df_scaled <- df %>%
    mutate(across(all_of(predictors), ~ scale(.x)[,1]))
  
  percentiles <- df_scaled %>%
    summarise(across(all_of(predictors), ~ quantile(.x, 0.95, na.rm = TRUE)))
  
  for (var in predictors) {
    df_scaled <- df_scaled %>% filter(.data[[var]] <= percentiles[[var]])
  }
  
  return(df_scaled)
}

predictors <- c("SRL", "RTD", "Avg_Diameter", "Fine_to_Coarse", 
                "C_root_percent", "N__root_percent", "Phosphorus")

responses <- c("PCoA1", "PCoA2")

data_bacteria <- metadata_bact_pcoa[, c(responses, predictors, "PlantShort")]
data_bacteria$PlantShort <- factor(
  data_bacteria$PlantShort, 
  levels = levels(metadata_bact_pcoa$PlantShort)
)

data_bacteria_filtered <- scale_center_filter(data_bacteria, predictors)
data_bacteria_filtered <- data_bacteria_filtered %>%
  filter(abs(C_root_percent) <= 2)

X_bacteria_filtered <- data_bacteria_filtered[, predictors]

# =============================================================================
# Fit ranger models
# =============================================================================
set.seed(123)
model_ranger_pcoa1 <- ranger(
  formula = PCoA1 ~ ., 
  data = data_bacteria_filtered[, c("PCoA1", predictors, "PlantShort")], 
  importance = "permutation",
  num.trees = 500
)

model_ranger_pcoa2 <- ranger(
  formula = PCoA2 ~ ., 
  data = data_bacteria_filtered[, c("PCoA2", predictors, "PlantShort")], 
  importance = "permutation",
  num.trees = 500
)

# =============================================================================
# Define prediction function
# =============================================================================
pfun <- function(object, newdata) {
  predict(object, data = newdata)$predictions
}

# =============================================================================
# Compute SHAP values
# =============================================================================
set.seed(123)
shap_pcoa1 <- fastshap::explain(
  object = model_ranger_pcoa1,
  feature_names = c(predictors, "PlantShort"),
  X = data_bacteria_filtered[, c(predictors, "PlantShort")],
  pred_wrapper = pfun,
  nsim = 100
)

set.seed(123)
shap_pcoa2 <- fastshap::explain(
  object = model_ranger_pcoa2,
  feature_names = c(predictors, "PlantShort"),
  X = data_bacteria_filtered[, c(predictors, "PlantShort")],
  pred_wrapper = pfun,
  nsim = 100
)

# =============================================================================
# Reshape and combine SHAP values
# =============================================================================
shap_long_pcoa1 <- data.frame(shap_pcoa1) %>%
  mutate(Sample = row_number(), PC = "PCoA1") %>%
  pivot_longer(cols = -c(Sample, PC), names_to = "Feature", values_to = "SHAP")

shap_long_pcoa2 <- data.frame(shap_pcoa2) %>%
  mutate(Sample = row_number(), PC = "PCoA2") %>%
  pivot_longer(cols = -c(Sample, PC), names_to = "Feature", values_to = "SHAP")

shap_combined <- bind_rows(shap_long_pcoa1, shap_long_pcoa2)

# =============================================================================
# Merge SHAP values with scaled predictor values
# =============================================================================
predictor_values <- data_bacteria_filtered %>%
  mutate(Sample = row_number()) %>%
  select(Sample, all_of(predictors), PlantShort)

shap_with_values <- shap_combined %>%
  left_join(predictor_values, by = "Sample")

# =============================================================================
# Save individual SHAP ~ Scaled Trait Value LOESS plots
# =============================================================================
shap_save_dir <- "/Users/justinstewart/Dropbox/Collaborations/FamilyExperiment/Data/Code/New/Feb2024/Figures/BetaDiversity/RF/SHAP/Bacteria/"
if (!dir.exists(shap_save_dir)) dir.create(shap_save_dir, recursive = TRUE)

features <- unique(shap_with_values$Feature)

for (feat in features) {
  p <- ggplot(shap_with_values %>% filter(Feature == feat), 
              aes(x = .data[[feat]], y = SHAP, color = PC)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = "loess", se = TRUE) +
    scale_color_manual(values = c("PCoA1" = "red", "PCoA2" = "cornflowerblue")) +
    theme_minimal() +
    labs(
      title = paste("SHAP Values vs Scaled", feat),
      x = paste(feat, "(scaled value)"),
      y = "SHAP value",
      color = "PCoA Axis"
    )
  
  ggsave(
    filename = file.path(shap_save_dir, paste0("SHAP_Scaled_", feat, ".pdf")),
    plot = p,device = "pdf",
    width = 6,
    height = 4,
    dpi = 300
  )
}

# =============================================================================
# Multi-panel SHAP plot
# =============================================================================
multi_panel_plot <- ggplot(shap_with_values, aes(x = .data[[Feature]], y = SHAP, color = PC)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "loess", se = TRUE) +
  facet_wrap(~ Feature, scales = "free_x") +
  scale_color_manual(values = c("PCoA1" = "red", "PCoA2" = "cornflowerblue")) +
  theme_minimal() +
  labs(
    title = "SHAP Values vs Scaled Trait Values (All Features)",
    x = "Scaled Trait Value",
    y = "SHAP value",
    color = "PCoA Axis"
  )

ggsave(
  filename = file.path(shap_save_dir, "SHAP_All_Traits_Scaled.pdf"),
  plot = multi_panel_plot,device = "pdf",
  width = 12,
  height = 10,
  dpi = 300
)

# =============================================================================
# Barplot of mean absolute SHAP values per PCoA axis
# =============================================================================

shap_bar_data <- shap_combined %>%
  mutate(Abs_SHAP = abs(SHAP)) %>%
  group_by(PC, Feature) %>%
  summarise(Mean_Abs_SHAP = mean(Abs_SHAP, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    PC = factor(PC, levels = c("PCoA1", "PCoA2")),
    Feature = fct_reorder(Feature, Mean_Abs_SHAP, .fun = max)
  )

ggplot(shap_bar_data, aes(y = Feature, x = Mean_Abs_SHAP, fill = PC)) +
  geom_col(position = "dodge", alpha = 0.9) +
  scale_fill_manual(values = c("PCoA1" = "red", "PCoA2" = "cornflowerblue")) +
  theme_minimal() +
  labs(
    title = "Mean Absolute SHAP Values per Feature",
    x = "Mean Absolute SHAP Value",
    y = "Feature",
    fill = "PCoA Axis"
  ) +
  theme(legend.position = "bottom")

# Save plot
ggsave(
  filename = file.path(shap_save_dir, "SHAP_Mean_Absolute_Barplot.pdf"),
  width = 8,
  height = 6,device = "pdf",
  dpi = 300,scale=0.5
)
