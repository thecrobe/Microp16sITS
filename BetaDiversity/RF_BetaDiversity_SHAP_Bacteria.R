# =============================================================================
# Load Required Libraries
# =============================================================================
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
# Load and process bacteria data
# =============================================================================
load("VU16s_roottraits_updated.Rdata")
VU16s_bacteria <- VU_16s_root_architecture_chem_dry
VU16s_bacteria <- prune_samples(sample_sums(VU16s_bacteria) > 0, VU16s_bacteria)
metadata_bacteria_pcoa <- data.frame(sample_data(VU16s_bacteria))
metadata_bacteria_pcoa$N__root_percent <- as.numeric(metadata_bacteria_pcoa$N__root_percent)

metadata_bacteria_pcoa <- metadata_bacteria_pcoa %>%
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

physeq_obj <- prune_samples(rownames(metadata_bacteria_pcoa), VU16s_bacteria)
otu_matrix <- as.matrix(otu_table(physeq_obj))
if (taxa_are_rows(physeq_obj)) otu_matrix <- t(otu_matrix)
bray_dist <- vegdist(otu_matrix, method = "bray")
pcoa_res <- cmdscale(bray_dist, eig = TRUE, k = 5)
pcoa_scores <- as.data.frame(pcoa_res$points)
colnames(pcoa_scores) <- paste0("PCoA", 1:5)
metadata_bacteria_pcoa <- bind_cols(metadata_bacteria_pcoa, pcoa_scores)

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
feature_names <- c(predictors, "PlantShort")

data_bacteria <- metadata_bacteria_pcoa[, c("PCoA1", "PCoA2", predictors, "PlantShort")]
data_bacteria$PlantShort <- factor(data_bacteria$PlantShort, levels = levels(metadata_bacteria_pcoa$PlantShort))

data_bacteria_filtered <- scale_center_filter(data_bacteria, predictors) %>%
  filter(abs(C_root_percent) <= 2)

X_bacteria_filtered <- data_bacteria_filtered[, feature_names]

# =============================================================================
# Fit ranger models using same features
# =============================================================================
set.seed(123)
model_ranger_pcoa1 <- ranger(
  formula = as.formula(paste("PCoA1 ~", paste(feature_names, collapse = " + "))),
  data = data_bacteria_filtered,
  importance = "permutation",
  num.trees = 500
)

set.seed(123)
model_ranger_pcoa2 <- ranger(
  formula = as.formula(paste("PCoA2 ~", paste(feature_names, collapse = " + "))),
  data = data_bacteria_filtered,
  importance = "permutation",
  num.trees = 500
)

# -----------------------------------------------------------------------------
# Fit models to show OOB error, easier with this package.
# -----------------------------------------------------------------------------
rf_data <- data_bacteria_filtered[, c("PCoA1", "PCoA2", predictors, "PlantShort")]

rf_pcoa1 <- randomForest(
  PCoA1 ~ .,
  data = rf_data,
  ntree = 500,
  importance = TRUE
)

rf_pcoa2 <- randomForest(
  PCoA2 ~ .,
  data = rf_data,
  ntree = 500,
  importance = TRUE
)

# -----------------------------------------------------------------------------
# Extract OOB MSE per tree
# -----------------------------------------------------------------------------
err_df <- bind_rows(
  data.frame(
    Trees = seq_len(rf_pcoa1$ntree),
    OOB_Error = rf_pcoa1$mse,
    PC = "PCoA1"
  ),
  data.frame(
    Trees = seq_len(rf_pcoa2$ntree),
    OOB_Error = rf_pcoa2$mse,
    PC = "PCoA2"
  )
)

# -----------------------------------------------------------------------------
# Plot combined convergence figure
# -----------------------------------------------------------------------------
ggplot(err_df, aes(x = Trees, y = OOB_Error)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ PC, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Random Forest convergence - bacteria",
    x = "Number of trees",
    y = "Out of bag mean squared error"
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
  feature_names = feature_names,
  X = X_bacteria_filtered,
  pred_wrapper = pfun,
  nsim = 100
)

set.seed(123)
shap_pcoa2 <- fastshap::explain(
  object = model_ranger_pcoa2,
  feature_names = feature_names,
  X = X_bacteria_filtered,
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
  select(Sample, all_of(feature_names))

shap_with_values <- shap_combined %>%
  left_join(predictor_values, by = "Sample")

# =============================================================================
# Save individual SHAP ~ Scaled Trait Value LOESS plots
# =============================================================================
shap_save_dir <- "RF/SHAP/bacteria/"
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
      title = paste("SHAP Values vs Scaled", feat, "(bacteria)"),
      x = paste(feat, "(scaled value)"),
      y = "SHAP value",
      color = "PCoA Axis"
    )
  
  ggsave(
    filename = file.path(shap_save_dir, paste0("SHAP_Scaled_", feat, "_bacteria.pdf")),
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
    title = "SHAP Values vs Scaled Trait Values (bacteria)",
    x = "Scaled Trait Value",
    y = "SHAP value",
    color = "PCoA Axis"
  )

ggsave(
  filename = file.path(shap_save_dir, "SHAP_All_Traits_Scaled_bacteria.pdf"),
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
  filename = file.path(shap_save_dir, "SHAP_Mean_Absolute_Barplot_bacteria.pdf"),
  width = 8,
  height = 6,device = "pdf",
  dpi = 300,scale=0.5
)

