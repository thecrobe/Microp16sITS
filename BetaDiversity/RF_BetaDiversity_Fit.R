# =============================================================================
# Load Required Libraries
# =============================================================================
library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(caret)
library(ape)
library(phytools)
library(patchwork)

set.seed(123)

# =============================================================================
# Load Data
# =============================================================================
load("VU16s_roottraits_updated.Rdata")
load("VU_ITS_roottraits_updated.Rdata")
tree <- read.nexus("Phylogeny/host_phylogeny_microp_legume_guest.nex")
save_dir <- "RF"

# =============================================================================
# Prepare Metadata Function
# =============================================================================
prepare_metadata <- function(physeq_obj) {
  metadata <- data.frame(sample_data(physeq_obj)) %>%
    mutate(
      SRL = (Total.Root.Length.mm)/10 / g_dry,
      RTD = (Volume.mm3)/10 / g_dry,
      Avg_Diameter = Average.Diameter.mm,
      FineRoots = Root.Length.mm0.0.5/10 + Root.Length.mm0.5.2.0/10,
      CoarseRoots = Root.Length.mm2.0.max/10,
      Fine_to_Coarse = FineRoots / (CoarseRoots + 1),
      Phosphorus = ifelse(Phosphorus == "High", 1, 0),
      Plant_Long = as.factor(Plant_Long),
      N__root_percent = as.numeric(N__root_percent)
    ) %>% na.omit()
  return(metadata)
}

scale_and_filter <- function(df, predictors) {
  df_scaled <- df
  df_scaled[predictors] <- scale(df[predictors], center = TRUE, scale = TRUE)
  percentile_95 <- apply(df_scaled[predictors], 2, function(x) quantile(x, 0.95, na.rm = TRUE))
  df_filtered <- df_scaled %>%
    filter(rowSums(sapply(predictors, function(var) df_scaled[[var]] > percentile_95[var])) == 0)
  return(df_filtered)
}

get_pcoa <- function(physeq_obj, metadata) {
  otu_mat <- as.matrix(otu_table(physeq_obj))
  if (taxa_are_rows(physeq_obj)) otu_mat <- t(otu_mat)
  bray <- vegdist(otu_mat, method = "bray")
  pcoa <- cmdscale(bray, eig = TRUE, k = 2)
  metadata <- metadata[match(rownames(pcoa$points), rownames(metadata)), ]
  metadata$PCoA1 <- pcoa$points[, 1]
  metadata$PCoA2 <- pcoa$points[, 2]
  return(metadata)
}

root_traits <- c("SRL", "RTD", "Avg_Diameter", "Fine_to_Coarse", "C_root_percent", "N__root_percent", "Phosphorus")

train_rf_models <- function(metadata) {
  folds <- groupKFold(metadata$Greenhouse_Block, k = 10)
  ctrl <- trainControl(method = "cv", index = folds, savePredictions = "final", allowParallel = TRUE)
  models <- list(
    Roots_Only_PCoA1 = train(PCoA1 ~ ., data = metadata[, c("PCoA1", root_traits)], method = "rf", trControl = ctrl),
    Roots_Only_PCoA2 = train(PCoA2 ~ ., data = metadata[, c("PCoA2", root_traits)], method = "rf", trControl = ctrl),
    All_Predictors_PCoA1 = train(PCoA1 ~ ., data = metadata[, c("PCoA1", root_traits, "Plant_Long")], method = "rf", trControl = ctrl),
    All_Predictors_PCoA2 = train(PCoA2 ~ ., data = metadata[, c("PCoA2", root_traits, "Plant_Long")], method = "rf", trControl = ctrl)
  )
  return(models)
}

metadata_bact <- prepare_metadata(VU_16s_root_architecture_chem_dry)
VU_16s_root_architecture_chem_dry <- prune_samples(rownames(metadata_bact), VU_16s_root_architecture_chem_dry)
metadata_bact <- get_pcoa(VU_16s_root_architecture_chem_dry, metadata_bact)
metadata_bact <- scale_and_filter(metadata_bact, root_traits)
models_bact <- train_rf_models(metadata_bact)

metadata_fungi <- prepare_metadata(VU_ITS_filtered_joined)
VU_ITS_filtered_joined <- prune_samples(rownames(metadata_fungi), VU_ITS_filtered_joined)
metadata_fungi <- get_pcoa(VU_ITS_filtered_joined, metadata_fungi)
metadata_fungi <- scale_and_filter(metadata_fungi, root_traits)
models_fungi <- train_rf_models(metadata_fungi)



compute_centroids <- function(metadata, model1, model2, dataset_name) {
  metadata <- metadata %>%
    mutate(
      PCoA1_pred = predict(model1, metadata),
      PCoA2_pred = predict(model2, metadata)
    )
  
  metadata %>%
    group_by(PlantShort) %>%
    summarise(
      PCoA1_obs = mean(PCoA1, na.rm = TRUE),
      PCoA2_obs = mean(PCoA2, na.rm = TRUE),
      PCoA1_pred = mean(PCoA1_pred, na.rm = TRUE),
      PCoA2_pred = mean(PCoA2_pred, na.rm = TRUE),
      Model = dataset_name
    )
}

centroids_all <- bind_rows(
  compute_centroids(metadata_bact, models_bact$Roots_Only_PCoA1, models_bact$Roots_Only_PCoA2, "Bacteria - Roots Only"),
  compute_centroids(metadata_bact, models_bact$All_Predictors_PCoA1, models_bact$All_Predictors_PCoA2, "Bacteria - Roots + Plant"),
  compute_centroids(metadata_fungi, models_fungi$Roots_Only_PCoA1, models_fungi$Roots_Only_PCoA2, "Fungi - Roots Only"),
  compute_centroids(metadata_fungi, models_fungi$All_Predictors_PCoA1, models_fungi$All_Predictors_PCoA2, "Fungi - Roots + Plant")
)


# =============================================================================
# Predicted vs Observed Plot Function with R²
# =============================================================================
plot_pred_vs_obs <- function(model, metadata, response, dataset) {
  preds <- predict(model, metadata)
  obs <- metadata[[response]]
  r2 <- round(summary(lm(obs ~ preds))$r.squared, 2)
  df <- data.frame(Predicted = preds, Observed = obs)
  ggplot(df, aes(x = Predicted, y = Observed)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", lwd = 1.5, color = "red") +
    annotate("text", x = min(df$Predicted), y = max(df$Observed),
             label = paste("R² =", r2), hjust = 0, vjust = 1, size = 4) +
    theme_minimal() +
    labs(
      title = dataset,
      x = "Predicted",
      y = "Observed"
    )
}

# =============================================================================
# Create Predicted vs Observed Plots for All Models
# =============================================================================
p1 <- plot_pred_vs_obs(models_bact$Roots_Only_PCoA1, metadata_bact, "PCoA1", "Bacteria - Roots Only PCoA1")
p2 <- plot_pred_vs_obs(models_bact$Roots_Only_PCoA2, metadata_bact, "PCoA2", "Bacteria - Roots Only PCoA2")
p3 <- plot_pred_vs_obs(models_bact$All_Predictors_PCoA1, metadata_bact, "PCoA1", "Bacteria - Roots + Plant PCoA1")
p4 <- plot_pred_vs_obs(models_bact$All_Predictors_PCoA2, metadata_bact, "PCoA2", "Bacteria - Roots + Plant PCoA2")
p5 <- plot_pred_vs_obs(models_fungi$Roots_Only_PCoA1, metadata_fungi, "PCoA1", "Fungi - Roots Only PCoA1")
p6 <- plot_pred_vs_obs(models_fungi$Roots_Only_PCoA2, metadata_fungi, "PCoA2", "Fungi - Roots Only PCoA2")
p7 <- plot_pred_vs_obs(models_fungi$All_Predictors_PCoA1, metadata_fungi, "PCoA1", "Fungi - Roots + Plant PCoA1")
p8 <- plot_pred_vs_obs(models_fungi$All_Predictors_PCoA2, metadata_fungi, "PCoA2", "Fungi - Roots + Plant PCoA2")

# Combine and Save Multipanel Plot
multi_panel_pred_obs <- (p1 + p2 + p3 + p4) / (p5 + p6 + p7 + p8)

ggsave(
  filename = file.path(save_dir, "RF_Pred_vs_Obs_All_Models.pdf"),
  plot = multi_panel_pred_obs,
  width = 14,
  height = 10,
  dpi = 300
)

# =============================================================================
# Compute Phylogenetic Signal on Model Residuals (Grouped by Plant_Long)
# =============================================================================
compute_phylo_signal <- function(model, metadata, response, tree) {
  preds <- predict(model, metadata)
  residuals <- metadata[[response]] - preds
  residuals_df <- data.frame(Plant_Long = metadata$Plant_Long, Residual = residuals) %>%
    group_by(Plant_Long) %>%
    summarise(Mean_Residual = mean(Residual, na.rm = TRUE)) %>%
    filter(Plant_Long %in% tree$tip.label)
  
  tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, residuals_df$Plant_Long))
  residual_vector <- residuals_df$Mean_Residual
  names(residual_vector) <- residuals_df$Plant_Long
  residual_vector <- residual_vector[tree_pruned$tip.label]
  
  K <- phylosig(tree_pruned, residual_vector, method = "K", test = TRUE)
  return(data.frame(K = round(K$K, 2), P = signif(K$P, 3)))
}

phylo_signal_results <- bind_rows(
  cbind(Model = "Bacteria - Roots Only", PCoA = "PCoA1", compute_phylo_signal(models_bact$Roots_Only_PCoA1, metadata_bact, "PCoA1", tree)),
  cbind(Model = "Bacteria - Roots Only", PCoA = "PCoA2", compute_phylo_signal(models_bact$Roots_Only_PCoA2, metadata_bact, "PCoA2", tree)),
  cbind(Model = "Bacteria - Roots + Plant", PCoA = "PCoA1", compute_phylo_signal(models_bact$All_Predictors_PCoA1, metadata_bact, "PCoA1", tree)),
  cbind(Model = "Bacteria - Roots + Plant", PCoA = "PCoA2", compute_phylo_signal(models_bact$All_Predictors_PCoA2, metadata_bact, "PCoA2", tree)),
  cbind(Model = "Fungi - Roots Only", PCoA = "PCoA1", compute_phylo_signal(models_fungi$Roots_Only_PCoA1, metadata_fungi, "PCoA1", tree)),
  cbind(Model = "Fungi - Roots Only", PCoA = "PCoA2", compute_phylo_signal(models_fungi$Roots_Only_PCoA2, metadata_fungi, "PCoA2", tree)),
  cbind(Model = "Fungi - Roots + Plant", PCoA = "PCoA1", compute_phylo_signal(models_fungi$All_Predictors_PCoA1, metadata_fungi, "PCoA1", tree)),
  cbind(Model = "Fungi - Roots + Plant", PCoA = "PCoA2", compute_phylo_signal(models_fungi$All_Predictors_PCoA2, metadata_fungi, "PCoA2", tree))
)

phylo_labels <- phylo_signal_results %>%
  group_by(Model) %>%
  summarise(
    Label = paste0(
      "PCoA1: K=", K[PCoA == "PCoA1"], ", P=", P[PCoA == "PCoA1"], "\n",
      "PCoA2: K=", K[PCoA == "PCoA2"], ", P=", P[PCoA == "PCoA2"]
    )
  )

centroids_all <- left_join(centroids_all, phylo_labels, by = "Model")

# =============================================================================
# Updated Centroid Plot with Phylogenetic Signal Labels
# =============================================================================

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
plot_centroids <- ggplot() +
  geom_point(data = centroids_all, aes(x = PCoA1_obs, y = PCoA2_obs, fill = PlantShort), size = 3, shape = 21) +
  geom_point(data = centroids_all, aes(x = PCoA1_pred, y = PCoA2_pred, fill = PlantShort), size = 3, shape = 24) +
  facet_wrap(~ Model, scales = "free") +
  geom_text(
    data = centroids_all %>% group_by(Model) %>% slice(1),
    aes(x = -Inf, y = Inf, label = Label),
    hjust = -0.1, vjust = 1.1, inherit.aes = FALSE, size = 4
  ) +
  labs(
    title = "Observed (circles) vs Predicted (triangles) Microbiome Centroids\nwith Phylogenetic Signal (K, P)",
    x = "PCoA1",
    y = "PCoA2"
  ) +
  theme_classic() +
  scale_fill_manual(values = custom_colors)

ggsave(
  filename = file.path(save_dir, "RF_PCoA_Centroids_PhyloSignal.pdf"),
  plot = plot_centroids,
  width = 14,
  height = 10,
  dpi = 300,scale = 0.5
)
