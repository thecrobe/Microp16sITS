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
library(patchwork)
library(progress)
library(progressr)
library(car)

set.seed(123)

# =============================================================================
# Define save directory
# =============================================================================
save_dir <- "RF"
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

# =============================================================================
# FUNCTION: Process Data for RF 
# =============================================================================
process_data <- function(physeq_obj) {
  metadata <- data.frame(sample_data(physeq_obj))
  metadata$N__root_percent <- as.numeric(metadata$N__root_percent)
  
  metadata <- metadata %>%
    mutate(
      SRL = (Total.Root.Length.mm)/10 / g_dry,
      RTD = (Volume.mm3)/10 / g_dry,
      Avg_Diameter = Average.Diameter.mm,
      N__root_percent = N__root_percent,
      FineRoots = Root.Length.mm0.0.5/10 + Root.Length.mm0.5.2.0/10,
      CoarseRoots = Root.Length.mm2.0.max/10,
      Phosphorus = ifelse(Phosphorus == "High", 1, 0),  # Convert to binary 0/1
      Greenhouse_Block = as.factor(Greenhouse_Block),
      PlantShort = as.factor(PlantShort),
      Fine_to_Coarse = FineRoots / (CoarseRoots + 1)
    ) %>% na.omit()
  
  physeq_obj <- prune_samples(rownames(metadata), physeq_obj)
  
  otu_matrix <- as.matrix(otu_table(physeq_obj))
  if (taxa_are_rows(physeq_obj)) otu_matrix <- t(otu_matrix)
  bray_dist <- vegdist(otu_matrix, method = "bray")
  pcoa_res <- cmdscale(bray_dist, eig = TRUE, k = 5)
  pcoa_scores <- as.data.frame(pcoa_res$points)
  colnames(pcoa_scores) <- paste0("PCoA", 1:5)
  
  return(bind_cols(metadata, pcoa_scores))
}

# =============================================================================
# LOAD & PROCESS DATA
# =============================================================================
load("/VU16s_roottraits_updated.Rdata")
metadata_bact_pcoa <- process_data(VU_16s_root_architecture_chem_dry) 
metadata_bact_pcoa<-metadata_bact_pcoa %>% filter(C_root_percent >0)
load("/VU_ITS_roottraits_updated.Rdata")
metadata_fungi_pcoa <- process_data(VU_ITS_filtered_joined) 

# =============================================================================
# Apply Outlier Removal
# =============================================================================
remove_outliers <- function(data, predictors) {
  data_filtered <- data
  
  for (var in predictors) {
    threshold <- quantile(data[[var]], 0.95, na.rm = TRUE)  # Compute 95th percentile
    data_filtered <- data_filtered %>% filter(.data[[var]] <= threshold)
  }
  
  return(data_filtered)
}
pred_incl <- c("SRL", "RTD", "Avg_Diameter","Fine_to_Coarse", "C_root_percent", "N__root_percent")

metadata_bact_pcoa <- remove_outliers(metadata_bact_pcoa, pred_incl)
metadata_fungi_pcoa <- remove_outliers(metadata_fungi_pcoa, pred_incl)

# =============================================================================
# Compute VIF for Each Root Trait
# =============================================================================
calculate_vif_each <- function(data, predictors) {
  results <- list()
  for (response_var in predictors) {
    pred_vars <- setdiff(predictors, response_var)
    if (length(pred_vars) > 1) {
      formula <- as.formula(paste(response_var, "~", paste(pred_vars, collapse = " + ")))
      vif_model <- lm(formula, data = data)
      vif_values <- vif(vif_model)
      results[[response_var]] <- data.frame(Response = response_var, Predictor = names(vif_values), VIF = vif_values)
    }
  }
  return(bind_rows(results))
}

vif_results_each <- calculate_vif_each(
  metadata_bact_pcoa, 
  c("SRL", "RTD", "Avg_Diameter", "Fine_to_Coarse", "C_root_percent", "N__root_percent")
)

write.csv(vif_results_each, file = file.path(save_dir, "VIF_Results.csv"), row.names = FALSE)

# =============================================================================
# FUNCTION: Scale & Center Predictors
# =============================================================================
scale_and_filter <- function(df, predictors) {
  df_scaled <- df  # Copy input dataframe
  
  # Scale and center numeric variables
  df_scaled[predictors] <- scale(df[predictors], center = TRUE, scale = TRUE)
  
  # Compute 95th percentile threshold for each predictor
  percentile_95 <- apply(df_scaled[predictors], 2, function(x) quantile(x, 0.95, na.rm = TRUE))
  
  # Remove rows where any predictor exceeds its 95th percentile
  df_filtered <- df_scaled %>%
    filter(rowSums(sapply(predictors, function(var) df_scaled[[var]] > percentile_95[var])) == 0)
  
  return(df_filtered)
}

# Apply function to your dataset
predictors <- c("SRL", "RTD", "Avg_Diameter", "Fine_to_Coarse", "C_root_percent", "N__root_percent")
metadata_bact_pcoa <- scale_and_filter(metadata_bact_pcoa, predictors)
metadata_fungi_pcoa <- scale_and_filter(metadata_fungi_pcoa, predictors)
summary(metadata_bact_pcoa$C_root_percent)

# =============================================================================
# Run models
# =============================================================================
train_models <- function(metadata, root_traits, k_folds = 10) {
  
  # Model 1: Root Traits Only (No PlantShort)
  formula_roots <- as.formula(paste("PCoA1 ~", paste(root_traits, collapse = " + ")))
  formula_roots_pcoa2 <- as.formula(paste("PCoA2 ~", paste(root_traits, collapse = " + ")))
  
  # Model 2: PlantShort Only
  formula_plant <- as.formula("PCoA1 ~ PlantShort")
  formula_plant_pcoa2 <- as.formula("PCoA2 ~ PlantShort")
  
  # Model 3: Root Traits + PlantShort (All Predictors)
  formula_all <- as.formula(paste("PCoA1 ~", paste(c(root_traits, "PlantShort"), collapse = " + ")))
  formula_all_pcoa2 <- as.formula(paste("PCoA2 ~", paste(c(root_traits, "PlantShort"), collapse = " + ")))
  
  # Generate cross-validation folds by Greenhouse_Block
  set.seed(123)  # Ensures reproducibility
  folds <- groupKFold(metadata$Greenhouse_Block, k = k_folds)  
  
  train_ctrl <- trainControl(
    method = "cv",
    index = folds,  # Use grouped folds based on Greenhouse_Block
    savePredictions = "final",
    allowParallel = TRUE
  )
  
  models <- list(
    Roots_Only_PCoA1 = train(formula_roots, data = metadata, method = "rf", trControl = train_ctrl, importance = TRUE),
    Roots_Only_PCoA2 = train(formula_roots_pcoa2, data = metadata, method = "rf", trControl = train_ctrl, importance = TRUE),
    
    PlantShort_Only_PCoA1 = train(formula_plant, data = metadata, method = "rf", trControl = train_ctrl, importance = TRUE),
    PlantShort_Only_PCoA2 = train(formula_plant_pcoa2, data = metadata, method = "rf", trControl = train_ctrl, importance = TRUE),
    
    All_Predictors_PCoA1 = train(formula_all, data = metadata, method = "rf", trControl = train_ctrl, importance = TRUE),
    All_Predictors_PCoA2 = train(formula_all_pcoa2, data = metadata, method = "rf", trControl = train_ctrl, importance = TRUE)
  )
  
  return(models)
}

# Define root traits (excluding PlantShort)
root_traits <- c("SRL", "RTD", "Avg_Diameter", "Fine_to_Coarse", "C_root_percent", "N__root_percent", "Phosphorus")

# Train models for Bacteria and Fungi datasets
models_bact <- train_models(metadata_bact_pcoa, root_traits)
models_fungi <- train_models(metadata_fungi_pcoa, root_traits)

extract_performance <- function(models, dataset) {
  results <- bind_rows(lapply(names(models), function(pc) {
    res <- models[[pc]]$results
    best_idx <- which.max(res$Rsquared)
    data.frame(Dataset = dataset, Model_Type = pc, RMSE = res$RMSE[best_idx], Rsquared = res$Rsquared[best_idx])
  }))
  return(results)
}

performance_results <- bind_rows( 
  extract_performance(models_bact, "Bacteria"),
  extract_performance(models_fungi, "Fungi")
)

write.csv(performance_results, file = file.path(save_dir, "RF_Model_Performance.csv"), row.names = FALSE)
