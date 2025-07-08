# ---------------------------------------------------
# Load Required Libraries
# ---------------------------------------------------
library(brms)
library(dplyr)
library(tidyverse)

# ---------------------------------------------------
# Define Model Paths
# ---------------------------------------------------
model_components <- c("PC1", "PC2", "Chao1")  # Richness included
model_results <- list()

# Base directory where models are stored
model_dir <- "/Users/justinstewart/Dropbox/Collaborations/FamilyExperiment/Data/Code/New/Feb2024/Host_phylogenetic_signal/Models/Bacteria/"

# ---------------------------------------------------
# Load All Models Per Component
# ---------------------------------------------------
for (component in model_components) {
  model_files <- list.files(model_dir, pattern = paste0(component, "_Bacteria_.*variance.RData"), full.names = TRUE)
  
  for (file in model_files) {
    load(file)
    model_results[[component]] <- model_list
  }
}

# ---------------------------------------------------
# Compute WAIC and Bayesian R² for Each Model
# ---------------------------------------------------
waic_results <- data.frame()

for (component in names(model_results)) {
  models <- model_results[[component]]
  
  for (model_name in names(models)) {
    model <- models[[model_name]]
    waic_info <- waic(model)$estimates
    r2_info <- bayes_R2(model)  # Compute Bayesian R²
    
    waic_results <- rbind(waic_results, data.frame(
      Component = component,
      Model = model_name,
      WAIC = waic_info["waic", "Estimate"],
      SE = waic_info["waic", "SE"],
      Bayes_R2_Mean = mean(r2_info[,1]),  # Extract mean Bayesian R²
      Bayes_R2_SD = (r2_info[,2])  # Extract SD of Bayesian R²
    ))
  }
}

# Order WAIC results per component (lowest to highest WAIC)
waic_results <- waic_results %>%
  group_by(Component) %>%
  arrange(Component, WAIC) %>%
  ungroup()

# Save WAIC results with Bayesian R²
write.csv(waic_results, file = paste0(model_dir, "Bacteria_WAIC_R2_results.csv"), row.names = FALSE)

print("WAIC and Bayesian R² values calculated, ordered, and saved.")
print(waic_results)

# ---------------------------------------------------
# Compute ICC for Species and Phylogeny for Each Component
# ---------------------------------------------------
icc_draws <- data.frame()

for (component in names(model_results)) {
  models <- model_results[[component]]
  
  relevant_models <- names(models)[grepl("Species", names(models)) | grepl("Phylogeny", names(models))]
  
  for (model_name in relevant_models) {
    model <- models[[model_name]]
    
    has_species <- grepl("Species", model_name)
    has_phylogeny <- grepl("Phylogeny", model_name)
    
    if (has_species & has_phylogeny) {
      hyp_phylo <- "sd_Phylogeny__Intercept^2 / (sd_Phylogeny__Intercept^2 + sd_Species__Intercept^2 + sigma^2) = 0"
      hyp_species <- "sd_Species__Intercept^2 / (sd_Phylogeny__Intercept^2 + sd_Species__Intercept^2 + sigma^2) = 0"
      
      hyps_phylo <- hypothesis(model, hyp_phylo, class = NULL)$samples
      hyps_species <- hypothesis(model, hyp_species, class = NULL)$samples
      
      df_phylo <- data.frame(Component = component, Model = model_name, Effect = "Phylogeny", Draws = hyps_phylo$H1)
      df_species <- data.frame(Component = component, Model = model_name, Effect = "Species", Draws = hyps_species$H1)
      
      icc_draws <- rbind(icc_draws, df_phylo, df_species)
      
    } else if (has_species) {
      hyp_species <- "sd_Species__Intercept^2 / (sd_Species__Intercept^2 + sigma^2) = 0"
      hyps_species <- hypothesis(model, hyp_species, class = NULL)$samples
      
      df_species <- data.frame(Component = component, Model = model_name, Effect = "Species", Draws = hyps_species$H1)
      icc_draws <- rbind(icc_draws, df_species)
      
    } else if (has_phylogeny) {
      hyp_phylo <- "sd_Phylogeny__Intercept^2 / (sd_Phylogeny__Intercept^2 + sigma^2) = 0"
      hyps_phylo <- hypothesis(model, hyp_phylo, class = NULL)$samples
      
      df_phylo <- data.frame(Component = component, Model = model_name, Effect = "Phylogeny", Draws = hyps_phylo$H1)
      icc_draws <- rbind(icc_draws, df_phylo)
    }
  }
}

# Save ICC Draws
write.csv(icc_draws, file = paste0(model_dir, "Bacteria_ICC_Draws.csv"), row.names = FALSE)
icc_draws<-read.csv(file = paste0(model_dir, "Bacteria_ICC_Draws.csv"))

# ---------------------------------------------------
# Compute Overall Mean and SD Across All Models Per Component & Effect
# ---------------------------------------------------
icc_overall_summary <- icc_draws %>%
  group_by(Component, Effect) %>%
  summarise(
    Mean_ICC = mean(Draws, na.rm = TRUE),
    SD_ICC = sd(Draws, na.rm = TRUE),
    .groups = "drop"
  )

# ---------------------------------------------------
# Generate Barplot with Individual Model Estimates
# ---------------------------------------------------
icc_plot <- ggplot() +
  geom_bar(data = icc_overall_summary, aes(x = Component, y = Mean_ICC, fill = Effect),
           stat = "identity", position = "dodge", alpha = 0.6) +
  geom_errorbar(data = icc_overall_summary, aes(x = Component, ymin = Mean_ICC - SD_ICC, ymax = Mean_ICC + SD_ICC, fill = Effect),
                position = position_dodge(0.9), width = 0.2, color = "black") +
  labs(x = "Component", y = "Mean ICC", title = "Mean ICC Across All Models with Individual Model Estimates") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("Species" = "blue", "Phylogeny" = "red")) +
  scale_fill_manual(values = c("Species" = "blue", "Phylogeny" = "red")) +
  ylim(0,1)

ggsave(plot = icc_plot, filename = paste0(model_dir, "../Figures/ICC_Bacteria_barplot.pdf"), device = "pdf", scale = 0.7)

print("ICC values calculated, saved, and plotted.")
print(icc_plot)


# ---------------------------------------------------
# Assess changes in WAIC
# ---------------------------------------------------
waic_results <- read.csv(paste0(model_dir, "Bacteria_WAIC_R2_results.csv"))

waic_summary <- waic_results %>%
  mutate(Model_Type = case_when(
    grepl("Species", Model) & grepl("Phylogeny", Model) ~ "Species + Phylogeny",
    grepl("Species", Model) & !grepl("Phylogeny", Model) ~ "Species Only",
    grepl("Phylogeny", Model) & !grepl("Species", Model) ~ "Phylogeny Only",
    TRUE ~ "Other"
  )) %>%
  group_by(Model_Type) %>%
  summarise(
    Mean_WAIC = mean(WAIC),
    SD_WAIC = sd(WAIC),
    Mean_R2 = mean(Bayes_R2_Mean),
    SD_R2 = sd(Bayes_R2_Mean),
    .groups = "drop"
  )

# Save summary
write.csv(waic_summary, file = paste0(model_dir, "Fungi_WAIC_R2_ModelType_summary.csv"), row.names = FALSE)
