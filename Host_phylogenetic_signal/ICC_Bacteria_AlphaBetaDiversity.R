# ---------------------------------------------------
# Load Required Libraries
# ---------------------------------------------------
library(phyloseq)
library(brms)
library(ape)
library(dplyr)
library(vegan)
library(progress)
library(ggplot2)

# ---------------------------------------------------
# Load Data and Prepare
# ---------------------------------------------------
load("VU16s_roottraits_updated.Rdata")
VU16s <- VU_16s_root_architecture_chem_dry

tree <- read.nexus("Phylogeny/host_phylogeny_microp_legume_guest.nex")

# Subset metadata
phyloseq_obj <- subset_samples(VU16s, Control != "yes")

# Estimate richness and log-transform it
richness <- estimate_richness(phyloseq_obj, measures = "Chao1")
metadata <- data.frame(sample_data(phyloseq_obj))
metadata$Species <- metadata$PlantShort
metadata$Chao1 <- log1p(richness$Chao1)  # Log-transform richness

# Prepare phylogenetic covariance matrix
metadata$Phylogeny <- metadata$Plant_Long
metadata$Species <- metadata$Plant_Long
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, metadata$Species))
phylo_cov <- ape::vcv.phylo(tree_pruned)

# ---------------------------------------------------
# Perform PCoA on Bray-Curtis Dissimilarity Matrix
# ---------------------------------------------------
otu_table_raw <- as(otu_table(phyloseq_obj), "matrix")
bc_dissim <- vegdist(t(otu_table_raw), method = "bray")

pcoa_results <- ape::pcoa(bc_dissim)

metadata$PC1 <- pcoa_results$vectors[, 1]
metadata$PC2 <- pcoa_results$vectors[, 2]

variance_explained_pcoa <- round(pcoa_results$values$Relative_eig[1:2] * 100, 1)

# Print variance explained for PC1 and PC2
print(data.frame(
  Axis = paste0("PC", 1:2),
  Variance_Explained = variance_explained_pcoa
))

# ---------------------------------------------------
# Define Model Structures
# ---------------------------------------------------
model_structures <- list(
  "P" = "~ Phosphorus",
  "Species" = "~ (1 | Species)",
  "Greenhouse" = "~ (1 | Greenhouse_Block)",
  "Phylogeny" = "~ (1 | gr(Phylogeny, cov = phylo_cov))",
  "P_Greenhouse" = "~ Phosphorus + (1 | Greenhouse_Block)",
  "Species_P" = "~ Phosphorus + (1 | Species)",
  "Species_Greenhouse" = "~ (1 | Species) + (1 | Greenhouse_Block)",
  "Species_P_Greenhouse" = "~ Phosphorus + (1 | Species) + (1 | Greenhouse_Block)",
  "Phylogeny_P" = "~ Phosphorus + (1 | gr(Phylogeny, cov = phylo_cov))",
  "Species_Phylogeny" = "~ (1 | Species) + (1 | gr(Phylogeny, cov = phylo_cov))",
  "Species_Phylogeny_P" = "~ Phosphorus + (1 | Species) + (1 | gr(Phylogeny, cov = phylo_cov))",
  "Phylogeny_Greenhouse" = "~ (1 | gr(Phylogeny, cov = phylo_cov)) + (1 | Greenhouse_Block)",
  "Phylogeny_P_Greenhouse" = "~ Phosphorus + (1 | gr(Phylogeny, cov = phylo_cov)) + (1 | Greenhouse_Block)",
  "Phylogeny_Species_Greenhouse" = "~ (1 | gr(Phylogeny, cov = phylo_cov)) + (1 | Species) + (1 | Greenhouse_Block)",
  "Phylogeny_Species_Greenhouse_P" = "~ Phosphorus + (1 | gr(Phylogeny, cov = phylo_cov)) + (1 | Species) + (1 | Greenhouse_Block)"
)

# ---------------------------------------------------
# Fit Models for PC1, PC2, and Richness in a Loop with Progress Bar
# ---------------------------------------------------
components <- c("PC1", "PC2", "Chao1")  # Richness added
model_results <- list()

for (i in seq_along(components)) {
  component <- components[i]
  variance <- ifelse(component %in% c("PC1", "PC2"), variance_explained_pcoa[i], "Richness")  # Use variance for PCoA, otherwise "Richness"
  
  total_models <- length(model_structures)
  pb <- progress_bar$new(
    format = paste0("  Fitting Models (", component, ") [:bar] :percent ETA: :eta"),
    total = total_models, clear = FALSE, width = 60
  )
  
  model_list <- list()
  for (model_name in names(model_structures)) {
    formula_str <- paste0(component, " ", model_structures[[model_name]])
    
    model_formula <- as.formula(formula_str)
    
    use_phylo_cov <- grepl("Phylogeny", formula_str)
    
    model <- brm(
      formula = model_formula,
      data = metadata,
      family = gaussian(),
      data2 = if (use_phylo_cov) list(phylo_cov = phylo_cov) else NULL,
      cores = 4, chains = 4, iter = 5000, warmup = 500, silent = TRUE
    )
    
    model_list[[model_name]] <- model
    pb$tick()
  }
  
  # Store models in the results list
  model_results[[component]] <- model_list
  
  # Save models for this component with variance explained in the filename
  save(
    model_list, 
    file = paste0("/Users/justinstewart/Dropbox/Collaborations/FamilyExperiment/Data/Code/New/Feb2024/Host_phylogenetic_signal/Models/", component, "_Bacteria_", variance, "variance.RData")
  )
}

print("All models fitted and saved with variance explained in filenames.")
