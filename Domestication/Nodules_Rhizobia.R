# =============================================================================
# Libraries
# =============================================================================
library(phyloseq)
library(speedyseq)
library(dplyr)
library(ggplot2)
library(forcats)
library(patchwork)
library(tibble)

# =============================================================================
# Settings
# =============================================================================
save_dir <- "/"
dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Load data
# =============================================================================
load("VU16s_roottraits_updated.Rdata")

ps <- VU_16s_root_architecture_chem_dry
ps <- prune_samples(sample_sums(ps) > 0, ps)

meta <- data.frame(sample_data(ps)) %>%
  rownames_to_column("Sample")

# =============================================================================
# Rhizobial genera
# =============================================================================
rhizobia_genera <- c(
  "Rhizobium",
  "Bradyrhizobium",
  "Sinorhizobium",
  "Ensifer",
  "Mesorhizobium",
  "Azorhizobium",
  "Allorhizobium",
  "Neorhizobium",
  "Pararhizobium",
  "Cupriavidus",
  "Burkholderia"
)

# =============================================================================
# Transform to relative abundance
# =============================================================================
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# =============================================================================
# Melt once with speedyseq
# =============================================================================
df <- speedyseq::psmelt(ps_rel)

# =============================================================================
# Rhizobia summary by plant species and genus
# =============================================================================
rhizobia_by_plant <- df %>%
  filter(!is.na(Genus), Genus %in% rhizobia_genera) %>%
  group_by(PlantShort, Genus) %>%
  summarise(
    MeanRelAbundance = mean(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

# =============================================================================
# Rhizobia plot
# =============================================================================
p_rhizobia <- ggplot(
  rhizobia_by_plant,
  aes(
    x = fct_reorder(PlantShort, MeanRelAbundance, .fun = sum),
    y = MeanRelAbundance,
    fill = Genus
  )
) +
  geom_col() +
  coord_flip() +
  theme_classic(base_size = 14) +
  labs(
    title = "Rhizobial genera across plant species",
    x = "Plant species",
    y = "Mean relative abundance"
  )

print(p_rhizobia)

# =============================================================================
# Nodulation summary for Fabaceae only
# Ignore NAs and count only Yes replicates
# =============================================================================
nodule_summary_legumes <- meta %>%
  filter(Plant_Family == "Fabaceae", !is.na(Nodules)) %>%
  group_by(PlantShort) %>%
  summarise(
    NodulatedReplicates = sum(Nodules == "Yes"),
    .groups = "drop"
  )

print(nodule_summary_legumes)

# =============================================================================
# Nodulation plot
# =============================================================================
p_nodules <- ggplot(
  nodule_summary_legumes,
  aes(
    x = fct_reorder(PlantShort, NodulatedReplicates),
    y = NodulatedReplicates
  )
) +
  geom_col(fill = "black") +
  coord_flip() +
  theme_classic(base_size = 14) +
  labs(
    title = "Nodulation across Fabaceae species",
    x = "Plant species",
    y = "Number of nodulated replicates"
  )

print(p_nodules)

# =============================================================================
# Combined figure
# =============================================================================
combined <- p_rhizobia / p_nodules +
  plot_layout(heights = c(1.2, 1))

print(combined)

# =============================================================================
# Save outputs
# =============================================================================
write.csv(
  rhizobia_by_plant,
  file.path(save_dir, "rhizobia_by_plant.csv"),
  row.names = FALSE
)

write.csv(
  nodule_summary_legumes,
  file.path(save_dir, "nodule_summary_fabaceae.csv"),
  row.names = FALSE
)

ggsave(
  filename = file.path(save_dir, "rhizobia_plot.pdf"),
  plot = p_rhizobia,
  width = 9,
  height = 6
)

ggsave(
  filename = file.path(save_dir, "nodules_plot.pdf"),
  plot = p_nodules,
  width = 7,
  height = 6
)

ggsave(
  filename = file.path(save_dir, "combined.pdf"),
  plot = combined,
  width = 9,
  height = 10
)