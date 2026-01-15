Plant Rhizosphere Microbiomes and Root Architecture (In Review)

This repository contains data and analysis workflows examining how plant root architecture and host identity shape rhizosphere microbial community composition. 
The project integrates root phenotyping with bacterial and fungal community profiling to test the relative importance of plant traits, species identity, and phylogenetic relatedness in structuring rhizosphere microbiomes.

BetaDiversity/
Contains scripts for analyzing rhizosphere microbial beta diversity, including PCoA ordinations, PERMANOVA, and Random Forest models. This folder also includes SHAP analyses to interpret bacterial and fungal machine learning models.

Data/
Houses processed data objects used across all analyses, including bacterial 16S and fungal ITS community data, root architectural trait measurements, and associated experimental metadata.

Domestication/
Includes analyses comparing domesticated plant species with their wild relatives, focusing on differences in root architecture and rhizosphere microbiome composition.

Host_phylogenetic_signal/
Contains scripts testing for phylogenetic signal in rhizosphere microbiome composition and in model residuals, including ICC based and mixed model approaches.

Root examples/
Provides example root images and reference material used to illustrate and validate root architecture traits.
