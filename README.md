Plant Rhizosphere Microbiomes and Root Architecture

This repository contains data and analysis workflows examining how plant root architecture and host identity shape rhizosphere microbial community composition. 
The project integrates root phenotyping with bacterial and fungal community profiling to test the relative importance of plant traits, species identity, and phylogenetic relatedness in structuring rhizosphere microbiomes.

This sutdy is currently in review. 

.
├── BetaDiversity/
│   PCoA, beta diversity metrics, PERMANOVA, and Random Forest models
│   Includes SHAP analyses for bacterial and fungal communities
│
│   ├── Pcoa_BetaDiversity.R
│   ├── RF_BetaDiversity.R 
│   ├── RF_BetaDiversity_SHAP_Bacteria.R
│   └── RF_BetaDiversity_SHAP_Fungi.R
│
├── Data/
│   Input data objects used across analyses
│   Includes microbiome data, root trait data, and experimental metadata
│
│   ├── VU16s_roottraits_updated.RData
│   ├── VU16s_w_control.RData
│   └── VU_ITS_roottraits_updated.RData
│
├── Domestication/
│   Analyses comparing domesticated plants and wild relatives
│   Focuses on both root traits and rhizosphere microbiomes
│
│   ├── Domestication.R
│   ├── Domestication_Roots.R
│   └── Pcoa_Domestication.R
│
├── Host_phylogenetic_signal/
│   Tests for phylogenetic signal in microbiome composition and model residuals
│
│   ├── Models/
│   │   └── Bacteria/
│   ├── Hypotheses_ICC_Bacteria_AlphaBetaDiversity.R
│   └── ICC_Bacteria_AlphaBetaDiversity.R
│
├── Root examples/
│   Example root images and reference material for root architecture
│
└── README.md
