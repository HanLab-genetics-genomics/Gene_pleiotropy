
# Tracing time’s imprint: unveiling the pervasive pleiotropy of evolutionarily old genes in human diseases and traits

## Overview
This repository contains the data and source code used to produce the manuscript’s main analysis. It includes:
### Scripts
Seven scripts (one per main figure) that generate the figures from the input data. R scripts and functions demonstrating how the analysis was performed.
* [Gene pleiotropy across the human genome](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Scripts/Gene_pleiotropy_across_human_genome.R): Gene-level pleiotropy scores reveal pervasive pleiotropy in the human genome.  
* [Gene pleiotropy and gene age](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Scripts/Gene_pleiotropy_and_gene_age.R): Gene evolutionary age is associated with pleiotropy, genes from older evolutionary groups exhibiting higher pleiotropy.  
* [Structural and evolutionary characteristics](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Scripts/Structural_evolutionary_characteristics.R): Gene structural and evolutionary characteristics across pleiotropy groups and age categories. 
* [Gene regulatory circuitry](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Scripts/Gene_regulatory_circuitry.R): Complex gene regulatory circuitry and interaction architectures in highly pleiotropic and evolutionarily older genes.
* [Epigenomic modifications](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Scripts/Epigenomic_modifications.R): Active epigenomic modifications in highly pleiotropic and evolutionarily older genes.
* [Functional characterization](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Scripts/Functional_characterization.R): Functional characterization of highly pleiotropic and evolutionarily older genes.
* [Drug development](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Scripts/Drug_development.R): Pleiotropy evidence contributes to drug development.

### Data:
Input datasets and processed files required to reproduce the results.
| Data description | Data source |
| :-------- | :-------: |
| **Datasets for [Gene pleiotropy across the human genome]** |  | 
| Pleiotropy main dataset |  [pleiotropy_maindata.RData](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/pleiotropy_maindata.RData) |
| Gene locations (NCBI 37.3) |  [NCBI37.3.gene.loc](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/NCBI37.3.gene.loc) |
| **Datasets for [Structural and evolutionary characteristics]** |  | 
| Tissue specific τ index and gene expression data |  [Tau_gene_V8.csv](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/Tau_gene_V8.csv) |
| Processed results of Random forest analysis (unrar first) |  [Randomforest](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/forest_13var_seed123_250114.RData) |
| Processed results of single mediation analysis (unrar first) |  [Single_mediation.rar](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/mediation_list_pn.rar) |
| Processed results of multiple mediation analysis |  [Multiple_mediation.RData](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/lavaan_multiple_mediation_singlesig_pn.RData) |
| **Datasets for [Gene regulatory circuitry]** |  | 
| Proportion data of cis-eGenes and cis-sGenes |  [cis_esqtl.RData](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/cis_esqtl.RData) |
| Trans-eQTL data |  [Trans-eQTL](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/eqtlgen_transeqtl.txt) |
| DNA-binding RNA polymerase II TFs |  [TFC2_16102023b.tsv](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/TFC2_16102023b.tsv) |
| Enhancer data (unrar first) |  [Enhancer](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/genehancer_all.rar) |
| Enhancer annotaion data |  [Enhancer_annotation](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/genehancer_annotation.csv) |
| Protein-protein interactions |  [connectivity.txt](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/connectivity.txt) |
| **Datasets for [Epigenomic modifications]** |  |
| Gene region coordinates |  [bed_data](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/bed_data) |
| Chromatin 3D freature source |  [three_d_data.xlsx](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/three_d_data.xlsx) |
| Chromatin 3D freature dataset for LOLA analysis |  [LOLA3D](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/LOLA3D) |
| UCSC features for LOLA analysis (download first) |  [LOLACore - UCSC features](http://big.databio.org/regiondb/LOLACoreCaches_180412.tgz) |
| Roadmap epigenomic data for LOLA analysis (download first) |  [LOLACore - Roadmap epigenomics](http://big.databio.org/regiondb/LOLARoadmap_180423.tgz) |
| **Datasets for [Functional characterization]** |  | 
| Gene effect score (unrar first) |  [CRISPRGeneEffect](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/CRISPRGeneEffect) |
| Metabolic processes in Last Universal Common Ancestor (LUCA) |  [LUCA](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/41559_2024_2461_MOESM4_ESM.tsv) |
| **Datasets for [Drug development]** |  | 
| FDA approved drug targets |  [FDA approved drug targets](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/fda_2023.xlsx) |
| Data for durg Target-Indication pairs |  [genetic_support-main](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/genetic_support-main) |

## Usage
- **Reproduce figures**: Run the scripts provided for each main figure. These scripts will load the relevant datasets from the Data folder and output the final plots. You may need to revise the working directory using the R function setwd().
- **Review the Scripts**: See the [**Scripts**](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Scripts/) folder for detailed scripts and function definitions.
- **Explore the Data**: Browse the [**Data**](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/) folder to view the input datasets and intermediate processed files.

## Dependencies
The scripts were developed and tested in R (version 4.3.1).
