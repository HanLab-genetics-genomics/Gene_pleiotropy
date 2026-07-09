
# A gene-level framework for statistical pleiotropy in humans reveals evolutionary architecture of pleiotropic genes and aging-related biology

## Overview
This repository contains the data and source code used to produce the manuscript’s main analysis. It includes:
### Scripts
Seven scripts (one per main figure) that generate the figures from the input data. R scripts and functions demonstrating how the analysis was performed.
* [Gene pleiotropy across the human genome](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Scripts/2.Gene_pleiotropy_across_human_genome.R): Gene-level pleiotropy scores reveal pervasive pleiotropy in the human genome.  
* [Gene pleiotropy and gene age](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Scripts/3.Gene_pleiotropy_and_gene_age.R): Gene evolutionary age is associated with pleiotropy, genes from older evolutionary groups exhibiting higher pleiotropy.  
* [Structural and evolutionary characteristics](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Scripts/4.Structural_evolutionary_characteristics.R): Gene structural and evolutionary characteristics across pleiotropy groups and age categories. 
* [Gene regulatory circuitry](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Scripts/5.Gene_regulatory_circuitry.R): Complex gene regulatory circuitry and interaction architectures in highly pleiotropic and evolutionarily older genes.
* [Epigenomic modifications](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Scripts/6.Epigenomic_modifications.R): Active epigenomic modifications in highly pleiotropic and evolutionarily older genes.
* [Functional characterization](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Scripts/7.Functional_characterization.R): Functional characterization of highly pleiotropic and evolutionarily older genes.
* [AgingGPS](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Scripts/8.AgingGPS.R): Aging-specific pleiotropy scores.
* [Duplication](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Scripts/9.Duplication.R): Gene pleiotropy and gene dupliation (supporting analyses).
* [Other](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Scripts/10.Other.R): Supporting analyses in Supplementary Notes

### Data:
Input datasets and processed files required to reproduce the results.
| Data description | Data source |
| :-------- | :-------: |
| **Datasets for [Gene pleiotropy across the human genome]** |  | 
| Pleiotropy main dataset |  [pleiotropy_maindata.RData](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/pleiotropy_maindata.RData) |
| Gene locations (NCBI 37.3) |  [NCBI37.3.gene.loc](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/NCBI37.3.gene.loc) |
| **Datasets for [Structural and evolutionary characteristics]** |  | 
| Tissue specific τ index and gene expression data |  [Tau_gene_V8.csv](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/Tau_gene_V8.csv) |
| Processed results of Random forest analysis (unrar first) |  [RandomForest](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/RandomForest) |
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
| UCSC features for LOLA analysis (unrar all file in the LOLACore directory) |  [LOLACore_UCSC features](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/LOLACore) |
| Roadmap epigenomic data for LOLA analysis (unrar all file in the LOLACore directory) |  [LOLACore_Roadmap epigenomics](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/LOLACore) |
| **Datasets for [Functional characterization]** |  | 
| Gene effect score (unrar first) |  [CRISPRGeneEffect](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/CRISPRGeneEffect) |
| Metabolic processes in Last Universal Common Ancestor (LUCA) |  [LUCA](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/41559_2024_2461_MOESM4_ESM.tsv) |
| **GPSs for sensitivity analysis** |  [GPS](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/GPS) |
| **GPS using a restricted panel of reproduction and pre-reproductive survival traits** |  [hopsgene_repro_cut.RData](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/hopsgene_repro_cut.RData) |
| **Datasets for [AgingGPS]** |  [AgingGPS](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/AgingGPS) |
| **Datasets for [Duplication]** |  [Duplication](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/Duplication) |
| **Datasets for [matched-gene permutation analyses]** |  [Matched-gene permutation analyses](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/permutation) |

### Output:
All figures are saved as PDF files with informative names that map directly to the manuscript figures (e.g. figure1_*, figure2_*) in [**Output.**](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Output/)


## System requirements
### Operating systems
Linux (Ubuntu 20.04), Windows 10/11, or macOS (12+ Monterey tested)
### Dependencies
R ≥4.3.1

The scripts were developed and tested in R (version 4.3.1).
### Hardware
No non-standard hardware required. A typical laptop/desktop (≥16 GB RAM recommended) is sufficient for these scripts.


##  Instructions for use
### Directory conventions
- All analysis scripts are stored in "/path/to/Gene_pleiotropy-main/Scripts".
- All input files are stored in "/path/to/Gene_pleiotropy-main/Data". The required input files are loaded from the "Data" directory.
- Each script writes PDF files to "/path/to/Gene_pleiotropy-main/Output/". 

### Installation guide and preparation before running the scripts:
- **Download the source code and data archive [Gene_pleiotropy-main.zip]**
- **Unzip the archive**
- **Change to the working directory and unrar the necessary files in the native directory as indicated in the Data Table above:**
```bash
wget -c https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/archive/refs/heads/main.zip
unzip main.zip
cd Gene_pleiotropy-main
unrar x Data/LOLACore/LOLACore.part01.rar Data/LOLACore
unrar x Data/RandomForest/forest_13var_seed123_250114.part01.rar Data/RandomForest
unrar x Data/mediation_list_pn.rar Data
unrar x Data/genehancer_all.rar Data
unrar x Data/CRISPRGeneEffect/CRISPRGeneEffect.part01.rar Data/CRISPRGeneEffect
unrar x Data/AgingGPS/AgingGPS.part1.rar Data/AgingGPS
unrar x Data/Duplication/Duplication.part01.rar Data/Duplication
```
- **Install the required R packages by running "[1.Install_packages.R](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Scripts/1.Install_packages.R)**
  
Typical install time on a “normal” desktop computer: approximately 0.5-1h, although the actual time will vary depending on your hardware and internet connection. 

- **Set the Working directory to the [**Data**](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/).**
  
Each script begins with a line such as: setwd("/path/to/Gene_pleiotropy-main/Data"), which you should modify to the path where you saved the Gene_pleiotropy-main/Data folder.

- **Set the Output directory to the [**Output**](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Output/).**
  
Each script also defines an output path, for example: figure_file <- "/path/to/Gene_pleiotropy-main/Output/Figure1". You can edit this to change the directory where output files will be saved.

## Running on your data
Set the working directory to the [**Data**](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/) folder at the top of each script [setwd("/path/to/Gene_pleiotropy-main/Data")].

- **Reproduce figures**: Run the scripts provided for each main figure. These scripts will load the relevant datasets from the Data folder and output the final plots.
- **Review the Scripts**: See the [**Scripts**](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Scripts/) folder for detailed scripts and function definitions.
- **Explore the Data**: Browse the [**Data**](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Data/) folder to view the input datasets and intermediate processed files.
- **Explore the Output**: All figures are saved as PDF files with informative names that map directly to the manuscript figures (e.g. figure1_*, figure2_*) in [**Output.**](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Output/)

Runtime: approximately 2-3 hours for all scripts combined.  
Note: Runtimes refer to a standard laptop/desktop with required packages installed and the provided input files.

## Notes for matched-gene permutation analyses
The matched-gene empirical null analyses are implemented in two scripts:

1. [11.1 run_matched_permutation.R](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Scripts/11.1%20run_matched_permutation.R) generates matched control matrices for GPS-N and GPS-M.
2. [11.2 matched_permutation_analysis.R](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Scripts/11.2%20matched_permutation_analysis.R) uses these matched control matrices to re-evaluate gene-level feature enrichments and generate the summary table.

The first step is computationally intensive because it performs matched sampling across the gene universe and 1,000 permutations. We tested this step on a Linux server using parallel workers, and we recommend running it as a background job on Linux. Additional running notes are provided at the beginning of [11.1 run_matched_permutation.R](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Scripts/11.1%20run_matched_permutation.R).

For convenience, we provide the output files generated by step 11.1: [matched_ctrl_mat_pn.rds](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Output/permutation/matched_ctrl_mat_pn.rds) and [matched_ctrl_mat_pm.rds](https://github.com/HanLab-genetics-genomics/Gene_pleiotropy/blob/main/Output/permutation/matched_ctrl_mat_pm.rds).

Therefore, users who want to reproduce the final matched-permutation summary table can skip step 11.1 and directly run step 11.2.

###  Reproduction instructions (end-to-end)
1. Install R (version ≥ 4.3.1) and all dependencies (CRAN + Bioconductor).
2. Place all required input files into "/path/to/Gene_pleiotropy-main/Data".
3. Run the scripts in "/path/to/Gene_pleiotropy-main/Scripts".
4. Collect outputs from "/path/to/Gene_pleiotropy-main/Output". Each file name maps to the figure referenced in the manuscript.

### Troubleshooting
Packages fail to install or load  
- Rerun the installation script; the built-in safe installers will retry and print a message identifying the problematic package.  
- Ensure that you have a working internet connection and that CRAN/Bioconductor mirrors are reachable.

If you encounter any issues, please reach out, and we will resolve them promptly.
