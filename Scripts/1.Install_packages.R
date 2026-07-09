## pakages use in our scripts--------------------------
pkgs <- c("data.table", "dplyr", "ggplot2", "corrplot", "ggpubr", "Gmisc", "openxlsx", "readxl", "paletteer", "MuMIn",
          "ggsci", "scales", "RColorBrewer", "gridExtra", "tidyr", "stringr", "colorspace", "cowplot",
          "ggbreak", "rstatix", "randomForest", "forestploter", "sem", "lavaan", "bruceR", "tidyverse", "patchwork",
          "ggridges", "dorothea", "OmnipathR", "decoupleR", "viridis", "ggpattern", "biomaRt", "dbplyr", 
          "GenomicRanges", "LOLA", "simpleCache", "ggbeeswarm",  
          "genekitr", "igraph", "ggraph", "rrvgo", "clusterProfiler", "enrichplot", 
          "janitor", "binom", "glue", "lawstat", "weights", "epitools", "DescTools", "optparse", "MASS","mediation", "org.Hs.eg.db",
          "GOSemSim", "DOSE", "tidygraph", "scatterpie", "e1071", "ggrepel", "broom", "purrr", "effectsize", "tibble",
          "STRINGdb", "ReactomePA", "forcats", "ggstatsplot", "progeny"，
         "future", "future.apply")

miss_pkgs <- setdiff(pkgs, installed.packages()[, "Package"])
bio_pkgs <- c("biomaRt", "GenomicRanges", "LOLA", "clusterProfiler", "enrichplot", 
              "dorothea", "OmnipathR", "decoupleR", "org.Hs.eg.db", "GOSemSim", "DOSE", "STRINGdb",
              "ReactomePA", "progeny")


## install the missing packages--------------------------
for(i in miss_pkgs){
  
  tryCatch({
    if(i %in% bio_pkgs) {
      
      # If it's a Bioconductor package, install using BiocManager
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(i)
    } else {
      
      # If it's a CRAN package, install using install.packages
      install.packages(i)
    }
  }, error = function(e) {
    
    cat("Error installing package:", i, "\n")
    message("Error message:", e$message, "\n")
  })
  
}
