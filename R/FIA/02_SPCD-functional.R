#'------------------------------------------------------------------------------
#' @title A script to get functional species attributes
#'------------------------------------------------------------------------------

#' @description It compiles functional attributes of species from different
#' sources (FIA database, Niinemets & Valladares, and BIEN), and imputes
#' missing traits using phylogenetic information (V.PhyloMaker2, PVR, and
#' missForest). Some steps require a manual cleaning of the observations.
#'
#' @return .csv files with the functional traits per species (SPCD) after each
#' step of the compilation, being 06_fia-traits.csv the final imputed traits.

#'------------------------------------------------------------------------------
#' @Library

library(data.table)
library(phytools)
library(V.PhyloMaker2)
library(PVR)
library(missForest)
library(bit64)
library(BIEN)

#'------------------------------------------------------------------------------
#' @Paths-and-files

root_path <- "/media/antonio/Work/MMDB/FIA"
fia <- fread(paste0(root_path, "/00_FIA/trees_eastern_20230427.csv"))
SPCD <- fread(paste0(root_path, "/00_FIA/SPCD.csv"), na.strings="NA")
domain <- fread(paste0(root_path, "/00_FIA/plot_domain.csv"))

#'------------------------------------------------------------------------------
#' @Get-traits-from-the-FIA-database

fia_traits <- function(fia,
                       SPCD) {
  
  ### Prepare plots (Cleaning) -------------------------------------------------
  
  # Read fia
  fia <- fia[order(PLT_CN, SUBP)]
  
  # Inventories from 2018 until now
  fia <- fia[MEASYEAR >= 2018]
  
  # Remove dead trees
  fia <- fia[STATUSCD == 1]
  
  # Get forest metrics by TPA
  fia$BA <- pi*(fia$DIA/2)^2
  fia$BA_TPA <- fia$BA*fia$TPA_UNADJ
  fia$DRYBIO_AG_TPA <- fia$DRYBIO_AG*fia$TPA_UNADJ
  fia$VOL <- pi*(round(((fia$DIA*0.0833333333)/2), 6)^2)*fia$HT
  fia[VOL >= 600, VOL := NA,] #Potential errors
  fia$VOL_TPA <- fia$VOL*fia$TPA_UNADJ
  
  # Select plots with complete volume information
  complete_plots <- fia[, sum(is.na(VOL_TPA)), by = "PLT_CN"]
  complete_plots <- complete_plots[V1 == 0]
  
  fia <- merge(complete_plots[, 1], 
               fia, 
               by = "PLT_CN", 
               all.x = TRUE, 
               all.y = FALSE)
  
  # Select plots with complete tree identification
  fia <- merge(fia, 
               SPCD[, c("SPCD", "Taxon2")], 
               by = "SPCD", 
               all.x = TRUE, 
               all.y = TRUE)
  
  fia[is.na(Taxon2), Remove := 1]
  fia[is.na(Remove), Remove := 0]
  
  complete_plots <- fia[, sum(Remove), by = "PLT_CN"]
  complete_plots <- complete_plots[V1 == 0]
  
  fia <- merge(complete_plots[, 1], 
               fia, 
               by = "PLT_CN", 
               all.x = TRUE, 
               all.y = FALSE)
  
  # Get traits 
  fia$DIA <- fia$DIA*0.0254
  fia$HT <- fia$HT*0.3048
  fia$SL <- fia$HT/fia$DIA
  
  # Loop by SPCD
  
  unique_SPCD <- unique(fia$SPCD)
  traits <- data.table()
  
  for(i in 1:length(unique_SPCD)) {
    
    sub_fia <- fia[SPCD == unique_SPCD[i],]
    HTq <- quantile(sub_fia$HT, 0.95)
    SLq <- quantile(sub_fia$SL, 0.50)
    
    sub <- data.table(SPCD = unique_SPCD[i],
                      HTq = HTq,
                      SLq = SLq)
    
    traits <- rbind(traits, sub)
    
  }
  
  
  # Export 
  return(traits)
  
}

traits <- fia_traits(fia, SPCD)
traits <- merge(SPCD, traits, by = "SPCD", all.x = TRUE, all.y = TRUE)

### Export for manual cleaning of observations
fwrite(traits, paste0(root_path, "/02_FIA-traits/00_fia-traits.csv")) 

#'------------------------------------------------------------------------------
#' @Merge-traits-from-the-existing-FIA-database

traits <- fread(paste0(root_path, "/02_FIA-traits/00_fia-traits.csv")) 
res <- fread(paste0(root_path, "/00_FIA/REF_SPECIES.csv"))
colnames(res)[2] <- "WD"
traits <- merge(traits, res[, 1:2], by = "SPCD", all.x = TRUE, all.y = FALSE)

### Export for manual cleaning of observations
fwrite(traits, paste0(root_path, "/02_FIA-traits/01_fia-traits.csv")) 

#'------------------------------------------------------------------------------
#' @Merge-other-traits-from-Niinemets-and-Valladares

traits <- fread(paste0(root_path, "/02_FIA-traits/01_fia-traits.csv")) 
NV <- fread(paste0(root_path, "/02_FIA-traits/Niinemets_Valladares.csv"))
NV <- NV[, c("Species", "shade", "drought", "waterlogging")]
traits <- merge(traits, NV, by.x = "Taxon", by.y = "Species", all.x = TRUE, all.y = FALSE)

### Manual cleaning of observations
fwrite(traits, paste0(root_path, "/02_FIA-traits/03_fia-traits.csv")) 

#'------------------------------------------------------------------------------
#' @Get-traits-from-BIEN

traits <- fread(paste0(root_path, "/02_FIA-traits/03_fia-traits.csv"))

# Download BIEN (date: 2024-02-07)
SLA <- BIEN_trait_mean(trait = "leaf area per leaf dry mass", species = traits$Taxon)
N <- BIEN_trait_mean(trait = "leaf nitrogen content per leaf dry mass", species = traits$Taxon)

# Export BIEN records
bien <- rbind(SLA, N)
fwrite(bien, paste0(root_path, "/02_FIA-traits/FIA_BIEN.csv")) 

# Merge BIEN traits with existing traits
SLA <- as.data.table(SLA)
#SLA <- SLA[level_used == "Species", ]
SLA <- SLA[, 1:2]
LMA <- SLA
LMA$mean_value <- as.numeric(LMA$mean_value)
LMA$mean_value <- as.numeric((1/LMA$mean_value)*1000)
colnames(LMA) <- c("Taxon", "LMA")
traits <- merge(traits, LMA, by = "Taxon", all.x = TRUE, all.y = FALSE)

N <- as.data.table(N)
#N <- N[level_used == "Species", ]
N <- N[, 1:2]
colnames(N) <- c("Taxon", "N")
traits <- merge(traits, N, by = "Taxon", all.x = TRUE, all.y = FALSE)

### Manual cleaning of observations
fwrite(traits, paste0(root_path, "/02_FIA-traits/04_fia-traits.csv"))

################################################################################
## Do a cheeking to the traits to detect potential outliers
################################################################################

#'------------------------------------------------------------------------------
#' @Create-the-phylogenetic-tree-to-impute-traits

# Read species
traits <- fread(paste0(root_path, "/02_FIA-traits/04_fia-traits.csv"))
hypotheses <- phylo.maker(traits[, c("Taxon2", "Genus", "Family")], scenarios = "S3")
phylo <- multi2di(hypotheses$scenario.3)
is.binary.phylo(phylo)
is.ultrametric(phylo)

# Decompose the phylogenetic tree
decompose <- PVRdecomp(phylo)
pc <- decompose@Eigen$vectors[,1:10]
pc <- as.data.table(pc)
pc$Taxon2 <- phylo$tip.label

traits_and_phylo <- merge(traits, pc, by = "Taxon2", all.x = TRUE, all.y = TRUE)
fwrite(traits_and_phylo, paste0(root_path, "/02_FIA-traits/05_fia-traits.csv")) 

# Imputation of traits
imputation_traits <- missForest(traits_and_phylo[, .SD, .SDcols = c(10:27)], 
                                maxiter = 500, ntree = 100, verbose = TRUE)
imputation_traits_and_phylo <- cbind(traits_and_phylo[, 1:9], round(imputation_traits$ximp[, 1:8], 4))
fwrite(imputation_traits_and_phylo, paste0(root_path, "/02_FIA-traits/06_fia-traits.csv"))

################################################################################
## Do a cheeking and insert other categorical traits manually selected
################################################################################