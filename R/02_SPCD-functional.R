#'------------------------------------------------------------------------------
#' @title A script to get functional species attributes
#'------------------------------------------------------------------------------

#' @description It compile a series of functional attributes of species
#' 
#' @return Return a frame as a .csv file

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(phytools)
library(V.PhyloMaker2)
library(PVR)
library(missForest)
library(bit64)
library(BIEN)

#'------------------------------------------------------------------------------
#' Paths and files

root_path <- "/media/antonio/Work/MMDB"
#root_path <- "E:/Projects/DESIS_mapping"
root_path <- "/media/antonio/work/Projects/DESIS_mapping"
fia <- fread(paste0(root_path, "/00_FIA/trees_eastern_20230427.csv"))
SPCD <- fread(paste0(root_path, "/00_FIA/SPCD.csv"), na.strings="NA")
FIA_SPCD <- fread(paste0(root_path, "/00_FIA/FIA_SPCD.csv"), na.strings="NA")
domain <- fread(paste0(root_path, "/00_FIA/plot_domain.csv"))


#'------------------------------------------------------------------------------
# Mode function
mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#'------------------------------------------------------------------------------
#' Get traits from the FIA database

fia_traits <- function(fia,
                       SPCD) {
  
  ### Prepare plots (Cleaning) -------------------------------------------------
  
  # Read fia
  fia <- fia[order(PLT_CN, SUBP)]
  
  # Inventories from 2018 until now
  fia <- fia[MEASYEAR >= 2018]
  
  # Merge species code
  fia <- merge(fia, 
               SPCD[, c("SPCD", "Remove")], 
               by = "SPCD", 
               all.x = TRUE, 
               all.y = TRUE)
  
  # Remove dead trees
  fia <- fia[STATUSCD == 1]
  
  # Get forest metrics by TPA
  fia$BA <- pi*(fia$DIA/2)^2
  fia$BA_TPA <- fia$BA*fia$TPA_UNADJ
  fia$DRYBIO_AG_TPA <- fia$DRYBIO_AG*fia$TPA_UNADJ
  fia$VOL <- pi*(round(((fia$DIA*0.0833333333)/2), 6)^2)*fia$HT
  fia[VOL >= 600, VOL := NA,] #Potential errors
  fia$VOL_TPA <- fia$VOL*fia$TPA_UNADJ
  
  # Get plots with complete data
  complete_plots <- fia[, sum(is.na(VOL_TPA)), by = "PLT_CN"]
  complete_plots <- complete_plots[V1 == 0]
  
  fia <- merge(complete_plots[, 1], 
               fia, 
               by = "PLT_CN", 
               all.x = TRUE, 
               all.y = FALSE)
  
  # Select plots with complete information
  fia[is.na(Remove), Remove := 1]
  
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
    HTq <- quantile(sub_fia$HT, 0.50)
    SLq <- quantile(sub_fia$SL, 0.50)
    sub_fia <- sub_fia[HT >= quantile(sub_fia$HT, 0.75), ]
    CCLCD_mode <- mode(sub_fia$CCLCD)
    
    sub <- data.table(SPCD = unique_SPCD[i],
                      HTq = HTq,
                      SLq = SLq,
                      CCLCD_mode = CCLCD_mode)
    
    traits <- rbind(traits, sub)
    
  }
  
  
  # Export 
  return(traits)
  
}

traits <- fia_traits(fia, SPCD)
traits <- merge(traits, FIA_SPCD, by = "SPCD", all.x = TRUE, all.y = FALSE)

### Export for manual cleaning of observations
fwrite(traits, paste0(root_path, "/00_FIA/traits/fia_traits.csv")) 

#'------------------------------------------------------------------------------
# Merge other traits from FIA
traits <- fread(paste0(root_path, "/00_FIA/traits/fia_traits.csv")) 
res <- fread(paste0(root_path, "/00_FIA/traits/REF_SPECIES.csv"))
traits <- merge(traits, res, by = "SPCD", all.x = TRUE, all.y = FALSE)

### Export for manual cleaning of observations
fwrite(traits, paste0(root_path, "/00_FIA/traits/fia_traits_res.csv")) 

#'------------------------------------------------------------------------------
# Merge other traits from Niinemets and Valladares
traits <- fread(paste0(root_path, "/00_FIA/traits/fia_traits_res.csv")) 
NV <- fread(paste0(root_path, "/00_FIA/traits/Niinemets_Valladares.csv"))
NV <- NV[, c("Species", "shade", "drought", "waterlogging")]
traits <- merge(traits, NV, by.x = "Taxon", by.y = "Species", all.x = TRUE, all.y = FALSE)

### Manual cleaning of observations
fwrite(traits, paste0(root_path, "/00_FIA/traits/fia_traits_res_nv.csv")) 

#'------------------------------------------------------------------------------
# Get traits from BIEN (Botanical Information and Ecology Network)

traits <- fread(paste0(root_path, "/00_FIA/traits/fia_traits_res_nv.csv"))

# Download BIEN (date: 2024-02-07)
leaf_area <- BIEN_trait_mean(trait = "leaf area", species = traits$Taxon)
SLA <- BIEN_trait_mean(trait = "leaf area per leaf dry mass", species = traits$Taxon)
N <- BIEN_trait_mean(trait = "leaf nitrogen content per leaf dry mass", species = traits$Taxon)

# Export BIEN records
bien <- rbind(leaf_area, SLA, N)
fwrite(bien, paste0(root_path, "/00_FIA/traits/FIA_BIEN.csv")) 

# Merge BIEN traits with existing traits
leaf_area <- as.data.table(leaf_area)
#leaf_area <- leaf_area[level_used == "Species", ]
leaf_area <- leaf_area[, 1:2]
colnames(leaf_area) <- c("Taxon", "leaf_area")
traits <- merge(traits, leaf_area, by = "Taxon", all.x = TRUE, all.y = FALSE)

SLA <- as.data.table(SLA)
#SLA <- SLA[level_used == "Species", ]
SLA <- SLA[, 1:2]
colnames(SLA) <- c("Taxon", "SLA")
traits <- merge(traits, SLA, by = "Taxon", all.x = TRUE, all.y = FALSE)

N <- as.data.table(N)
#N <- N[level_used == "Species", ]
N <- N[, 1:2]
colnames(N) <- c("Taxon", "N")
traits <- merge(traits, N, by = "Taxon", all.x = TRUE, all.y = FALSE)

### Manual cleaning of observations
fwrite(traits, paste0(root_path, "/00_FIA/traits/fia_traits_res_nv_bien.csv")) 

#'------------------------------------------------------------------------------
# Create the phylogenetic tree

# Read species
traits <- fread(paste0(root_path, "/00_FIA/traits/fia_traits_res_nv_bien.csv"))
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
fwrite(traits_and_phylo, paste0(root_path, "/00_FIA/traits/traits_and_phylogeny.csv")) 

# Imputation of traits
imputation_traits <- missForest(traits_and_phylo[, .SD, .SDcols = c(13, 14, 16:34)], 
                                maxiter = 500, ntree = 100, verbose = TRUE)
imputation_traits_and_phylo <- cbind(traits_and_phylo[, 1:12], imputation_traits$ximp)
fwrite(imputation_traits_and_phylo, paste0(root_path, "/00_FIA/traits/imputation_traits_and_phylo (2024-10-08).csv")) 
