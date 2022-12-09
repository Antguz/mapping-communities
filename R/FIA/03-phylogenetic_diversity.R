#'------------------------------------------------------------------------------
#' @title Metrics phylogenetic diversity
#'------------------------------------------------------------------------------

#' @description It estimates metrics of the phylogenetic diversity
#' @return .txt file with the stack of results

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(picante)
library(Taxonstand)
library(V.PhyloMaker)
library(dplyr)
library(tidyr)
library(phytools)

#'------------------------------------------------------------------------------
#' Source code
source("https://raw.githubusercontent.com/jesusNPL/BiodiversityScience/master/Spring2021/R-Functions/qDp.R")

#'------------------------------------------------------------------------------
#' Arguments
#' @param master
#' @param ID_species
#' @param out_path
#' @param threads

root_path <- "F:/DESIS"
root_path <- "/media/antonio/Extreme SSD/DESIS"
fia <- fread(paste0(root_path, "/00_FIA/trees_WI_MN_20221102.csv"))
SPCD <- fread(paste0(root_path, "/00_FIA/SPCD.csv"), na.strings="NA")
out_path <- paste0(root_path, "/03_diversity")
threads <- 30

#'------------------------------------------------------------------------------
#' @example 

#'------------------------------------------------------------------------------
#' Function

phylogenetic_diversity <- function(master, threads) {
  
  #Prepare plots----------------------------------------------------------------
  #Clean FIA
  fia <- subset(fia, STATECD == 27 | STATECD == 55)
  fia <- fia[order(PLT_CN, SUBP)]
  
  #Merge species code
  fia <- merge(fia, SPCD[, c("SPCD", "PLANTS-Code", "Taxon2")], by = "SPCD", all.x = TRUE, all.y = TRUE)
  
  #Get Basal Area
  fia$BA <- pi*(fia$DIA/2)^2
  
  #Remove dead trees
  fia <- fia[is.na(STANDING_DEAD_CD)]
  
  #Summarize
  fia_summary <- fia[, .(ntrees = .N,
                         BA = sum(BA, na.rm = TRUE), 
                         HT = mean(HT, na.rm = TRUE), 
                         TPA_UNADJ = sum(TPA_UNADJ, na.rm = TRUE), 
                         DRYBIO_BG = sum(DRYBIO_BG, na.rm = TRUE),
                         CARBON_AG = sum(CARBON_AG, na.rm = TRUE), 
                         DRYBIO_AG = sum(DRYBIO_AG, na.rm = TRUE), 
                         CARBON_BG = sum(CARBON_BG, na.rm = TRUE)), 
                     by = c("PLT_CN", "PLANTS-Code", "Taxon2")]
  fia_summary <- fia_summary[order(PLT_CN)]
  
  #Check taxonomic
  #ID_species <- paste0(SPCD$Genus, " ", SPCD$Species)
  #ID_species <- ID_species[SPCD$FLAG != 1]
  #spp_check <- TPL(ID_species, infra = TRUE, corr = TRUE)
  
  #taxonomy <- spp_check %>% 
  #  select(Taxon, Family, New.Genus, New.Species) #, Tax_res
  
  #taxonomy$SP <- ID_species #Provide the ID
  #colnames(taxonomy) <- c("species", "family", "genus", "sp", "SP")
  #taxonomy <- taxonomy[, c("SP", "species", "genus", "sp", "family")]
  #taxonomy$SPP <- paste0(taxonomy$genus, "_", taxonomy$sp)
  #taxonomy$species[12] <- "Populus_X_deltoides"
  #taxonomy$SPP[12] <- "Populus_X_deltoides"
  
  #Prepare phylo  --------------------------------------------------------------
  #Makes phylogenetic hypotheses
  sppPhylo <- data.table(species = SPCD$Taxon2, genus = SPCD$Genus, family = SPCD$Family)
  sppPhylo <- na.exclude(sppPhylo)
  hypotheses <- phylo.maker(sppPhylo, scenarios = "S3")
  
  #Resolve Multichotomies
  phylo <- multi2di(hypotheses$scenario.3)
  is.binary.phylo(phylo) #Test for Binary Tree
  is.ultrametric(phylo) #Test if a Tree is Ultrametric
  
  #Create matrix
  master_matrix <- sample2matrix(fia_summary[, c("PLT_CN", "BA", "Taxon2")])
  master_matrix <- round(master_matrix, 2)
  
  matched <- picante::match.phylo.comm(phy = phylo, comm = master_matrix)
  matched$phy
  matched$comm
  
  MDP <- comdist(matched$comm, cophenetic(matched$phy), abundance.weighted=TRUE)
  
  #Get Phylogenetic diversity
  PD <- pd(matched$comm, matched$phy, include.root = TRUE) # Faith's PD
  
  #Get Mean pairwise distance (MPD) and mean nearest-pairwise distance (MNTD)
  MPD <- mpd(matched$comm, cophenetic(matched$phy)) # MPD
  MNTD <- mntd(matched$comm, cophenetic(matched$phy)) # MNTD
  
  #Phylogenetic diversity in a community - PD
  CPD <- ses.pd(matched$comm, matched$phy, runs = 99)
  CPD <- CPD[, c(1, 2, 6, 7)]
  
  #Add rownames
  CPD$STAND_KEY <- rownames(CPD)
  CPD <- CPD[, c(5, 1, 2, 3, 4)]
  CPD <- as.data.table(CPD)
  
  #Phylogenetic Rao’s quadratic entropy - RaoD
  CPD$RaoD <- raoD(matched$comm, force.ultrametric(matched$phy))$Dkk
  
  #Mean pairwise distance separating taxa in a community - MPD
  SESMPD <- ses.mpd(matched$comm, cophenetic(matched$phy), null.model = "taxa.labels", runs = 99)
  CPD$mpd <- SESMPD[, c(2)]
  CPD$mpd.obs.z <- SESMPD[, c(6)]
  CPD$mpd.obs.p <- SESMPD[, c(7)]
  
  #Mean nearest taxon distance for taxa in a community - MNTD
  SESMNTD <- ses.mntd(matched$comm, cophenetic(matched$phy), 
                      null.model = "taxa.labels", runs = 99)
  CPD$mntd <- SESMNTD[, c(2)]
  CPD$mntd.obs.z <- SESMNTD[, c(6)]
  CPD$mntd.obs.p <- SESMNTD[, c(7)]
  
  #Phylogenetic species variability - PSV
  PSV <- psv(matched$comm, matched$phy, compute.var = TRUE)
  CPD$PSV <- PSV[, 1]
  
  #Phylogenetic species richness - PSR
  PSR <- psr(matched$comm, matched$phy, compute.var = TRUE)
  CPD$PSR <- PSR[, 1]
  
  #Phylogenetic species evenness - PSE
  PSE <- pse(matched$comm, matched$phy)
  CPD$PSE <- PSE[, 1]
  
  #qDp
  qDp <- qDp(matched$phy, matched$comm, q = 2)
  CPD$qDP <- qDp
  
  new <- CPD[, c("STAND_KEY", "ntaxa", "pd.obs", "RaoD", "mpd", "mntd", "PSV", "PSE", "PSR", "qDP")]
  
  fwrite(new, paste0(out_path, "/master_phylogenetic-diversity.txt"), sep = "\t")
  
}