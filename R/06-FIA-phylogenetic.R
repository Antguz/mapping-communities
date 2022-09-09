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

fia <- fread("data/FIA_SPCD.csv")
ID_species <- fread("data/SPCDindex.csv")
out_path <- "/media/antonio/antonio_ssd/mapping_species/results"
threads <- 30

#'------------------------------------------------------------------------------
#' @example 

#'------------------------------------------------------------------------------
#' Function

phylogenetic_diversity <- function(master, threads) {
  
  #Clean FIA
  fia <- subset(fia, STATECD == 27 | STATECD == 55)
  fia <- fia[order(PLOTID, SUBP)]
  
  #Summarize
  fia_summary <- fia[, .(BA = sum(BA, na.rm = TRUE), 
                         TPA = sum(TPA, na.rm = TRUE), 
                         SAPLING = sum(SAPLING, na.rm = TRUE), 
                         BIGTREE = sum(BIGTREE, na.rm = TRUE),
                         TREETALLY = sum(TREETALLY, na.rm = TRUE), 
                         TOTALBIOMASS = sum(TOTALBIOMASS, na.rm = TRUE)), 
                     by = c("PLOTID", "SPCD")]
  
  #Merge fia and species
  fia_summary <- merge(fia_summary, ID_species[, c(1:3)], by = "SPCD", all.x = TRUE, all.y = FALSE)
  
  #Check taxonomic
  ID_species <- unique(fia_summary$scientific)
  spp_check <- TPL(ID_species, infra = TRUE, corr = TRUE)
  
  taxonomy <- spp_check %>% 
    drop_na(New.Genus, New.Species) %>% 
    select(Taxon, Family, New.Genus, New.Species) #, Tax_res
  
  taxonomy$SP <- ID_species$id #Provide the ID
  colnames(taxonomy) <- c("species", "family", "genus", "sp", "SP")
  taxonomy <- taxonomy[, c("SP", "species", "genus", "sp", "family")]
  taxonomy$SPP <- paste0(taxonomy$genus, "_", taxonomy$sp)
  taxonomy$species[12] <- "Populus_X_deltoides"
  taxonomy$SPP[12] <- "Populus_X_deltoides"
  
  #Makes phylogenetic hypotheses
  sppPhylo <- taxonomy[, c(2, 3, 5)]
  hypotheses <- phylo.maker(sppPhylo, scenarios = "S3")
  
  phylo <- multi2di(hypotheses$scenario.3)
  is.binary.phylo(phylo)
  is.ultrametric(phylo)
  
  #Merge databases
  master <- merge(master, taxonomy, by = c("SP"), all.x = TRUE, all.y = TRUE)
  #Missing 7 20 29 90
  
  #Create matrix
  master_matrix <- sample2matrix(master[, c("STAND_KEY", "VOL", "SPP")])
  master_matrix <- round(master_matrix, 1)
  
  matched <- picante::match.phylo.comm(phy = phylo, comm = master_matrix)
  matched$phy
  matched$comm
  
  #Get Phylogenetic diversity
  PD <- pd(matched$comm, matched$phy, include.root = FALSE) # Faith's PD
  
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