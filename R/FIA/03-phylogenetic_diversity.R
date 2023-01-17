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
source("R/FIA/comdist_rcpp.R")

#'------------------------------------------------------------------------------
#' Arguments
#' @param fia
#' @param SPCD
#' @param out_path
#' @param threads

root_path <- "F:/DESIS"
root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
fia <- fread(paste0(root_path, "/00_FIA/trees_WI_MN_20221102.csv"))
SPCD <- fread(paste0(root_path, "/00_FIA/SPCD.csv"), na.strings="NA")
out_path <- paste0(root_path, "/03_diversity")
threads <- 30

#'------------------------------------------------------------------------------
#' @example 

#'------------------------------------------------------------------------------
#' Function

phylogenetic_diversity <- function(master, threads) {
  
  ### Prepare plots ------------------------------------------------------------
  #Clean FIA
  fia <- subset(fia, STATECD == 27 | STATECD == 55)
  fia <- fia[order(PLT_CN, SUBP)]
  
  #Merge species code
  fia <- merge(fia, 
               SPCD[, c("SPCD", "PLANTS-Code", "Taxon2")], 
               by = "SPCD", 
               all.x = TRUE, 
               all.y = TRUE)
  
  #Get Basal Area
  fia$BA <- pi*(fia$DIA/2)^2
  
  #Basal area multiplied by TPA
  fia$BA_TPA <- fia$BA*fia$TPA_UNADJ
  
  #Remove dead trees
  fia <- fia[STATUSCD == 1]
  
  #Get plots with complete information
  complete_plots <- fia[, sum(is.na(BA_TPA)), by = "PLT_CN"]
  complete_plots <- complete_plots[V1 == 0]
  
  #Select plots with complete information
  fia <- merge(complete_plots[, 1], 
               fia, 
               by = "PLT_CN", 
               all.x = TRUE, 
               all.y = FALSE)
  
  #Summarize
  fia_summary <- fia[, .(ntrees = .N,
                         HT = mean(HT, na.rm = TRUE), 
                         BA_TPA = sum(BA, na.rm = TRUE)), 
                     by = c("PLT_CN", "PLANTS-Code", "Taxon2")]
  fia_summary <- fia_summary[order(PLT_CN)]
  
  ### Prepare phylo  -----------------------------------------------------------
  #Makes phylogenetic hypotheses
  sppPhylo <- data.table(species = SPCD$Taxon2, 
                         genus = SPCD$Genus, 
                         family = SPCD$Family)
  sppPhylo <- na.exclude(sppPhylo)
  hypotheses <- phylo.maker(sppPhylo, scenarios = "S3")
  
  #Resolve Multichotomies
  phylo <- multi2di(hypotheses$scenario.3)
  is.binary.phylo(phylo) #Test for Binary Tree
  is.ultrametric(phylo) #Test if a Tree is Ultrametric
  
  ### Community matrix ---------------------------------------------------------
  
  #dcast
  master_matrix <- sample2matrix(fia_summary[, c("PLT_CN", "BA_TPA", "Taxon2")])
  
  # Remove plots without species
  nsp <- apply(master_matrix, 1, function(x) sum(x > 0))
  master_matrix <- master_matrix[nsp != 0, ]
  
  #Log transform 
  master_matrix <- log10(master_matrix + 1)
  
  ### Get phylogenetic indices -------------------------------------------------
  
  matched <- picante::match.phylo.comm(phy = phylo, comm = master_matrix)
  saveRDS(matched$comm, paste0(out_path, "/master_community.rds"))
 
  #Calculates inter-community mean pairwise distance
  # this function uses the Rcpp modification
  ic_MDP <- comdist_rcpp(matched$comm, 
                         cophenetic(matched$phy), 
                         abundance.weighted = TRUE, 
                         threads = 28)
  saveRDS(ic_MDP, paste0(out_path, "/MDP_distance.rds")) #Export beta distance
  
  #Phylogenetic diversity in a community - Faith's PD
  CPD <- ses.pd(matched$comm, matched$phy, runs = 99)
  CPD <- CPD[, c(1, 2, 6, 7)]
  
  #Add PLT_CN and arrange
  CPD$PLT_CN <- rownames(CPD)
  CPD <- CPD[, c(5, 1, 2, 3, 4)]
  CPD <- as.data.table(CPD)
  
  #Mean pairwise distance separating taxa in a community - MPD
  SESMPD <- ses.mpd(matched$comm, 
                    cophenetic(matched$phy), 
                    null.model = "taxa.labels", 
                    runs = 99)
  CPD$mpd <- SESMPD[, c(2)]
  CPD$mpd.obs.z <- SESMPD[, c(6)]
  CPD$mpd.obs.p <- SESMPD[, c(7)]
  
  #Mean nearest taxon distance for taxa in a community - MNTD
  SESMNTD <- ses.mntd(matched$comm, 
                      cophenetic(matched$phy), 
                      null.model = "taxa.labels", 
                      runs = 99)
  CPD$mntd <- SESMNTD[, c(2)]
  CPD$mntd.obs.z <- SESMNTD[, c(6)]
  CPD$mntd.obs.p <- SESMNTD[, c(7)]
  
  #Phylogenetic Rao’s quadratic entropy - RaoD
  CPD$RaoD <- raoD(matched$comm, force.ultrametric(matched$phy))$Dkk
  
  #Phylogenetic species variability - PSV
  PSV <- psv(matched$comm, 
             matched$phy, 
             compute.var = TRUE)
  CPD$PSV <- PSV[, 1]
  CPD$PSV_vars <- PSV[, 3]
  
  #Phylogenetic species richness - PSR
  PSR <- psr(matched$comm, 
             matched$phy, 
             compute.var = TRUE)
  CPD$PSR <- PSR[, 1]
  CPD$PSR_vars <- PSR[, 3]
  
  #Phylogenetic species evenness - PSE
  PSE <- pse(matched$comm, 
             matched$phy)
  CPD$PSE <- PSE[, 1]
  
  #Phylogenetic species clustering - PSC
  PSC <- psc(matched$comm, matched$phy)
  CPD$PSC <- PSC[, 1]
  
  #qDp
  rqDp <- qDp(matched$phy, 
              matched$comm, 
              q = 2)
  CPD$qDP <- rqDp
  
  fwrite(CPD, paste0(out_path, "/master_phylogenetic-diversity.txt"), sep = "\t")
  
}