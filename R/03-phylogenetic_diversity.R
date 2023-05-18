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
library(geiger)
library(bit64)

#'------------------------------------------------------------------------------
#' Source code
#source("https://raw.githubusercontent.com/jesusNPL/BiodiversityScience/master/Spring2021/R-Functions/qDp.R")
#source("R/comdist_rcpp.R")
source("R/beta_rcpp.R")

#'------------------------------------------------------------------------------
#' Arguments
#' @param fia
#' @param SPCD
#' @param out_path
#' @param threads

root_path <- "F:/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
fia <- fread(paste0(root_path, "/00_FIA/trees_eastern_20230427.csv"))
SPCD <- fread(paste0(root_path, "/00_FIA/SPCD_new.csv"), na.strings="NA")
out_path <- paste0(root_path, "/03_diversity")
threads <- 30

#'------------------------------------------------------------------------------
#' @example 

#'------------------------------------------------------------------------------
#' Function

phylogenetic_diversity <- function(master, threads) {
  
  ### Prepare plots (Cleaning) -------------------------------------------------
  
  fia <- fia[order(PLT_CN, SUBP)]
  
  #Inventories from 2017 until now
  fia <- fia[MEASYEAR >= 2017]
  
  #Merge species code
  fia <- merge(fia, 
               SPCD[, c("SPCD", "PLANTS-Code-new", 
                        "Taxon2", "Genus", "Family", 
                        "Remove")], 
               by = "SPCD", 
               all.x = TRUE, 
               all.y = TRUE)
  
  #Remove dead trees
  fia <- fia[STATUSCD == 1]
  
  #Get forest metrics by TPA
  fia$BA <- pi*(fia$DIA/2)^2
  fia$BA_TPA <- fia$BA*fia$TPA_UNADJ
  fia$DRYBIO_AG_TPA <- fia$DRYBIO_AG*fia$TPA_UNADJ
  
  #Get plots with complete information (no missing BA)
  complete_plots <- fia[, sum(is.na(BA_TPA)), by = "PLT_CN"]
  complete_plots <- complete_plots[V1 == 0]
  
  fia <- merge(complete_plots[, 1], 
               fia, 
               by = "PLT_CN", 
               all.x = TRUE, 
               all.y = FALSE)
  
  #Select plots with complete information
  fia[is.na(Remove), Remove := 1]
  
  complete_plots <- fia[, sum(Remove), by = "PLT_CN"]
  complete_plots <- complete_plots[V1 == 0]
  
  fia <- merge(complete_plots[, 1], 
               fia, 
               by = "PLT_CN", 
               all.x = TRUE, 
               all.y = FALSE)
  
  #Summarize
  fia_summary <- fia[, .(ntrees = .N,
                         HT = mean(HT, na.rm = TRUE), 
                         BA_TPA = sum(BA, na.rm = TRUE), 
                         DRYBIO_AG_TPA = sum(DRYBIO_AG_TPA, na.rm = TRUE)),
                     by = c("PLT_CN", "PLANTS-Code-new", 
                            "Taxon2", "Genus", "Family")]
  fia_summary <- fia_summary[order(PLT_CN)]
  
  ### Prepare phylo  -----------------------------------------------------------
  #Summary of species and families
  species <- unique(fia_summary[, c("Taxon2", "Genus", "Family")])
  colnames(species) <- c("species", "genus", "family")
  species <- na.exclude(species)
  hypotheses <- phylo.maker(species, scenarios = "S3")
  
  #Resolve Multichotomies
  phylo <- multi2di(hypotheses$scenario.3)
  is.binary.phylo(phylo) #Test for Binary Tree
  is.ultrametric(phylo) #Test if a Tree is Ultrametric
  
  #Transform phylo
  tol = 1e-9
  phylo$edge.length[phylo$edge.length <= 0] <- tol
  is.ultrametric(phylo)
  
  #Scale phylo
  rtrns <- rescale(phylo, "nrate")
  phylo <- rtrns(time=c(0.25, 0.50, 0.75), rate=c(2.5, 5, 7.5))
  plot(phylo)
  is.ultrametric(phylo)
  
  ### Community matrix ---------------------------------------------------------
  
  #dcast
  master_matrix <- sample2matrix(fia_summary[, c("PLT_CN", "DRYBIO_AG_TPA", "Taxon2")])
  
  # Remove plots without species
  nsp <- apply(master_matrix, 1, function(x) sum(x > 0))
  master_matrix <- master_matrix[nsp != 0, ]
  
  #Standardization
  master_matrix <- decostand(master_matrix, method = "hellinger") 
  
  ### Get phylogenetic beta diversity ------------------------------------------
  
  matched <- picante::match.phylo.comm(phy = phylo, comm = master_matrix)
  saveRDS(matched$comm, paste0(out_path, "/master_matched_comm-DRYBIO.rds"))
  saveRDS(matched$phy, paste0(out_path, "/master_matched_phylo-tree.rds"))
 
  beta <- beta_rcpp(comm = as.matrix(matched$comm), 
                    tree = matched$phy, 
                    threads = 4)
  object.size(beta)
  saveRDS(beta, paste0(out_path, "/beta_DRYBIO.rds"))
  
}