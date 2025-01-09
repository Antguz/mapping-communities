#'------------------------------------------------------------------------------
#' @title Community landmarks for Landmark Multi-Dimensional Scaling
#'------------------------------------------------------------------------------

#' Define communities to be considered as landmarks for Landmark Multi-Dimensional 
#' Scaling

#'------------------------------------------------------------------------------
#' @Source-code

source("R/auxiliary/lmds.R")

#'------------------------------------------------------------------------------
#' @Librarys

library(data.table)
library(bit64)

#'------------------------------------------------------------------------------
#' @Select-landmarks

# Define root path
root_path <- "/media/antonio/Work/MMDB"

# Read and get communities used for beta diversity estimations
comm <- readRDS(paste0(root_path, "/Diversity/01_beta/master_matched_comm_phylo_VOL.rds"))
PLT_CN <- as.integer64(rownames(comm))

# Read meta data
meta <- fread(paste0(root_path, "/FIA/00_FIA/plots_eastern_20230427.csv"))

# Match meta and actual communities used
matching_vector <- match(PLT_CN, meta$PLT_CN)
meta <- meta[matching_vector,]

# Define uniques
meta$STATECD_COUNTYCD <- paste0(meta$STATECD, "_", meta$COUNTYCD)
meta$ID_to_select <- 1:nrow(meta)
statecd_countycd <- unique(meta$STATECD_COUNTYCD)

# Define application
landmark <- function(X = 1:length(statecd_countycd),
                     statecd_countycd,
                     meta,
                     fraction = 0.005) {
  
  sub <- meta[STATECD_COUNTYCD == statecd_countycd[X],]
  samp <- sample(sub$ID_to_select, 1, replace = FALSE)
  return(samp)
  
}

# Apply application
comm_landmarks <- lapply(X = 1:length(statecd_countycd),
                         FUN = landmark,
                         statecd_countycd = statecd_countycd,
                         meta = meta,
                         fraction = 0.01)

comm_landmarks <- do.call(c, comm_landmarks)
comm_landmarks <- comm_landmarks[order(comm_landmarks)]

# Select and export landmarks
meta_landmarks <- meta[comm_landmarks, ]
nrow(meta_landmarks)/nrow(meta)
plot(meta_landmarks$LON, meta_landmarks$LAT)
comm_landmarks <- meta_landmarks$PLT_CN
saveRDS(comm_landmarks, paste0(root_path, "/Diversity/01_beta/comm_landmarks.rds"))
