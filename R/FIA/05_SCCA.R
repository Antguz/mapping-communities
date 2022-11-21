#'------------------------------------------------------------------------------
#' @title Sparse Canonical Correlation Analysis
#'------------------------------------------------------------------------------

#' @description Its match the spectral information with the reshape stands
#' @return .txt file with the stack of results

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(sRDA)
library(sgdm)
library(PMA)

options(scipen = 100)

#'------------------------------------------------------------------------------
#' Arguments
#' @param fia
#' @param out_path
#' @param threads

root_path <- "E:/DESIS/04_match_spectra-stand"
root_path <- "/media/antonio/Extreme SSD/DESIS/04_match_spectra-stand"

stand <- fread(paste0(root_path, "/species.csv"))
spectra <- fread(paste0(root_path, "/spectra.csv"), header = TRUE)
spectra <- spectra[, c(1:3, 13:238)]
out_path <- paste0(root_path, "/04_match_spectra-stand")

#'------------------------------------------------------------------------------
#' Analysis

library(PMA)

set.seed(22)

perm.out <- CCA.permute(x= p, 
                        z= r,
                        typex = "standard",
                        typez = "standard",
                        nperms = 10,
                        penaltyx = seq(.2,.9,len=100),
                        penaltyz = seq(.2,.9,len=100))

out <- CCA(x = r,
           z = p,
           typex="standard",
           typez="standard",
           K = ncol(r),
           penaltyx = perm.out$bestpenaltyx,
           penaltyz = perm.out$bestpenaltyz,
           v = perm.out$v.init)


## We run CCA using all gene exp. data, but CGH data on chrom 1 only.
print(perm.out)
plot(perm.out)
out <- CCA(x=rna,
           z=dna[,chrom==1], 
           typex="standard", 
           typez="ordered",
           penaltyx=perm.out$bestpenaltyx,
           v=perm.out$v.init, 
           penaltyz=perm.out$bestpenaltyz,
           xnames=substr(genedesc,1,20),
           znames=paste("Pos", sep="", nuc[chrom==1]))
  # Save time by inputting  lambda and v
  print(out) # could do print(out,verbose=TRUE)
  print(genechr[out$u!=0]) # Cool! The genes associated w/ gain or loss
  ## on chrom 1 are located on chrom 1!!
  par(mfrow=c(1,1))
  PlotCGH(out$v, nuc=nuc[chrom==1], chrom=chrom[chrom==1],
          main="Regions of gain/loss on Chrom 1 assoc'd with gene expression")
  detach(breastdata)

## End(Not run)

