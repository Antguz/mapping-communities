#'------------------------------------------------------------------------------
#' @title Quality assesement of spectra
#'------------------------------------------------------------------------------

#' @description It helps to do a quality assesement of the spectra stracted
#' @return .txt file with the stack of results

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(spectrolab)

#'------------------------------------------------------------------------------
#' Arguments

root_path <- "E:/Projects/DESIS_mapping/02_compile-files"
spectra <- fread(paste0(root_path, "/02-spectra_SG.txt"))

#Product filters
spectra <- spectra[Shadow != 1]
spectra <- spectra[`Haze (land)` != 1]
spectra <- spectra[`Cloud (land)` != 1]
spectra <- spectra[Water != 1]

#Filters by NDVI
NDVI <- (spectra$`801` - spectra$`629.6`) / (spectra$`801` + spectra$`629.6`)
spectra <- spectra[NDVI >= 0.5,]

NDVI <- (spectra$`801` - spectra$`629.6`) / (spectra$`801` + spectra$`629.6`)
NDVI <- round(NDVI, 2)
CCI <- (spectra$`532.2` - spectra$`629.6`) / (spectra$`532.2` + spectra$`629.6`)
CCI <- round(CCI, 2)
spectra$sample <- 1:nrow(spectra)


#sample <- spectra[, .SD, .SDcols = 14:208]
#sample <- as_spectra(as.matrix(sample), name_idx = NULL, meta_idxs = NULL)

Wavelength <- as.numeric(colnames(spectra[, .SD, .SDcols = 14:208]))

for(i in 1:nrow(sample)) {
  
  Reflectance <- as.numeric(pixel[i, .SD, .SDcols = 14:208])
  
  plot <- ggplot() +
    geom_point(aes(x = Wavelength, y = Reflectance), size = 2) +
    coord_cartesian(xlim = c(420, 930), ylim = c(0, 10000)) +
    theme_bw() +
    labs(title = paste0("Sample:  ", pixel$sample[i]),
         subtitle = paste0("NDVI: ", NDVI[i],"  CCI: ", CCI[i]),
         caption = paste0("Min: ", round(min(Reflectance), 2),
                          "  Mean: ", round(mean(Reflectance), 2),
                          "  Max: ", round(max(Reflectance), 2)))
  
  ggsave(filename = paste0(root_path, "/quality/sample_", pixel$sample[i], ".png"), 
         plot,
         device = "png")
  
}

