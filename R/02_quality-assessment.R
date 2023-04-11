#'------------------------------------------------------------------------------
#' @title Quality assesement of spectra
#'------------------------------------------------------------------------------

#' @description It helps to do a quality assesement of the spectra stracted
#' @return .txt file with the stack of results

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(ggplot2)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Arguments

root_path <- "/media/antonio/Work/Projects/DESIS_mapping/02_compile-files"
#root_path <- "E:/Projects/DESIS_mapping/02_compile-files"
spectra <- fread(paste0(root_path, "/02-spectra_SG.txt"))

#Product filters
spectra <- spectra[Shadow != 1]
spectra <- spectra[`Haze (land)` != 1]
spectra <- spectra[`Cloud (land)` != 1]
spectra <- spectra[Water != 1]

#Filters by NDVI
NDVI <- (spectra$`801` - spectra$`629.6`) / (spectra$`801` + spectra$`629.6`)
hist(NDVI)
spectra <- spectra[NDVI >= 0.6 & NDVI <= 1.0,]
NIR <- spectra$'804'
hist(NIR)
spectra <- spectra[NIR >= 1500]
CCI <- (spectra$`532.2` - spectra$`629.6`) / (spectra$`532.2` + spectra$`629.6`)
hist(CCI)
spectra <- spectra[CCI >= 0.0 & CCI <= 1.0]


NDVI <- (spectra$`801` - spectra$`629.6`) / (spectra$`801` + spectra$`629.6`)
NDVI <- round(NDVI, 2)
CCI <- (spectra$`532.2` - spectra$`629.6`) / (spectra$`532.2` + spectra$`629.6`)
CCI <- round(CCI, 2)
NIR <- round(spectra$'804', 2)
spectra$sample <- 1:nrow(spectra)
Wavelength <- as.numeric(colnames(spectra[, .SD, .SDcols = 14:208]))

master <- spectra[, c("folder", "scene", "PLT_CN")]
fwrite(master, paste0(root_path, "/master_spectra.csv"))

#nrow(sample)

for(i in 1:nrow(spectra)) {
  
  pixel <- spectra[i, ]
  
  Reflectance <- as.numeric(pixel[1, .SD, .SDcols = 14:208])/10000
  
  if(CCI[i] <= 0) {
    col <- "red"
  } else {
    col <- "green"
  }
  
  plot <- ggplot() +
    geom_point(aes(x = Wavelength, y = Reflectance), size = 2, shape = 21, fill = col) +
    coord_cartesian(xlim = c(420, 930), ylim = c(0, 0.75)) +
    theme_bw() +
    labs(title = paste0("Sample:  ", pixel$sample[1]),
         subtitle = paste0("NDVI: ", NDVI[i],"  CCI: ", CCI[i], "  NIR: ", round(pixel$'804'[1]/10000, 2)),
         caption = paste0("Min: ", round(min(Reflectance), 2),
                          "  Mean: ", round(mean(Reflectance), 2),
                          "  Max: ", round(max(Reflectance), 2)))
  
  ggsave(filename = paste0(root_path, "/quality/sample_", pixel$sample[1], ".png"), 
         plot,
         device = "png")
  
  rm("plot")
  
}

