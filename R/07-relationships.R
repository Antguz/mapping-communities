library(data.table)

masters_path <- "E:mapping_species/results"
masters_path <- "/media/antonio/antonio_ssd/mapping_species/results"

structural <- fread(paste0(masters_path, "/", "master_FIA-stand.txt"))
functional <- fread(paste0(masters_path, "/", "master_FIA-functional.txt"))
climate <- fread(paste0(masters_path, "/", "master_FIA-climate.txt"))  

compare <- merge(structural, climate, by = c("PLOTID", "ID", "LAT", "LON"), all.x = TRUE, all.y = FALSE)
compare <- compare[flag_sp != TRUE]
compare <- compare[flag_ba != TRUE]
compare <- compare[flag_tpa != TRUE]
compare <- compare[flag_tb != TRUE]

compare <- merge(compare, functional, by = "PLOTID", all.x = TRUE, all.y = TRUE)
compare <- compare[npixels >= 307]

cor <- compare[, c("date", "LAT", "LON", "nsp", "TBA", "TB", "FRic_TB", "FDis_TB", "Rao_TB", 
                 "BIO1", "BIO2", "BIO3", "BIO4", "BIO5", "BIO6", "BIO7", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", 
                 "BIO15", "BIO16", "BIO17", "BIO18", "BIO19", "elevation", "SR", "WVP",
                 "mean_KND", "mean_CCI", "mean_CRE", "mean_NDW", "FRic_all", "FDiv_all", "FEve_all", "FDis_all", "Pmd_all", "Rao_all")]

unique_date <- unique(cor$date)

for(i in 1:length(unique_date)) {
  
  cor_date <- subset(cor, date == unique_date[i])
  cor_date <- cor_date[, c(2:41)]
  cor_date <- na.exclude(cor_date)
  corr <- round(cor(cor_date), 6)
  
  name_file <- paste0("figures/", unique_date[i], ".jpeg")

  plot <- ggcorrplot(corr, outline.col = "white", title = as.character(unique_date[i]))
  
  ggsave(
    filename = name_file,
    plot = plot,
    width = 30,
    height = 30,
    units = "cm",
    dpi = 300,
    bg = "white")
  
  
}

plot(log10(cor_date$BIO4), log10(cor_date$FRic_all))


plot(cor_date$BIO5, log10(cor_date$FRic_all)) #Max temp
plot(cor_date$BIO6, log10(cor_date$FRic_all)) #Min temp
plot(log10(cor_date$BIO7), log10(cor_date$FRic_all)) #Range of temperature
plot(log10(cor_date$elevation), log10(cor_date$FRic_all))
plot(cor_date$LAT, log10(cor_date$FRic_all))
plot(cor_date$nsp, log10(cor_date$FRic_all))

plot(cor_date$BIO5, log10(cor_date$FDiv_all)) #Max temp
plot(cor_date$BIO6, log10(cor_date$FDiv_all)) #Min temp
plot(log10(cor_date$BIO7), log10(cor_date$FDiv_all)) #Range of temperature
plot(log10(cor_date$elevation), log10(cor_date$FDiv_all))
plot(cor_date$LAT, log10(cor_date$FDiv_all))

plot(log10(cor_date$BIO5), log10(cor_date$FEve_all)) #Max temp
plot(cor_date$BIO6, log10(cor_date$FEve_all)) #Min temp
plot(log10(cor_date$BIO7), cor_date$FEve_all) #Range of temperature
plot(log10(cor_date$elevation), log10(cor_date$FEve_all))
plot(cor_date$LAT, log10(cor_date$FEve_all))


plot(log10(cor_date$BIO15), log10(cor_date$FRic_all))
plot(log10(cor_date$elevation), log10(cor_date$FRic_all))
plot(log10(cor_date$SR), log10(cor_date$FRic_all))
plot(log10(cor_date$WVP), log10(cor_date$FRic_all))
plot(cor_date$nsp, cor_date$FRic_all)


plot(log10(cor_date$BIO4), log10(cor_date$FDiv_all))
plot(log10(cor_date$BIO5), log10(cor_date$FDiv_all))
plot(log10(cor_date$BIO15), log10(cor_date$FDiv_all))
plot(log10(cor_date$elevation), log10(cor_date$FDiv_all))
plot(log10(cor_date$SR), log10(cor_date$FDiv_all))
plot(log10(cor_date$WVP), log10(cor_date$FDiv_all))

plot(log10(cor_date$BIO4), log10(cor_date$FEve_all))
plot(log10(cor_date$BIO5), log10(cor_date$FEve_all))
plot(log10(cor_date$BIO15), log10(cor_date$FEve_all))
plot(log10(cor_date$elevation), log10(cor_date$FEve_all))
plot(log10(cor_date$SR), log10(cor_date$FEve_all))
plot(log10(cor_date$WVP), log10(cor_date$FEve_all))

plot(log10(cor_date$BIO4), log10(cor_date$FDis_all))
plot(log10(cor_date$BIO5), log10(cor_date$FDis_all))
plot(log10(cor_date$BIO15), log10(cor_date$FDis_all))
plot(log10(cor_date$elevation), log10(cor_date$FDis_all))
plot(log10(cor_date$SR), log10(cor_date$FDis_all))
plot(log10(cor_date$WVP), log10(cor_date$FDis_all))

plot(cor_date$BIO7, cor_date$FDiv_all)
plot(cor_date$BIO15, cor_date$FEve_all)
plot(log10(cor_date$BIO15), log10(cor_date$FDis_all))
abline(lm(log10(cor_date$BIO15) ~ log10(cor_date$BIO7)))
plot(log10(cor_date$BIO15), log10(cor_date$Pmd_all))
plot(log10(cor_date$BIO15), log10(cor_date$Rao_all))


cor <- frame[, c("MN_SIZE", "MN_DENS", "UTYPE", "USIZE", "UDENSE", "STAND_AGE", "TOTAL_ACRE", "TMBR_STAT", "TOPO", "PHYS", "ORIGIN", "BASAL_AREA",
                 "TOTAL_CORD", "TOTAL_MBF", "N_SPP", "TOTAL_VOL", "FRic_VOL", "FDiv_VOL", "FEve_VOL", "FDis_VOL", "Rao_VOL", "FRic_DBH", "FDiv_DBH",  
                 "FEve_DBH", "FDis_DBH", "Rao_DBH", "ntaxa", "pd.obs", "RaoD", "mpd", "mntd", "PSV", "PSE", "PSR", "qDP", "npixels", "mean_KND", 
                 "mean_CCI", "mean_CRE", "mean_NDW", "sd_KND", "sd_CCI", "sd_CRE", "sd_NDW", "FRic_CCI", "FRic_CRE", "FRic_NDW", "FRic_all", 
                 "FDiv_CCI", "FDiv_CRE", "FDiv_NDW", "FDiv_all", "FEve_CCI", "FEve_CRE", "FEve_NDW", "FEve_all", "FDis_CCI", "FDis_CRE", 
                 "FDis_NDW", "FDis_all", "Pmd_CCI", "Pmd_CRE", "Pmd_NDW", "Pmd_all", "Rao_all")]

#BIO4
#BIO5

library(CCA)
library(CCP)

cor_date <- cor_date[mean_KND >= tanh(0.7^2)]
cor_date <- cor_date[mean_KND >= tanh(0.6^2)]

cor_date <- na.exclude(cor_date)

X <- data.table(Temp_max = log10(cor_date$BIO5),
                Temp_seasonality = log10(cor_date$BIO4),
                Temp_range = log10(cor_date$BIO7),
                Precipitation_seasonality = log10(cor_date$BIO15),
                Elevation = log10(cor_date$elevation),
                Solar_radiation = log10(cor_date$SR),
                WVP = log10(cor_date$WVP))
 
Y <- data.table(FRic = log10(cor_date$FRic_all),
                FDiv = log10(cor_date$FDiv_all),
                FEve = log10(cor_date$FEve_all))

cc_results <- cancor(X, Y)
cc_results$cor

CC1_X <- as.matrix(x) %*% cc_results$xcoef[, 1:7]
CC1_Y <- as.matrix(y) %*% cc_results$ycoef[, 1:5]

plot(CC1_X[, 1], CC1_Y[, 2]) 

spearman <- maxCorGrid(X, Y, method = "kendall", consistent = TRUE)
perm <- permTest(x, y, R = 100, method = "kendall", nCores = 26, seed = 2016)



plot(log10(cor_date$BIO2), log10(cor_date$FRic_all))
plot(log10(cor_date$BIO5), log10(cor_date$FRic_all))
plot(log10(cor_date$BIO15), log10(cor_date$FRic_all))
plot(log10(cor_date$elevation), log10(cor_date$FRic_all))
plot(log10(cor_date$SR), log10(cor_date$FRic_all))
plot(log10(cor_date$WVP), log10(cor_date$FRic_all))


cormat <- matcor(X, Y)
img.matcor(cormat, type = 2)
can_cor1 <- cc(X, Y)


if (require(corrplot)) {
  M <- cor(cbind(X,Y))
  corrplot(M, method="ellipse", order="hclust", addrect=2, addCoef.col="black")
}

pca <- prcomp(X, center = TRUE, scale = TRUE)

plot(pca$x[,1], log10(cor_date$FRic_all))
lm(log10(cor_date$FRic_all) ~ pca$x[,1])

plot(pca$x[,1], log10(cor_date$FDiv_all))
plot(pca$x[,1], log10(cor_date$FEve_all))
plot(pca$x[,1], log10(cor_date$FDis_all))
lm(log10(cor_date$FDis_all) ~ pca$x[,1])

plot(pca$x[,1], log10(cor_date$Rao_all))
