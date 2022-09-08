
masters_path <- "/media/antonio/antonio_ssd/mapping_species/results"

structural <- fread(paste0(masters_path, "/", "master_stand-diversity.txt"))
functional <- fread(paste0(masters_path, "/", "master_functional-diversity.txt"))
phylogenetic <- fread(paste0(masters_path, "/", "master_phylogenetic-diversity.txt"))  

compare <- merge(structural, phylogenetic, by = "STAND_KEY", all.x = TRUE, all.y = TRUE)
compare <- merge(compare, functional, by = "STAND_KEY", all.x = TRUE, all.y = TRUE)
frame <- subset(compare, Unique_UOM == "Cords")

cor <- frame[, c("date", "MN_SIZE", "MN_DENS", "STAND_AGE", "TOTAL_ACRE", "TOTAL_VOL", "FRic_VOL", "FDiv_VOL", "FEve_VOL", "FDis_VOL", "Rao_VOL", "ntaxa", 
                 "pd.obs", "RaoD", "mpd", "mntd", "PSV", "PSE", "PSR", "qDP", "FRic_all", "FDiv_all", "FEve_all", "FDis_all","Pmd_all", "Rao_all")]
cor <- na.exclude(cor)

unique_date <- unique(cor$date)

for(i in 1:length(unique_date)) {
  
  cor_date <- subset(cor, date == unique_date[i])
  corr <- round(cor(log10(cor_date[, 2:26]), 4))
  
  name_file <- paste0("figures/", unique_date[i], ".jpeg")

  plot <- ggcorrplot(corr, outline.col = "white", title = as.character(unique_date[i]))
  
  ggsave(
    filename =name_file,
    plot = plot,
    width = 18,
    height = 18,
    units = "cm",
    dpi = 300,
    bg = "white")
  
  
}






cor <- frame[, c("MN_SIZE", "MN_DENS", "UTYPE", "USIZE", "UDENSE", "STAND_AGE", "TOTAL_ACRE", "TMBR_STAT", "TOPO", "PHYS", "ORIGIN", "BASAL_AREA",
                 "TOTAL_CORD", "TOTAL_MBF", "N_SPP", "TOTAL_VOL", "FRic_VOL", "FDiv_VOL", "FEve_VOL", "FDis_VOL", "Rao_VOL", "FRic_DBH", "FDiv_DBH",  
                 "FEve_DBH", "FDis_DBH", "Rao_DBH", "ntaxa", "pd.obs", "RaoD", "mpd", "mntd", "PSV", "PSE", "PSR", "qDP", "npixels", "mean_KND", 
                 "mean_CCI", "mean_CRE", "mean_NDW", "sd_KND", "sd_CCI", "sd_CRE", "sd_NDW", "FRic_CCI", "FRic_CRE", "FRic_NDW", "FRic_all", 
                 "FDiv_CCI", "FDiv_CRE", "FDiv_NDW", "FDiv_all", "FEve_CCI", "FEve_CRE", "FEve_NDW", "FEve_all", "FDis_CCI", "FDis_CRE", 
                 "FDis_NDW", "FDis_all", "Pmd_CCI", "Pmd_CRE", "Pmd_NDW", "Pmd_all", "Rao_all")]

