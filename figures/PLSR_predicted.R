data <- evaluation$Predicted
predicted <- data[, .SD, .SDcols = c(3:ncol(predicted))]
predicted_mean <- apply(predicted, 1, mean)
predicted_sd <- apply(predicted, 1, sd)

frame <- data[, .SD, .SDcols = c(1:2)]
frame$predicted_mean <- predicted_mean
frame$predicted_sd <- predicted_sd

frame[MDS == 1, label := "NMDS1"]
frame[MDS == 2, label := "NMDS2"]
frame[MDS == 3, label := "NMDS3"]

reg <- evaluation$Performance
summary <- reg[, lapply(.SD, mean, na.rm=TRUE), by= "MDS"]



plot <- ggplot(frame, aes(x = predicted_mean, y = observed)) +
  geom_abline(intercept = 0, slope = 1, colour = "grey", linetype = "solid", linewidth = 0.5) +
  geom_errorbar(aes(xmin= predicted_mean-predicted_sd, xmax= predicted_mean+predicted_sd), colour = "#134f5cff") + 
  geom_point(shape = 21, colour = "white", size = 2, fill = "#134f5cff") +
  geom_smooth(method='lm', se = FALSE, colour = "black", linewidth = 0.5) +
  theme_bw(base_size = 12) + 
  xlab("Predicted") + ylab("Observed") +
  facet_wrap(~ label, scales = "free") +
  scale_x_continuous(limits = c(-1.25, 1.25), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1.25, 1.25), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jpeg("figures/PLSR_taxonomic.jpg", width = 210, height = 90, units = "mm", res = 300)
plot
dev.off()

ci <- function(x, conf = 0.95, na.rm = FALSE){
  xbar <- mean(x, na.rm = na.rm)
  s <- sd(x, na.rm = na.rm)
  p <- c((1 - conf)/2, 1 - (1 - conf)/2)
  qq <- qnorm(p, mean = xbar, sd = s)
  setNames(qq, c("lower", "upper"))
}
stats <- function(x, na.rm = FALSE){
  CI <- ci(x, na.rm = na.rm)
  c(
    Mean = mean(x, na.rm = na.rm),
    SD = sd(x, na.rm = na.rm),
    Lower = CI[1],
    Upper = CI[2]
  )
}


VIP <- model$VIP
melt_VIP <- melt(VIP, id.vars = c("MDS", "repeated"), 
                 variable.name = "wavelength",
                 value.name = "value")
melt_VIP <- melt_VIP[, c(1, 3, 4)]
VIP_summary <- summarySE(data = melt_VIP, measurevar = "value", groupvars = c("MDS",  "wavelength"),
          na.rm = FALSE, conf.interval = 0.95, .drop = TRUE)
VIP_summary <- as.data.table(VIP_summary)
VIP_summary$wavelength <- as.numeric(as.character(VIP_summary$wavelength))

VIP_summary[MDS == 1, label := "NMDS1"]
VIP_summary[MDS == 2, label := "NMDS2"]
VIP_summary[MDS == 3, label := "NMDS3"]

plot <- ggplot(VIP_summary, aes(x = wavelength, y = value)) +
  geom_line(colour = "#134f5cff", linewidth = 1) +
  geom_ribbon(aes(ymin = value-sd, ymax= value+sd), fill = "#134f5cff", alpha = 0.25) + 
  theme_bw(base_size = 12) + 
  xlab("Wavelength (nm)") + ylab("Variable of Importance") +
  facet_wrap(~ label, scales = "free") +
  scale_x_continuous(limits = c(460, 900), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jpeg("figures/VIP_taxonomic.jpg", width = 210, height = 90, units = "mm", res = 300)
plot
dev.off()


coeff <- model$Coefficients
melt_coef <- melt(coeff, id.vars = c("MDS", "repeated"), 
                 variable.name = "wavelength",
                 value.name = "value")
melt_coef <- subset(melt_coef, wavelength != "intercept")
melt_coef <- melt_coef[, c(1, 3, 4)]
coef_summary <- summarySE(data = melt_coef, measurevar = "value", groupvars = c("MDS",  "wavelength"),
                         na.rm = FALSE, conf.interval = 0.95, .drop = TRUE)
coef_summary <- as.data.table(coef_summary)
coef_summary$wavelength <- as.numeric(as.character(coef_summary$wavelength))

coef_summary[MDS == 1, label := "NMDS1"]
coef_summary[MDS == 2, label := "NMDS2"]
coef_summary[MDS == 3, label := "NMDS3"]

plot <- ggplot(coef_summary, aes(x = wavelength, y = value)) +
  geom_line(colour = "#134f5cff", linewidth = 1) +
  geom_ribbon(aes(ymin = value-sd, ymax= value+sd), fill = "#134f5cff", alpha = 0.25) + 
  theme_bw(base_size = 12) + 
  xlab("Wavelength (nm)") + ylab("PLS coefficients") +
  facet_wrap(~ label, scales = "free") +
  scale_x_continuous(limits = c(460, 900), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.0013, 0.0013), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jpeg("figures/coef_taxonomic.jpg", width = 210, height = 90, units = "mm", res = 300)
plot
dev.off()

