mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

optimal_figure <- function(models, rAxes, rMetric, rEstimate) {
  
  sub <- models[Metric == rMetric & Estimate == rEstimate & Axes == rAxes, ]
  onesigma <- mode(sub$onasigma)
  randomization <- mode(sub$randomization)
  min_value <- mode(sub$abs_min)
  
  sub_mean <- apply(sub[, 9:ncol(sub)], 2, mean)
  sub_sd <- apply(sub[, 9:ncol(sub)], 2, sd)
  
  p <- ggplot() +
    geom_ribbon(aes(x = 1:length(sub_mean),
                    ymin = sub_mean - sub_sd,
                    ymax = sub_mean + sub_sd),
                colour = "grey", alpha = 0.5) +
    geom_line(aes(x = 1:length(sub_mean), y = sub_mean)) +
    geom_vline(xintercept = c(min_value, onesigma, randomization), colour = c("green", "red", "blue")) +
    scale_x_continuous(limits = c(1, length(sub_mean)), expand = c(0, 0), breaks = 1:length(sub_mean)) +
    xlab("Components") + ylab(paste0(rMetric, " - ", rEstimate))
  
  return(p)
  
}
