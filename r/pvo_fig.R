pvo_fig <- function(
    data_rec,
    data_pre
    ){
  rmse_rec <- sqrt(mean((data_rec$preds - data_rec$stem)^2, na.rm = TRUE))
  g_rec <- data_rec %>%
    ggplot(aes(preds, stem)) +
    geom_point(alpha = 0.25) +
    geom_abline(col = "red", linetype = "dashed") +
    ggtitle("Recovery", paste("RMSE =", round(rmse_rec, 2))) +
    theme_bw() +
    xlab("Predicted") +
    ylab("Observed")
  data_pre$preds <- fit_ltp$summary(c("mu_pre"), median)$median
  rmse_pre <- sqrt(mean((data_pre$preds - data_pre$stem)^2, na.rm = TRUE))
  g_pre <- data_pre %>%
    ggplot(aes(preds, stem, col = site)) +
    geom_point(alpha = 0.25) +
    geom_abline(col = "red", linetype = "dashed") +
    ggtitle("Prelogging", paste("RMSE =", round(rmse_pre, 2))) +
    theme_bw() +
    scale_color_discrete("") +
    theme(legend.position = "bottom") +
    xlab("Predicted") +
    ylab("Observed")
  return(g_pre + g_rec)
}
