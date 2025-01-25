pvo_fig <- function(
    data_rec,
    data_pre,
    fit) {
  data_rec$preds <- fit$summary(c("mu_rec"), median)$median
  data_pre$preds <- fit$summary(c("mu_pre"), median)$median
  rmse_rec <- sqrt(mean((data_rec$preds - data_rec$y)^2, na.rm = TRUE))
  rmse_pre <- sqrt(mean((data_pre$preds - data_pre$y)^2, na.rm = TRUE))
  g_rec <- data_rec %>%
    ggplot(aes(preds, y)) +
    geom_point(alpha = 0.25) +
    geom_abline(col = "red", linetype = "dashed") +
    ggtitle("Recovery", paste("RMSE =", round(rmse_rec, 2))) +
    theme_bw() +
    xlab("Predicted") +
    ylab("Observed")
  g_pre <- data_pre %>%
    ggplot(aes(preds, y, col = site)) +
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
