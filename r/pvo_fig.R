pvo_fig <- function(data) {
  data$data_rec %>%
    arrange(sitenum, plotnum, rel_year) %>%
    mutate(preds = data$fit_stp$summary(c("mu"), median)$median) %>%
    ggplot(aes(preds, y)) +
    geom_point(alpha = 0.25) +
    geom_abline(col = "red", linetype = "dashed") +
    ggtitle(paste(
      "RMSE =",
      round(sqrt(mean((data$fit_stp$summary(c("mu"), median)$median -
        arrange(
          data$data_rec, sitenum, plotnum,
          rel_year
        )$y)^2, na.rm = TRUE)), 2)
    )) +
    theme_bw() +
    xlab("Predicted") +
    ylab("Observed")
}
