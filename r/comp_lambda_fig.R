comp_lambda_fig <- function(
    data_rec,
    fit_stp,
    fit_ltp) {
  ind_rec <- data_rec %>%
    select(site, plot, sitenum, plotnum_rec) %>%
    unique() %>%
    arrange(plotnum_rec)
  fit_ltp$summary("lambda_p", "mean") %>%
    rename(ltp = mean) %>%
    left_join(
      fit_stp$summary("lambda_p", "mean") %>%
        rename(stp = mean)
    ) %>%
    separate(variable, c("lambda", "p", "plotnum_rec"), convert = TRUE) %>%
    left_join(ind_rec) %>%
    ggplot(aes(ltp, stp, col = site)) +
    geom_abline(linetype = "dashed") +
    geom_point() +
    theme_bw() +
    scale_y_log10() +
    scale_x_log10() +
    ggpubr::stat_cor(aes(group = NA)) +
    scale_color_discrete("") +
    xlab("LTP") +
    ylab("STP")
}
