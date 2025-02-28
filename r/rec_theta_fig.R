rec_theta_fig <- function(data) {
  bind_rows(
    data$fit_stp$summary("theta_s") %>%
      separate(variable, c("mu", "s", "sitenum"), convert = TRUE) %>%
      left_join(data$ind_rec) %>%
      select(site, mu, median, sd, q5, q95) %>%
      mutate(type = "posterior"),
    data$fit_ogf$summary("mu") %>%
      separate(variable, c("mu", "sitenum"), convert = TRUE) %>%
      left_join(data$ind_ogf) %>%
      select(site, mu, median, sd, q5, q95) %>%
      mutate(type = "prior")
  ) %>%
    ggplot(aes(median, site, col = type)) +
    geom_segment(aes(x = median - sd, xend = median + sd)) +
    geom_segment(aes(x = q5, xend = q95), size = .1) +
    geom_point() +
    theme_bw() +
    xlab(lab) +
    ylab(expression(theta)) +
    scale_color_discrete("") +
    theme(legend.position = "bottom")
}
