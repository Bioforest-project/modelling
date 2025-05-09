rec_theta_fig <- function(data) {
  bind_rows(
    data$fit_rec$summary("theta_s") %>%
      separate(variable, c("mu", "s", "sitenum"), convert = TRUE) %>%
      left_join(data$ind_rec) %>%
      select(site, mu, median, sd, q5, q95) %>%
      mutate(type = "posterior"),
    data$fit_equ$summary("mu") %>%
      separate(variable, c("mu", "sitenum"), convert = TRUE) %>%
      left_join(data$ind_equ) %>%
      select(site, mu, median, sd, q5, q95) %>%
      mutate(type = "prior")
  ) %>%
    na.omit() %>%
    ggplot(aes(median, site, col = type)) +
    geom_segment(aes(x = median - sd, xend = median + sd)) +
    geom_segment(aes(x = q5, xend = q95), size = .1) +
    geom_point() +
    theme_bw() +
    xlab(lab) +
    ylab(expression(theta[recovery])) +
    scale_color_discrete("") +
    theme(legend.position = "bottom")
}
