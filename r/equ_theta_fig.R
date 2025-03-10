equ_theta_fig <- function(data) {
  t <- data$data_equ %>%
    group_by(site, sitenum) %>%
    summarise(
      median = median(value),
      sd = sd(value),
      q5 = quantile(value, .0),
      q95 = quantile(value, 1)
    ) %>%
    mutate(type = "data")
  data$fit_equ$summary("mu") %>%
    separate(variable, c("mu", "sitenum"), convert = TRUE) %>%
    left_join(data$ind_equ) %>%
    mutate(type = "posterior") %>%
    bind_rows(t) %>%
    ggplot(aes(median, paste(site, type), col = type)) +
    geom_segment(aes(x = median - sd, xend = median + sd)) +
    geom_segment(aes(x = q5, xend = q95), size = .3) +
    geom_point() +
    theme_bw() +
    xlab(data$lab) +
    ylab(expression(theta[equilibirum])) +
    scale_color_discrete(guide = "none")
}
