ogf_theta_fig <- function(data) {
  data$fit_ogf$summary("mu") %>%
    separate(variable, c("mu", "sitenum"), convert = TRUE) %>%
    left_join(data$ind_ogf) %>%
    na.omit() %>% 
    ggplot(aes(median, site)) +
    geom_segment(aes(x = median - sd, xend = median + sd)) +
    geom_segment(aes(x = q5, xend = q95), size = .1) +
    geom_point() +
    theme_bw() +
    xlab(data$lab) +
    ylab(expression(theta[OGF]))
}
