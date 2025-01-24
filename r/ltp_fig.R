ltp_fig <- function(
    data_rec,
    data_pre,
    data_old,
    fit_ltp,
    lab) {
  pars_s <- fit_ltp$summary(c("thetaInf_s"), "median") %>%
    separate(variable, c("par", "s", "sitenum"), convert = TRUE) %>%
    select(-s) %>%
    pivot_wider(names_from = par, values_from = median)
  pars_p <- fit_ltp$summary(c("dist_p", "lambda_p"), "mean") %>%
    separate(variable, c("par", "p", "plotnum"), convert = TRUE) %>%
    select(-p) %>%
    pivot_wider(names_from = par, values_from = mean)
  projs <- data_rec %>%
    select(site, plot, sitenum, plotnum) %>%
    left_join(pars_s, by = join_by(sitenum)) %>%
    left_join(pars_p, by = join_by(plotnum)) %>%
    mutate(rel_year = list(1:60)) %>%
    unnest(rel_year) %>%
    mutate(stem = ltp_pred(rel_year, dist, thetaInf, lambda))
  g <- projs %>%
    ggplot(aes(rel_year, stem)) +
    geom_boxplot(data = data_old %>% mutate(rel_year = 61),
                 width = 4, fill = NA) +
    geom_boxplot(data = data_pre %>% mutate(rel_year = -1),
                 width = 4, fill = NA) +
    geom_point(aes(group = paste(site, plot), col = plot),
               data = data_rec, alpha = .5) +
    geom_line(aes(group = paste(site, plot), col = plot)) +
    facet_wrap(~site, scales = "free") +
    theme_bw() +
    xlab("") +
    theme(legend.position = "bottom") +
    ylab(lab) +
    scale_color_discrete(guide = "none")
  return(g)
}

ltp_pred <- function(time, dist, thetainf, lambda) {
  theta0 <- thetainf * dist
  ltp <- 1 - exp(-lambda * time)
  return(theta0 + (thetainf - theta0) * ltp)
}
