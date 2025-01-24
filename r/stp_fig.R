stp_fig <- function(
    data_rec,
    data_pre,
    data_old,
    fit_stp,
    lab,
    max_ll = FALSE) {
  if (max_ll) {
    max_ll <- which(fit_stp$draws("lp__") == max(fit_stp$draws("lp__")),
                    arr.ind = TRUE)
    pars_s <- lapply(
      c("thetaInf_s", "tau_s"),
      function(par) {
        fit_stp$draws(par)[max_ll[1], max_ll[2], ] %>%
          as_tibble() %>%
          gather()
      }
    ) %>%
      bind_rows() %>%
      separate(key, c("chain", "par", "s", "sitenum"), convert = TRUE) %>%
      select(-s, -chain) %>%
      pivot_wider(names_from = par, values_from = value)
    pars_p <- lapply(
      c("dist_p", "lambda_p", "delta_p"),
      function(par) {
        fit_stp$draws(par)[max_ll[1], max_ll[2], ] %>%
          as_tibble() %>%
          gather()
      }
    ) %>%
      bind_rows() %>%
      separate(key, c("chain", "par", "p", "plotnum"), convert = TRUE) %>%
      select(-p, -chain) %>%
      pivot_wider(names_from = par, values_from = value)
  } else {
    pars_s <- fit_stp$summary(c("thetaInf_s", "tau_s"), "median") %>%
      separate(variable, c("par", "s", "sitenum"), convert = TRUE) %>%
      select(-s) %>%
      pivot_wider(names_from = par, values_from = median)
    pars_p <- fit_stp$summary(c("dist_p", "lambda_p", "delta_p"), "median") %>%
      separate(variable, c("par", "p", "plotnum"), convert = TRUE) %>%
      select(-p) %>%
      pivot_wider(names_from = par, values_from = median)
  }
  projs <- data_rec %>%
    select(site, plot, sitenum, plotnum) %>%
    left_join(pars_s, by = join_by(sitenum)) %>%
    left_join(pars_p, by = join_by(plotnum)) %>%
    mutate(rel_year = list(1:60)) %>%
    unnest(rel_year) %>%
    mutate(stem = stp_pred(rel_year, dist, thetaInf, lambda, tau, delta))
  g <- projs %>%
    ggplot(aes(rel_year, stem)) +
    geom_boxplot(
      data = data_old %>% mutate(rel_year = 61),
      width = 4, fill = NA
    ) +
    geom_boxplot(
      data = data_pre %>% mutate(rel_year = -1),
      width = 4, fill = NA
    ) +
    geom_point(aes(group = paste(site, plot), col = plot),
      data = data_rec, alpha = .5
    ) +
    geom_line(aes(group = paste(site, plot), col = plot)) +
    facet_wrap(~site, scales = "free") +
    theme_bw() +
    xlab("") +
    theme(legend.position = "bottom") +
    ylab(lab) +
    scale_color_discrete(guide = "none")
  return(g)
}

stp_pred <- function(time, dist, thetainf, lambda, tau, delta) {
  theta0 <- thetainf * dist
  ltp <- 1 - exp(-lambda * time)
  stp <- delta * (time / tau * exp(1 - time / tau))^2
  return(theta0 + (thetainf - theta0) * (ltp + stp))
}
