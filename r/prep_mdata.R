prep_mdata <- function(
    data_rec,
    data_pre,
    data_old,
    dist_bounds,
    delta_bounds,
    mu_thetaInf_bounds,
    thetaInf_s_bounds,
    time_max = 60
  ) {
  ind_rec <- data_rec %>%
    select(site, plot, sitenum, plotnum_rec) %>%
    unique() %>%
    arrange(plotnum_rec)
  mdata <- list(
    n_rec = nrow(data_rec),
    n_old = nrow(data_old),
    n_pre = nrow(data_pre),
    n_site = max(data_rec$sitenum),
    n_plot_rec = max(data_rec$plotnum_rec),
    y_rec = data_rec$y,
    y_old = data_old$y,
    y_pre = data_pre$y,
    time = data_rec$rel_year - 3,
    site_rec = data_rec$sitenum,
    site_old = data_old$sitenum,
    site_pre = data_pre$sitenum,
    plot_rec = data_rec$plotnum_rec,
    site_plot_rec = ind_rec$sitenum,
    time_max = time_max,
    dist_bounds = dist_bounds,
    delta_bounds = delta_bounds,
    mu_thetaInf_bounds = mu_thetaInf_bounds,
    thetaInf_s_bounds = thetaInf_s_bounds
  )
  return(mdata)
}
