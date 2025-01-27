prep_mdata <- function(
    var,
    time_max = 60) {
  data_rec <- read_tsv("data/derived_data/data_rec.tsv") %>%
    filter(variable == var)
  data_pre <- read_tsv("data/derived_data/data_pre.tsv") %>%
    filter(variable == var)
  data_old <- read_tsv("data/derived_data/data_old.tsv") %>%
    filter(variable == var)
  ind_rec <- data_rec %>%
    select(site, plot, sitenum, plotnum_rec) %>%
    unique() %>%
    arrange(plotnum_rec)
  bounds <- read_tsv("data/derived_data/bounds.tsv") %>%
    filter(attribute == var)
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
    mu_thetaInf_bounds = c(
      bounds$mu_thetaInf_min,
      bounds$mu_thetaInf_max
    ),
    thetaInf_bounds = c(
      bounds$thetaInf_min,
      bounds$thetaInf_max
    ),
    lambda_bounds = c(
      bounds$lambda_min,
      bounds$lambda_max
    ),
    dist_bounds = c(
      bounds$dist_min,
      bounds$dist_max
    ),
    delta_bounds = c(
      bounds$delta_min,
      bounds$delta_max
    ),
    tau_bounds = c(
      bounds$tau_min,
      bounds$tau_max
    )
  )
  return(mdata)
}
