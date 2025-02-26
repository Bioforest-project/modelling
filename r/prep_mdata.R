prep_mdata <- function(
    var,
    time_max = 60) {
  data_rec <- read_tsv("data/derived_data/data_rec.tsv", col_types = cols()) %>%
    filter(variable == var)
  data_pre <- read_tsv("data/derived_data/data_pre.tsv", col_types = cols()) %>%
    filter(variable == var)
  data_old <- read_tsv("data/derived_data/data_old.tsv", col_types = cols()) %>%
    filter(variable == var)
  ind_rec <- data_rec %>%
    select(site, plot, sitenum, plotnum_rec) %>%
    unique() %>%
    arrange(plotnum_rec)
  ind_old <- data_old %>%
    select(site, plot, sitenum, plotnum_old) %>%
    unique() %>%
    arrange(plotnum_old)
  bounds <- read_tsv("data/derived_data/bounds.tsv", col_types = cols()) %>%
    filter(attribute == var)
  mdata <- list(
    n_rec = nrow(data_rec),
    n_old = nrow(data_old),
    n_pre = nrow(data_pre),
    n_site = max(data_rec$sitenum),
    n_plot_rec = max(data_rec$plotnum_rec),
    n_plot_old = max(data_old$plotnum_old),
    y_rec = data_rec$y,
    y_old = data_old$y,
    y_pre = data_pre$y,
    time = data_rec$rel_year - 3,
    site_rec = data_rec$sitenum,
    site_old = data_old$sitenum,
    site_pre = data_pre$sitenum,
    plot_rec = data_rec$plotnum_rec,
    plot_old = data_old$plotnum_old,
    plot_pre = data_pre$plotnum_pre,
    site_plot_rec = ind_rec$sitenum,
    site_plot_old = ind_old$sitenum,
    time_max = time_max,
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
