load_all <- function(var, lab) {
  data_ogf <- read_tsv("data/derived_data/data_ogf.tsv") %>%
    filter(variable == var)
  ind_ogf <- data_ogf %>%
    select(sitenum, site) %>%
    unique()
  path <- file.path("chains", paste0("ogf_", var))
  fit_ogf <- as_cmdstan_fit(list.files(path,
    full.names = TRUE, pattern = "csv"
  ))
  data_rec <- read_tsv("data/derived_data/data_rec.tsv", col_types = cols()) %>%
    filter(variable == var)
  ind_rec <- data_rec %>%
    select(site, plot, sitenum, plotnum) %>%
    unique() %>%
    arrange(plotnum)
  path <- file.path("chains", paste0("stp_", var))
  fit_stp <- as_cmdstan_fit(list.files(path, full.names = TRUE))
  list(
    data_ogf = data_ogf,
    ind_ogf = ind_ogf,
    fit_ogf = fit_ogf,
    data_rec = data_rec,
    ind_rec = ind_rec,
    fit_stp = fit_stp,
    lab = lab
  )
}
