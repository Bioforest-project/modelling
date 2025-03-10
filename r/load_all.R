load_all <- function(var, lab) {
  data_equ <- read_tsv("data/derived_data/data_equ.tsv") %>%
    filter(variable == var) %>%
    mutate(plot = as.character(plot))
  ind_equ <- data_equ %>%
    select(sitenum, site) %>%
    unique()
  path <- file.path("chains", paste0("equ_", var))
  fit_equ <- as_cmdstan_fit(list.files(path,
    full.names = TRUE, pattern = "csv"
  ))
  data_rec <- read_tsv("data/derived_data/data_rec.tsv", col_types = cols()) %>%
    filter(variable == var) %>%
    mutate(plot = as.character(plot))
  ind_rec <- data_rec %>%
    select(site, plot, sitenum, plotnum) %>%
    unique() %>%
    arrange(plotnum)
  path <- file.path("chains", paste0("rec_", var))
  fit_rec <- as_cmdstan_fit(list.files(path, full.names = TRUE))
  list(
    data_equ = data_equ,
    ind_equ = ind_equ,
    fit_equ = fit_equ,
    data_rec = data_rec,
    ind_rec = ind_rec,
    fit_rec = fit_rec,
    lab = lab
  )
}
