sample_ogf <- function(var) {
  ogf <- cmdstan_model("models/ogf.stan")
  chains_path <- file.path("chains", paste0("ogf_", var))
  unlink(chains_path, recursive = TRUE)
  dir.create(chains_path)
  data_ogf <- read_tsv("data/derived_data/data_ogf.tsv") %>%
    filter(variable == var) %>%
    arrange(sitenum, plotnum, rel_year)
  ind_ogf <- data_ogf %>%
    select(sitenum, site) %>%
    unique()
  fit_ogf <- ogf$sample(
    data = list(
      n_obs = nrow(data_ogf),
      n_site = max(data_ogf$sitenum),
      n_plot = max(data_ogf$plotnum),
      y = data_ogf$y,
      site = data_ogf$sitenum,
      plot = data_ogf$plotnum
    ),
    chains = 4,
    parallel_chains = 4,
    save_warmup = FALSE,
    output_dir = chains_path,
    refresh = 100
  )
  fit_ogf$summary("mu") %>%
    separate(variable, c("theta", "sitenum"), convert = TRUE) %>%
    left_join(ind_ogf) %>%
    write_tsv(file.path(chains_path, "mu.tsv"))
  print(paste(var, "done!"))
}
