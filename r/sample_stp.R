sample_stp <- function(var) {
  stp <- cmdstan_model("models/stp.stan")
  chains_path <- file.path("chains", paste0("stp_", var))
  unlink(chains_path, recursive = TRUE)
  dir.create(chains_path)
  data_rec <- read_tsv("data/derived_data/data_rec.tsv", col_types = cols()) %>%
    filter(variable == var) %>%
    arrange(sitenum, plotnum, rel_year)
  ind_rec <- data_rec %>%
    select(site, plot, sitenum, plotnum) %>%
    unique()
  mu_priors <- read_tsv(file.path("chains", paste0("ogf_", var), "mu.tsv"))
  fit_stp <- stp$sample(
    data = list(
      n = nrow(data_rec),
      s = max(data_rec$sitenum),
      p = max(data_rec$plotnum),
      y = data_rec$y,
      t = data_rec$rel_year - 3,
      site = data_rec$sitenum,
      plot = data_rec$plotnum,
      site_plot = ind_rec$sitenum,
      mu_theta_s = mu_priors$median,
      sigma_theta_s = mu_priors$sd
    ),
    chains = 4,
    parallel_chains = 4,
    save_warmup = FALSE,
    output_dir = chains_path,
    refresh = 100
  )
  print(paste(var, "done!"))
}
