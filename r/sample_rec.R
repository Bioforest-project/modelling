sample_rec <- function(var, sit = NULL) {
  if(var %in% c("ba", "agb", "nstem"))
    model <- cmdstan_model("models/recovery_positive.stan")
  else
    model <- cmdstan_model("models/recovery.stan")
  chains_path <- file.path("chains", paste0("rec_", var))
  unlink(chains_path, recursive = TRUE)
  dir.create(chains_path)
  data <- read_tsv("data/derived_data/data_rec.tsv", col_types = cols()) %>%
    filter(variable == var) %>%
    arrange(sitenum, plotnum, rel_year)
  if(!is.null(sit))
    data <- data %>%
      filter(site == sit) %>%
      mutate(sitenum = 1) %>% 
      mutate(plotnum = as.numeric(as.factor(plotnum)))
  ind <- data %>%
    select(site, plot, sitenum, plotnum) %>%
    unique()
  mu_priors <- read_tsv(file.path("chains", paste0("equ_", var), "mu.tsv"))
  if(!is.null(sit))
    mu_priors <- filter(mu_priors, site == sit)
  fit <- model$sample(
    data = list(
      n = nrow(data),
      s = max(data$sitenum),
      p = max(data$plotnum),
      y = data$y,
      t = data$rel_year - 3,
      site = data$sitenum,
      plot = data$plotnum,
      site_plot = ind$sitenum,
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
