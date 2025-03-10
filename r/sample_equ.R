sample_equ <- function(var) {
  model <- cmdstan_model("models/equilibirum.stan")
  chains_path <- file.path("chains", paste0("equ_", var))
  unlink(chains_path, recursive = TRUE)
  dir.create(chains_path)
  data <- read_tsv("data/derived_data/data_equ.tsv") %>%
    filter(variable == var) %>%
    arrange(sitenum, plotnum, rel_year)
  ind <- data %>%
    select(sitenum, site) %>%
    unique()
  fit <- model$sample(
    data = list(
      n_obs = nrow(data),
      n_site = max(data$sitenum),
      n_plot = max(data$plotnum),
      y = data$y,
      site = data$sitenum,
      plot = data$plotnum
    ),
    chains = 4,
    parallel_chains = 4,
    save_warmup = FALSE,
    output_dir = chains_path,
    refresh = 100
  )
  fit$summary("mu") %>%
    separate(variable, c("theta", "sitenum"), convert = TRUE) %>%
    left_join(ind) %>%
    write_tsv(file.path(chains_path, "mu.tsv"))
  print(paste(var, "done!"))
}
