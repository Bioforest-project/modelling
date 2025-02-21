sample_model <- function(
    var,
    model) {
  model_path <- file.path("models", paste0(model, ".stan"))
  stan_model <- cmdstan_model(model_path)
  chains_path <- file.path("chains", var, model)
  unlink(chains_path, recursive = TRUE)
  dir.create(chains_path)
  mdata <- prep_mdata(var)
  fit <- stan_model$sample(
    data = mdata,
    chains = 4,
    parallel_chains = 4,
    save_warmup = FALSE,
    output_dir = chains_path,
    refresh = 100
  )
  print(paste(model, var, "done!"))
}
