ogf_data_fig <- function(data) {
  data$data_ogf %>%
    ggplot(aes(y, fill = plot)) +
    geom_histogram() +
    theme_bw() +
    facet_wrap(~site, scales = "free_y") +
    scale_fill_discrete(guide = "none") +
    xlab(data$lab) +
    ylab("")
}
