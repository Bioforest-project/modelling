proj_fig <- function(data) {
  data$fit_stp$summary("y_pred", "median") %>%
    separate(variable,
      c("y", "pred", "rel_year", "plotnum"),
      convert = TRUE
    ) %>%
    select(-y, -pred) %>%
    rename(y = median) %>%
    mutate(rel_year = rel_year + 3) %>%
    left_join(data$ind_rec) %>%
    ggplot(aes(rel_year, y)) +
    geom_boxplot(
      data = data$data_ogf %>%
        filter(treatment == "control") %>%
        mutate(rel_year = 61),
      width = 4, fill = NA
    ) +
    geom_boxplot(
      data = data$data_ogf %>%
        filter(treatment == "logged") %>%
        mutate(rel_year = -1),
      width = 4, fill = NA
    ) +
    geom_point(aes(group = paste(site, plot), col = plot),
      data = data$data_rec, alpha = .5
    ) +
    geom_line(aes(group = paste(site, plot), col = plot)) +
    facet_wrap(~site, scales = "free") +
    theme_bw() +
    xlab("") +
    theme(legend.position = "bottom") +
    ylab(data$lab) +
    scale_color_discrete(guide = "none")
}
