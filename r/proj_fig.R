proj_fig <- function(data) {
  data$fit_rec$summary("y_pred", "median") %>%
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
      data = data$data_equ %>%
        filter(treatment == "control") %>%
        mutate(rel_year = 61),
      width = 4, fill = NA
    ) +
    geom_boxplot(
      data = data$data_equ %>%
        filter(treatment == "logged") %>%
        mutate(rel_year = -1),
      width = 4, fill = NA
    ) +
    geom_point(aes(group = paste(site, plotnum), col = as.character(plotnum)),
      data = data$data_rec, alpha = .5
    ) +
    geom_line(aes(group = paste(site, plotnum), col = as.character(plotnum))) +
    facet_wrap(~site, scales = "free") +
    theme_bw() +
    xlab("") +
    theme(legend.position = "bottom") +
    ylab(data$lab) +
    scale_color_discrete(guide = "none")
}
