proj_fig <- function(
    data_rec,
    data_pre,
    data_old,
    fit,
    lab,
    sit = NA) {
  projs <- fit$summary("y_pred", "median") %>%
    separate(variable,
      c("y", "pred", "rel_year", "plotnum_rec"),
      convert = TRUE
    ) %>%
    select(-y, -pred) %>%
    rename(y = median) %>%
    mutate(rel_year = rel_year + 3) %>%
    left_join(data_rec %>%
                select(site, plot, plotnum_rec) %>%
                unique())
  if (!is.na(sit)) {
    projs <- filter(projs, site == sit)
    data_pre <- filter(data_pre, site == sit)
    data_rec <- filter(data_rec, site == sit)
    data_old <- filter(data_old, site == sit)
  }
  g <- projs %>%
    ggplot(aes(rel_year, y)) +
    geom_boxplot(
      data = data_old %>% mutate(rel_year = 61),
      width = 4, fill = NA
    ) +
    geom_boxplot(
      data = data_pre %>% mutate(rel_year = -1),
      width = 4, fill = NA
    ) +
    geom_point(aes(group = paste(site, plot), col = plot),
      data = data_rec, alpha = .5
    ) +
    geom_line(aes(group = paste(site, plot), col = plot)) +
    facet_wrap(~site, scales = "free") +
    theme_bw() +
    xlab("") +
    theme(legend.position = "bottom") +
    ylab(lab)
  if (is.na(sit)) {
    g <- g +
      scale_color_discrete(guide = "none")
  } else {
    g <- g +
      scale_color_discrete("Plot") +
      theme(legend.position = "right")
  }
  return(g)
}
