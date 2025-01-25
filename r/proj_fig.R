proj_fig <- function(
    data_rec,
    data_pre,
    data_old,
    fit,
    lab) {
  projs <- fit$summary("y_pred", "median") %>% 
    separate(variable, 
             c("y", "pred", "rel_year", "plotnum"), convert = TRUE) %>% 
    select(-y, -pred) %>% 
    rename(y = median) %>% 
    left_join(data_rec %>% 
                select(site, plot, plotnum) %>% 
                unique())
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
    ylab(lab) +
    scale_color_discrete(guide = "none")
  return(g)
}
