
# Low-DF explorations



```{r cutcode, include=FALSE, echo=FALSE, results='asis'}
fig6cap <- "Figure 6: Boxplot of estimated standard error inflation for matched studies for the fully interacted CRVE estimators. Standard error inflation is measured as the ratio of the estimated standard error over the true super population standard error. Values above 1 indicate overestimation, values below 1 indicate underestimation. The faint horizontal lines go out to the full range, the bars are the 1% to 99%, and the boxes are the 10% to 90%."

crve<-sruns.clean %>% 
  filter(method %in% c("LRcw-FIcw-crve","LR-FIpw-crve")) 

inf_plot_lowdf_crve <- create_custom_boxplot2(
  crve %>% filter(J.bar=="S"), 
  "inflate.s", 
  "", 
  "Ratio of Estimated SE over True Super SE",
  vline_at = 1,
  log_scale = TRUE
) + theme_minimal() +
  theme(
    text = element_text(family = "serif", size = 10),
    axis.text.x = element_text(size = 8, family = "serif"),
    axis.text.y = element_text(size = 8, hjust = 0.5, family = "serif"),
    plot.title = element_text(hjust = 0.5, size = 12, family = "serif"),
    plot.caption = element_text(hjust = 0, size = 9, family = "serif", color = "gray40"),
    legend.text = element_text(size = 9, family = "serif"),
    legend.title = element_text(size = 10, family = "serif"),
    panel.grid = element_blank()
  )

ggdraw() +
  draw_plot(inf_plot_lowdf_crve, 0, 0.1, 1, 0.9) +
  draw_label(
    str_wrap(fig6cap, width = 120),
    x = 0.05, y = 0.05,
    hjust = 0, vjust = 0.5,
    size = 9, color = "gray40", fontfamily = "serif"
  )
```
