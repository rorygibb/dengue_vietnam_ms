
# ================== Themes and materials for plotting and visualisation ===================

# vietnam extent
ext = extent(c(100, 110, 7.5, 24))

# mapping theme
maptheme = theme_classic() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(), 
        axis.ticks = element_blank(),
        plot.title = element_text(hjust=0.5, size=14),
        legend.title = element_text(size=10), 
        strip.background = element_blank(),
        strip.text = element_text(size=16))
