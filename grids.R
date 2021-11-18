library(cowplot)
library(magick)

lenupset <- image_read_svg('anlen_setsize_2.svg',
                           width = 3000)

p2 <- ggdraw() +
      draw_image(lenupset, scale = 1.0)

p2

plot_grid(A,
          p2,
          heat,
          ncol = 2,
          nrow = 2,
          rel_widths = c(1,1,20),
          labels = "auto")


library(gridExtra)

svg(file="grid.svg",
    width = 7.08, #just under 180 mm
    height = 7.08,
    pointsize = 12,
    onefile = FALSE)

grid.arrange(A,
             p2,
             heat,
             layout_matrix = rbind(c(1,2),
                                   c(3,3),
                                   c(3,3)))

dev.off()
