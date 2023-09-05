###
###
#' 
#' Authors: Pablo Capilla-Lasheras
#' 
#' Last update 2023-09-05
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### libraries and functions #####
##
##
pacman::p_load(ggplot2, cowplot, magick)

#####

##
##
##### Plots for figure #####
##
##

plot3b <- readRDS(file = './plots/Figure 3b.RDS')
plot3c <- readRDS(file = './plots/Figure 3c.RDS')
plot3d <- readRDS(file = './plots/Figure 3d.RDS')
plot3e <- readRDS(file = './plots/Figure 3e.RDS')
plot3f <- readRDS(file = './plots/Figure 3f.RDS')

#####

##
##
##### Panel Figure 1 #####
##
##
Figure3 <- ggdraw() +
  draw_image(image = "./plots/v6.0 - Encounter tag - wider.jpg", x = 0.00, y = 0.55, width = 0.25, height = 0.45, scale = 0.9) +
  draw_plot(plot3b, x = 0.25, y = 0.55, width = 0.25, height = 0.45) +
  draw_plot(plot3c, x = 0.50, y = 0.55, width = 0.25, height = 0.45) +
  draw_plot(plot3d, x = 0.75, y = 0.55, width = 0.25, height = 0.45) +
  draw_plot(plot3e, x = 0.00, y = 0.00, width = 0.50, height = 0.50) +
  draw_plot(plot3f, x = 0.50, y = 0.00, width = 0.50, height = 0.50) +
  draw_plot_label(label = c("A", "B", "C", "D", "E", 'F'), size = 15,
                  x = c(0.00, 0.25, 0.50, 0.75, 0.00, 0.50), y = c(1, 1, 1, 1, 0.50, 0.50)) 


ggsave(filename = "./plots/Figure 3.png", 
       plot = Figure3, 
       units = "mm",
       device = "png", 
       width = 250,
       height = 200)

