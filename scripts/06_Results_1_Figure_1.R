###
###
#' 
#' Authors: Pablo Capilla-Lasheras
#' 
#' Last update 2023-09-04
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
pacman::p_load(ggplot2, cowplot)
               
#####

##
##
##### Plots for figure #####
##
##
plot1a <- readRDS(file = './plots/Figure 1a.RDS')
plot1b <- readRDS(file = './plots/Figure 1b.RDS')
plot1c <- readRDS(file = './plots/Figure 1c.RDS')
plot1d <- readRDS(file = './plots/Figure 1d.RDS')
plot1e <- readRDS(file = './plots/Figure 1e.RDS')

#####

##
##
##### Panel Figure 1 #####
##
##
Figure1 <- ggdraw() +
  draw_plot(plot1a, x = 0.00, y = 0.50, width = 0.50, height = 0.50) +
  draw_plot(plot1b, x = 0.50, y = 0.50, width = 0.50, height = 0.50) +
  draw_plot(plot1c, x = 0.00, y = 0.05, width = 0.33, height = 0.45) +
  draw_plot(plot1d, x = 0.33, y = 0.05, width = 0.33, height = 0.45) +
  draw_plot(plot1e, x = 0.66, y = 0.05, width = 0.33, height = 0.45) +
  draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 15,
                  x = c(0, 0.5, 0, 0.33, 0.66), y = c(1, 1, 0.50, 0.50, 0.50)) +
  draw_label("Subordinate sex", colour = "black", size = 15, angle = 0, x = 0.5, y = 0.03)


ggsave(filename = "./plots/Figure 1.png", 
       plot = Figure1, 
       units = "mm",
       device = "png", 
       width = 170,
       height = 225)

