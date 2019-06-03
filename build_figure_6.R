# Script to build Figure 6 in "Alleviating hypoxia through induced downwelling"


#Reactive transport modeling using 'ReacTran' package to model plume fate
#Numerical alternative to analytical plume modeling

# Written by David Koweek on 14 May 2019

#----Initialize_workspace----

#Load necessary packages
library(tidyverse)
library(ReacTran)
library(viridis)
library(cowplot)
library(ggpubr)

#----Set_model_parameters----

x_range <- 20 #meters
y_range <- 20 

dx <- 0.05 #meters
dy <- 0.05 #meters

n_x <- x_range / dx #number of grid cells in x-direction
n_y <- y_range / dy #number of grid cells in y-direction

#Diffusion
D <- 5e-3 #m^2/s

#Define advective velocities below in each of the cases
# v <- 0 #m/s (advective velocity)

O2_initial <- 0 #intial concentration
O2_input <- 2 #mg/L
R <- 0.1 * O2_input # production/consumption rate  

# O2 input pipe coordinates
pipe_radius <- 0.5 #meters

#----Set_model_grid----
grid <- matrix(nrow = n_y,
               ncol = n_x)

y_grid <- 
  ((row(grid) * dy) - (dy/2))  #multiply each grid cell by the dy increment (and then center to get the midpoint value)
  
x_grid <-
  (col(grid) * dx) - (dx/2)  #multiply each grid cell by the dx increment (and then center to get the midpoint value)

radial_grid <- 
  sqrt((x_grid - (x_range / 2))^2 + (y_grid - (y_range / 2))^2) #center at x = 0, y = 0

pipe_location <- 
  which(radial_grid <= pipe_radius,
        arr.ind = TRUE)

outside_locations <- 
  which(radial_grid > pipe_radius,
        arr.ind = TRUE)

ICs <- grid
ICs[pipe_location] <- O2_input
ICs[outside_locations] <- O2_initial

#----Set_timescales
pumping_duration <- 3 #hrs
times <- seq(0, 3600 * pumping_duration, by = 60) #seconds

#----Define_reactive_transport_function----
 Diff2D <- function (t, y, parms, N) { 
   
   CONC <- matrix(nrow = n_y, ncol = n_x, y) 
   
    # Transport 
     Tran <-
       tran.2D(C = CONC, 
               D.x = D, 
               D.y = D, 
               dx = dx, 
               dy = dy,
               v.x = u,
               v.y = v)
 
   # transport + reaction (respiratory losses)
       dCONC <- Tran$dC - (R / 3600) #seconds-> hours

  # Hold pipe outflow concentration constant     
    dCONC[pipe_location] <- 
      O2_input - CONC[pipe_location]
       
       
  return (list(dCONC))  
     
}


#----Run_model----

u <- 0 #m/s (advective velocity)
v <- 0 #m/s (advective velocity)  
  
plume_out_mixing <- 
    ode.2D(func = Diff2D, 
           y = as.vector(ICs), 
           times = times, 
           N = n, 
           parms = NULL, 
           lrw = 2e7, 
           dimens = c(n_y, n_x))

velocity <- 0.01 #m/s (advective velocity)
flow_angle <- 0 

u <- -velocity * sin((pi/180) * flow_angle) #moving in -x direction
v <- velocity * cos((pi/180) * flow_angle)
  
plume_out_mixing_transport <- 
    ode.2D(func = Diff2D, 
           y = as.vector(ICs), 
           times = times, 
           N = n, 
           parms = NULL, 
           lrw = 2e7, 
           dimens = c(n_y, n_x))


#----Wrangle_model_results_for_plotting----

#For each model test case
model_cases <- c("plume_out_mixing", "plume_out_mixing_transport")   
plume_fate_df <- list()
    
  for (i in 1:length(model_cases)) {
  
    #Convert into a data frame
    df <- 
      matrix(nrow = n_y, 
             ncol = n_x, 
             subset(eval(parse(text = model_cases[i])), time == max(times))
            ) %>% 
      as.data.frame()
    
    #Name rows and columns with x and y values
    rownames(df) <- c(1:n_y) * dy
    colnames(df) <- c(1:n_x) * dx
    
    #Gather wide data into tall data and re-center at x = 0, y = 0  
    plume_fate_df[[i]] <-
      df %>% 
      rownames_to_column() %>% 
      dplyr::rename(x = rowname) %>% 
      gather(key = y, value = plume, -x) %>% 
      mutate(x = as.numeric(x) - 10,
             y = as.numeric(y) - 10,
             plume = case_when((plume < 0) ~ 0,
                               TRUE ~ plume))

  }

#Free up memory by clearing deSolve model objects (which retain all time steps)
rm(plume_out_mixing, plume_out_mixing_transport)

#----Plot_plume_model_results----

plume_plots_theme <- 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10))

colour_limits <- 
  c(0, O2_input)

circle_df <- function(r = 0.5, n_points = 360){
  
  points <- c(0:n_points)
  x <- r * cos(points * (2 * pi / 360))
  y <- r * sin(points * (2 * pi / 360))
  
  return(data.frame(x = x,
                    y = y))
  
}

plume_plot_mixing <-
  plume_fate_df[[1]] %>%
  ggplot(aes(x = x,
             y = y,
             z = plume)) +
  geom_raster(aes(fill = plume)) +
  geom_contour(colour = "white",
               breaks = 0.5 * O2_input) +
  geom_contour(colour = "white",
               breaks = 0.25 * O2_input,
               linetype = "dashed") +
  scale_fill_viridis(name = expression(O[2]~"`"~(mg~L^{-1})),
                     limits = colour_limits) +
  scale_x_continuous(name = expression(symbol('\254')~~Distance~(m)~~symbol('\256')),
                     breaks = seq(-x_range/2, x_range/2, by = 2),
                     limits = c(-8,8)) +
  scale_y_continuous(name = expression(symbol('\254')~~Distance~(m)~~symbol('\256')),
                     breaks = seq(-y_range/2, y_range/2, by = 2),
                     limits = c(-8,8)) +
  #Add in outline of pipe
  geom_path(inherit.aes = FALSE,
            data = circle_df(r = 0.5),
            aes(x = x,
                y = y),
            colour = "black",
            alpha = 0.25) +
  coord_fixed() +
  plume_plots_theme +
  theme_bw()

plume_plot_mixing_transport <-
  plume_fate_df[[2]] %>%
  ggplot(aes(x = x,
             y = y,
             z = plume)) +
  geom_raster(aes(fill = plume)) +
  geom_contour(colour = "white",
               breaks = 0.5 * O2_input) +
  geom_contour(colour = "white",
               breaks = 0.25 * O2_input,
               linetype = "dashed") +
  scale_fill_viridis(name = expression(O[2]~"`"~(mg~L^{-1})),
                     limits = colour_limits) +
  scale_x_continuous(name = expression(symbol('\254')~~Distance~(m)~~symbol('\256')),
                     breaks = seq(-x_range/2, x_range/2, by = 2),
                     limits = c(-8,8)) +
  scale_y_continuous(name = expression(symbol('\254')~~Distance~(m)~~symbol('\256')),
                     breaks = seq(-y_range/2, y_range/2, by = 2),
                     limits = c(-8,8)) +
  #Add in outline of pipe
  geom_path(inherit.aes = FALSE,
            data = circle_df(r = 0.5),
            aes(x = x,
                y = y),
            colour = "black",
            alpha = 0.25) +
  coord_fixed() +
  plume_plots_theme +
  theme_bw()

#----Create_figure----

plumes_figure <- 
  plot_grid(
    plume_plot_mixing + theme(legend.position = "none"),
    plume_plot_mixing_transport + theme(legend.position = "none"),
    ncol = 2,
    align = "v",
    labels = c("A) Idealized Mixing and Respiration",
               "B) Idealized Mixing, Transport, and Respiration"),
    label_size = 10,
    hjust = c(-0.25, -0.18)
  ) %>%
  annotate_figure(bottom = get_legend(
    plume_plot_mixing + theme(
      legend.title.align = 0.5,
      legend.direction = "horizontal",
      legend.key.width = unit(0.75, "in"),
      legend.justification = "center"
    )
  ))

#----Export_figure----

cowplot::ggsave(filename = "figures/figure_6.pdf",
                plot = plumes_figure,
                device = "pdf",
                height = 5,
                width = 8,
                units = "in")