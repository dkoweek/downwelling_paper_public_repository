# Script to build Figure S6 in "Alleviating hypoxia through induced downwelling"

# Build lake overturn figure for oxygenation experiment study
# Written by David Koweek on 20 March 2019

#----Initialize_workspace----

#Load necessary packages
library(tidyverse)
library(lubridate)
library(scales)
library(viridis)
library(cowplot)

#Set system clock to avoid Mac OS X / lubridate problems with timestamps
Sys.setenv(TZ = "America/Los_Angeles")

#Load primary data set 
source("load_primary_data.R")

#Load meteorological station data
source("load_met_station_data.R")

#Load data from field notebooks
source("load_field_notes.R")

#----Define-plot-parameters----

datetime_limits <- #keep only complete days 
  ymd_hms(c("2018-09-15 00:00:00",
            "2018-10-03 23:59:00"),
          tz = "America/Los_Angeles")

#----Build-temperature-time-series-plot----

air_water_temperature_ts <-
  #Plot time series of water temperatures
  oxygenation_expt_sensor_data %>%
  filter(Mooring == "FF") %>%
  ggplot(aes(x = local_datetime,
             y = Temperature)) +
  geom_line(aes(colour = factor(Meters_above_bottom,
                                levels = c(2.5, 1.5, 0.5, 0.21)))) +
  # Add air temperature
  geom_line(data=met_data,
            aes(x = local_datetime,
                y = Temperature),
            colour="black",
            alpha = 0.5) +
  #Set scale details
  scale_colour_viridis(name = element_blank(),
                       labels = c("2.5 mab",
                                  "1.5 mab",
                                  "0.5 mab",
                                  "0.21 mab"),
                       discrete = TRUE,
                       option = "C",
                       direction = -1) +
  scale_y_continuous(name = expression(Temperature~~(degree~C))) +
  scale_x_datetime(name = element_blank(),
                   limits = datetime_limits,
                   date_labels = "%b %d",
                   date_breaks = "2 days") +
  #Set plot details
  theme_bw() +
  theme(axis.text.x = element_text(angle = -45))

#----Build-temperature-composite-daily-plot----

temperature_daily_cycle_plot <- 
  #Plot daily cycle of surface water temperature
  oxygenation_expt_sensor_data %>% 
  filter(Mooring == "FF",
         Meters_above_bottom == 2.5) %>%
  mutate(start_of_day = floor_date(local_datetime,
                                   unit = "day"), 
         time_of_day = difftime(local_datetime,
                                start_of_day,
                                units = c("hours")),
         time_of_day = as.numeric(time_of_day)) %>% 
  ggplot(aes(x = time_of_day,
             y = Temperature)) +
  geom_path(aes(group = date(local_datetime),
                colour = as.factor(Meters_above_bottom)),
            show.legend = FALSE) +
  #Add daily cycle of air temperature
  geom_path(data = met_data %>% 
                    mutate(start_of_day = floor_date(local_datetime,
                                                     unit = "day"), 
                           time_of_day = difftime(local_datetime,
                                                  start_of_day,
                                                  units = c("hours")),
                           time_of_day = as.numeric(time_of_day)),
            aes(x = time_of_day,
                y = Temperature,
                group = date(local_datetime)),
            colour = "black",
            alpha = 0.5) +
  #Scale details
  scale_colour_viridis(discrete = TRUE,
                       option = "C",
                       direction = -1) +
  scale_x_continuous(name = element_blank(),
                     limits = c(0, 24.01),
                     breaks = c(0,6,12,18,24),
                     labels = c("12 AM",
                                "6 AM",
                                "12 PM",
                                "6 PM",
                                "12 AM")) +
  scale_y_continuous(name = expression(Temperature~~(degree~C))) +
  theme_bw()


#----Build-O2-profile-time-series-plot----

FF_O2_ts <- 
  #Plot time series of O2
  oxygenation_expt_sensor_data %>%
  filter(Mooring == "FF") %>%
  ggplot(aes(x = local_datetime,
             y = O2)) +
  geom_line(aes(colour = factor(Meters_above_bottom,
                                levels = c(2.5, 1.5, 0.5, 0.21)))) +
  #Set scale details
  scale_colour_viridis(name = element_blank(),
                       labels = c("2.5 mab",
                                  "1.5 mab",
                                  "0.5 mab",
                                  "0.21 mab"),
                       discrete = TRUE,
                       option = "C",
                       direction = -1) +
  scale_y_continuous(name = expression(O[2]~~(mg~L^{-1}))) +
  scale_x_datetime(name = element_blank(),
                   limits = datetime_limits,
                   date_labels = "%b %d",
                   date_breaks = "2 days") +
  #Set plot details
  theme_bw() +
  theme(axis.text.x = element_text(angle = -45))

#----Calculate-depth-gradients-of-temperature-and-oxygen----

#First, calculate delta T (T_air - T_surface water)
delta_T <- 
  oxygenation_expt_sensor_data %>% 
  filter(Mooring == "FF",
         Meters_above_bottom %in% c(2.5, 0.21)) %>% 
  select(c(local_datetime, Temperature, Meters_above_bottom)) %>% 
  spread(Meters_above_bottom, Temperature) %>% 
  rename(surface_T = "2.5",
         bottom_T = "0.21") %>% 
  mutate(T_air_interp = approx(x = met_data[["local_datetime"]],
                               y = met_data[["Temperature"]],
                               xout = local_datetime)[["y"]]) %>% 
  mutate(delta_T_air_water = T_air_interp - surface_T,
         delta_T_water = surface_T - bottom_T)

#Second, calculate delta O2 (O2_surface - O2_bottom)
delta_O2 <- 
  oxygenation_expt_sensor_data %>% 
  filter(Mooring == "FF",
         Meters_above_bottom %in% c(0.21, 2.5)) %>% 
  select(c(local_datetime, O2, Meters_above_bottom)) %>% 
  spread(Meters_above_bottom, O2) %>% 
  rename(surface_O2 = "2.5",
         bottom_O2 = "0.21") %>% 
  mutate(delta_O2 = surface_O2 - bottom_O2)

#Join delta_T and delta_O2 together
delta_delta_df <- 
  delta_T %>% 
  left_join(delta_O2,
            by = "local_datetime") %>% 
  mutate(start_of_day = floor_date(local_datetime,
                                   unit = "day"), 
         time_of_day = difftime(local_datetime,
                                start_of_day,
                                units = c("hours")),
         time_of_day = as.numeric(time_of_day),
         hour_of_day = floor(time_of_day)) %>% 
  #Keep only complete days
  filter(local_datetime %within% interval(datetime_limits[1],
                                          datetime_limits[2],
                                          tzone = "America/Los_Angeles")) %>% 
  #Hourly bin-average the data to reduce high frequency noise
  group_by(start_of_day,
           hour_of_day) %>% 
  summarize(delta_T_bar = mean(delta_T_water, na.rm = TRUE),
            delta_O2_bar = mean(delta_O2, na.rm = TRUE)) %>% 
  ungroup()

#----Build-delta-temperature-delta-oxygen-plot----

#Now build plot
delta_T_delta_O2_plot <- 
  delta_delta_df %>% 
  ggplot(aes(x = delta_T_bar,
             y = delta_O2_bar)) +
  geom_point(aes(colour = hour_of_day)) +
  #Set scale details
  scale_colour_viridis(name = element_blank(),
                       option = "C",
                       breaks = c(0,6,12,18,23.9),
                       labels = c("12 AM",
                                  "6 AM",
                                  "12 PM",
                                  "6 PM",
                                  "12 AM")) +
  scale_x_continuous(name = expression(T[surface] - T[bottom]~~(degree~C)),
                     limits = c(-0.5, 4)) +
  scale_y_continuous(name = expression(O[2[surface]] - O[2[bottom]] ~~ (mg~L^{-1}))) +
  #Set plot details
  theme_bw() +
  theme(legend.text.align = 1) +
  #Add ellipse to help guide viewers
  geom_curve(data = data.frame(x = 0.5, y = 7.5, xend = 0, yend = 0.5),
             aes(x = x,
                 y = y,
                 xend = xend,
                 yend = yend),
             arrow = arrow(length = unit(0.03,
                                         "npc"),
                           type = "closed",
                           angle = 50),
             curvature = 0.35,
             alpha = 0.65,
             size = 1.5) +
  geom_curve(data = data.frame(x = 0.15, y = 0.25, xend = 1.75, yend = 3.5),
             aes(x = x,
                 y = y,
                 xend = xend,
                 yend = yend),
             arrow = arrow(length = unit(0.03,
                                         "npc"),
                           type = "closed",
                           angle = 50),
             alpha = 0.65,
             size = 1.5) +
  geom_curve(data = data.frame(x = 1.75, y = 3.65, xend = 0.6, yend = 7.5),
             aes(x = x,
                 y = y,
                 xend = xend,
                 yend = yend),
             arrow = arrow(length = unit(0.03,
                                         "npc"),
                           type = "closed",
                           angle = 50),
             alpha = 0.65,
             size = 1.5)


#----Build-summary-figure----

#Build the figure in two pieces:

#1. The time series panels
time_series_panels <- 
  plot_grid(
    air_water_temperature_ts,
    FF_O2_ts,
    nrow = 2,
    labels = "AUTO",
    align = "hv"
  )

#2. The composite plots
composite_plots <- 
  plot_grid(
    delta_T_delta_O2_plot,
    temperature_daily_cycle_plot,
    ncol = 2,
    align = "h",
    labels = c("C", "D"),
    rel_widths = c(1,1)
  )

#Now combine into a summary figure
lake_overturn_figure <- 
  plot_grid(
    time_series_panels,
    composite_plots,
    nrow = 2,
    rel_heights = c(2,1),
    align = "hv"
  )

#Export
cowplot::ggsave(filename = "figures/figure_S6.pdf",
                plot = lake_overturn_figure,
                device = "pdf",
                height = 10,
                width = 8,
                units = "in")