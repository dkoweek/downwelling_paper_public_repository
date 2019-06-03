# Script to build Figure S4 in "Alleviating hypoxia through induced downwelling"

# Script to build figure of stabily stratified lake based on hydrographic profiles collected in July 2018
# Written by David Koweek on 21 March 2019

#----Initialize-workspace----

#Load necessary packages
library(tidyverse)
library(lubridate)
library(viridis)
library(scales)
library(cowplot)

#Set system clock to avoid Mac OS X / lubridate problems with timestamps
Sys.setenv(TZ = "America/Los_Angeles")

#Load custom functions to unpack MiniDOT data
source("downwelling_field_study_functions.R")

#----Load_MiniDOT_data----

minidot_files <- list.files(path = "data/",
                            pattern = "Cat",
                            full.names = TRUE)
minidot_data <- list()

for (i in 1:length(minidot_files)) {
  minidot_data[[i]] <- 
    read_minidot(minidot_files[i])
}

#Wrangle individual files into single cleaned data frame
minidots <- 
  plyr::ldply(minidot_data,
        data.frame) %>%
  mutate(local_datetime =  with_tz(Coordinated.Universal.Time,
                                   tzone = "America/Los_Angeles")) %>% 
  filter(local_datetime > ymd_hms("2018-07-05 12:45:00", 
                                  tz = "America/Los_Angeles"),
         local_datetime < ymd_hms("2018-07-11 10:30:00",
                                  tz = "America/Los_Angeles")) %>% 
  tbl_df()

#----Plot_MiniDOT_data----

#Temperature time series
minidot_temperature_plot <- 
  minidots %>% 
  ggplot(aes(x = local_datetime,
             y = Temperature)) +
  geom_line(aes(colour = serial_number)) +
  #Set scale details
  scale_colour_viridis(name = element_blank(), 
                       labels = c("2 mab", "0.5 mab", "0 mab"), #checked that serial numbers align to these depths
                       discrete = TRUE,
                       option = "C",
                       direction = -1) +
  scale_y_continuous(name = expression(Temperature~(degree~C)),
                     limits = c(21, 25)) +
  scale_x_datetime(name = element_blank(),
                   date_labels = "%b %d",
                   date_breaks = "1 day") +
  #Set plot details
  theme_bw() +
  theme(axis.text.x = element_text(angle = -45))

#O2 time series
minidot_O2_plot <- 
  minidots %>% 
  ggplot(aes(x = local_datetime,
             y = Dissolved.Oxygen)) +
  geom_line(aes(colour = serial_number)) +
  #Set scale details
  scale_colour_viridis(name = element_blank(), 
                       labels = c("2 mab", "0.5 mab", "0 mab"), #checked that serial numbers align to these depths
                       discrete = TRUE,
                       option = "C",
                       direction = -1) +
  scale_y_continuous(name = expression(O[2]~(mg~L^{-1}))) +
  scale_x_datetime(name = element_blank(),
                   date_labels = "%b %d",
                   date_breaks = "1 day") +
  #Set plot details
  theme_bw() +
  theme(axis.text.x = element_text(angle = -45))

july_2018_hydrography_plot <- 
  plot_grid(minidot_temperature_plot,
            minidot_O2_plot,
            nrow = 2,
            align = "hv",
            labels = "AUTO")

#----Export-plot----

cowplot::ggsave(filename = "figures/figure_S4.pdf",
                plot = july_2018_hydrography_plot,
                device = "pdf",
                width = 8,
                height = 6,
                units = "in")

