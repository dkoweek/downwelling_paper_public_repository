# Script to build Figure S9 in "Alleviating hypoxia through induced downwelling"

#Script to compile figure of O2 times series data to show raw experimental effect
#Written by David Koweek 16 May 2019

#----Initialize_workspace----

#Load relevant packages
library(tidyverse)
library(lubridate)

#Set system clock to avoid Mac OS X / lubridate problems with timestamps
Sys.setenv(TZ = "America/Los_Angeles")

#Load primary data
source("load_primary_data.R")

#Load sensor locations log
source("load_field_notes.R")



#----Prepare_for_plotting----

#Build a character vector lookup table to return sensor depth and angle based on serial number

sensor_lookup_table <- 
  setNames(
    paste0(sensor_locations_log$Meters_above_bottom,
           " mab ",
           sensor_locations_log$Angle_CW_from_N,
           "°"), 
    as.character(sensor_locations_log$Serial_Number))

pumping_times <- c(18.25, 21.25) #average start ~6:15pm, average finish 3 hours later at ~9:15pm

#----Build_plot----
O2_time_series_plot <-
  oxygenation_expt_sensor_data %>% 
  #Take bottom sensors near the pipe
  filter(Radius == 2,
         Meters_above_bottom < 1) %>%
  #Calculate the time of day
  mutate(local_date = date(local_datetime),
         local_time_of_day = interval(ymd_hms(str_c(local_date, "00:00:00"),
                                              tz = "PDT"),
                                      local_datetime) %>% int_length(),
         local_time_of_day = local_time_of_day / 3600 #seconds -> hrs
         ) %>% 
  ggplot(aes(x = local_time_of_day,
             y = O2)) +
  geom_path(aes(group = local_date,
                colour = as.factor(local_date)),
            show.legend = FALSE) +
  scale_y_continuous(name = expression(O[2]~(mg~L^{-1}))) +
  scale_x_continuous(name = "Hour of day (PDT)") +
  annotate(geom = "rect",
           xmin = pumping_times[1],
           xmax = pumping_times[2],
           ymin = 0,
           ymax = 10,
           colour = "black",
           fill = NA) +
  facet_wrap(~serial_number,
             labeller = labeller(serial_number = sensor_lookup_table))
  

#----Export_plot----
ggplot2::ggsave(filename = "figures/figure_S9.pdf",
                plot = O2_time_series_plot,
                device = "pdf",
                height = 5,
                width = 8,
                units = "in")
