# Script to build Figure S6 in "Alleviating hypoxia through induced downwelling"

#----Initialize_workspace----

#Load relevant packages
library(tidyverse)
library(lubridate)
library(magrittr)
library(broom)
library(caTools)
library(viridis)
library(cowplot)

#Load primary data
source("load_primary_data.R")

#Load data from field notebooks
source("load_field_notes.R")

#Set system clock to avoid Mac OS X / lubridate problems with timestamps
Sys.setenv(TZ = "America/Los_Angeles")


#----Exploratory_data_visualization_of_lake_overturn----

#Define the offset in time from local time when the lake overturns
overturn_offset <- -7 #hours into local day when lake overturns


oxygenation_expt_sensor_data %>% 
  filter(Mooring == "FF",
         Meters_above_bottom < 0.5) %>% 
  mutate(overturn_datetime = local_datetime + hours(overturn_offset),
         overturn_date = date(overturn_datetime),
         overturn_time_of_day = interval(ymd_hms(str_c(overturn_date, "00:00:00"),
                                                 tz = "PDT"),
                                         overturn_datetime) %>% 
           int_length()) %>% 
  ggplot(aes(x = overturn_time_of_day,
             y = O2)) +
  geom_line() +
  labs(x = str_c("PDT ",
                 overturn_offset,
                 " hours (seconds)"),
       y = expression(O[2]~(mg~L^{-1})),
       title = str_c("Far field bottom sensor")) +
  facet_wrap(~overturn_date) 

# Now let's see if we can remove the "edges" of the overturn daily cycle for a cleaner signal


start_of_overturn_day_hrs_to_eliminate <- 2 
end_of_overturn_day_hrs_to_eliminate <- 20


oxygenation_expt_sensor_data %>% 
  filter(Mooring == "FF",
         Meters_above_bottom < 0.5) %>% 
  mutate(overturn_datetime = local_datetime + hours(overturn_offset),
         overturn_date = date(overturn_datetime),
         overturn_time_of_day = interval(ymd_hms(str_c(overturn_date, "00:00:00"),
                                                 tz = "PDT"),
                                         overturn_datetime) %>% 
           int_length()) %>% 
  #Filter out the edges of each overturn day
  filter(hour(overturn_datetime) > start_of_overturn_day_hrs_to_eliminate &
           hour(overturn_datetime) < end_of_overturn_day_hrs_to_eliminate) %>% 
  #Re-generate the plot
  ggplot(aes(x = overturn_time_of_day,
             y = O2)) +
  geom_line() +
  labs(x = str_c("PDT ",
                 overturn_offset,
                 " hours (seconds)"),
       y = expression(O[2]~(mg~L^{-1})),
       title = str_c("Far field bottom sensor")) +
  facet_wrap(~overturn_date) 

#----Wrangle_new_data_frame_of_O2_drawdown_based_on_exploratory_plots----

respiration_df <- 
  oxygenation_expt_sensor_data %>% 
  filter(Mooring == "FF",
         Meters_above_bottom < 0.5) %>% 
  mutate(overturn_datetime = local_datetime + hours(overturn_offset),
         overturn_date = date(overturn_datetime),
         overturn_time_of_day = interval(ymd_hms(str_c(overturn_date, "00:00:00"),
                                                 tz = "PDT"),
                                         overturn_datetime) %>% 
           int_length()) %>% 
  filter(overturn_date <= "2018-09-29", #Filter out the un-stratified days
         #Filter out the edges of each overturn day
         hour(overturn_datetime) > start_of_overturn_day_hrs_to_eliminate & 
           hour(overturn_datetime) < end_of_overturn_day_hrs_to_eliminate) %>% 
  #Remove a bit more data on the "overturn" day of 15-September when overturn happened earlier
  mutate(O2 = case_when(overturn_date == "2018-09-15" & hour(overturn_datetime) > 18 ~ NA_real_,
                        TRUE ~ O2))

#----Generate_plot_of_daily_O2_drawdown----

bottom_water_O2_plot <- 
  respiration_df %>% 
  ggplot(aes(x = overturn_time_of_day / 3600, #convert time scale from seconds to hours
             y = O2)) +
  geom_line() +
  labs(x = str_c("PDT ",
                 overturn_offset,
                 " (hours)"),
       y = expression(O[2]~(mg~L^{-1})),
       title = str_c("Far field bottom sensor (0 mab)")) +
  facet_wrap(~overturn_date) 

#----Calculate_respiration_rates----

respiration_estimates <- 
  respiration_df %>% 
  #Apply a linear regression to O2 as a function of time of day for each overturn date
  nest(-overturn_date) %>% 
  #Tidy the data
  mutate(fit = map(data, ~lm(O2 ~ overturn_time_of_day,
                             data = .x)),
         cleaned_fits = map(fit, tidy)) %>% 
  unnest(cleaned_fits) %>% 
  slice(seq(2, nrow(.), by = 2)) %>% 
  #Convert the units
  mutate(respiration_per_hr = estimate * 3600, #mg O2/L/s -> mg O2/L/hr
         respiration_per_hr_e  = std.error * 3600)

#----Plot_respiration_rates----

respiration_rates <- 
  respiration_estimates %>% 
  mutate(respiration_per_hr = abs(respiration_per_hr)) %>% 
  ggplot(aes(x = overturn_date)) +
  geom_bar(aes(y = respiration_per_hr),
           stat = "identity",
           fill = "red") +
  geom_errorbar(aes(ymin = respiration_per_hr - respiration_per_hr_e,
                    ymax = respiration_per_hr + respiration_per_hr_e),
                width = 0.4) +
  scale_x_date(breaks = scales::date_breaks("1 day"),
               labels = scales::date_format("%b-%d"),
               expand = c(0,0)) +
  geom_hline(aes(yintercept =  mean(abs(respiration_per_hr),
                                    na.rm = TRUE)),
             linetype = "dashed",
             alpha = 0.5) +
  #Plot details
  labs(x = "Start date of measurement (PDT)",
       y = expression(Respiration~Rate~(mg~O[2]~L^{-1}~hr^{-1})),
       title = expression(Far~~Field~~Respiration~~Rate~(mean~"\U00B1"~S.E.))) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = -45,
                                   vjust = -0.2))

#----Produce_multipanel_figure----


lake_respiration_figure <- 
  plot_grid(
    bottom_water_O2_plot,
    respiration_rates,
    nrow = 2,
    labels = "AUTO",
    align = "v"
  )


cowplot::ggsave(filename = "figures/figure_S7.pdf",
                plot = lake_respiration_figure,
                device = "pdf",
                height = 10,
                width = 8,
                units = "in")
