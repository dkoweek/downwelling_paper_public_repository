# Script to build Figure 5 in "Alleviating hypoxia through induced downwelling"

#----Initialize_workspace----

#Load relevant packages
library(tidyverse)
library(lubridate)
library(viridis)
library(gridExtra)
library(ggpubr)

#Set system clock to avoid Mac OS X / lubridate problems with timestamps
Sys.setenv(TZ = "America/Los_Angeles")

#Load primary data set to match up treatment ID to local dates
source("load_primary_data.R")

#Load data from field notebooks
source("load_field_notes.R")


#----Load_data_sets_of_interpolated_results----
#Load data files of temporal experimental effects considering 30, 60, and 120 minutes before/after pumping periods

interpolation_results <- list()

interpolation_results[[1]] <- 
  read_csv(file = "data/mapped_pumping_data_half_hour.csv")

interpolation_results[[2]] <- 
  read_csv(file = "data/mapped_pumping_data_1h.csv")

interpolation_results[[3]] <- 
  read_csv(file = "data/mapped_pumping_data_2h.csv")

#----Wrangle_data_sets_for_analysis_and_plotting----

cleaned_results <- list()
for (i in 1:length(interpolation_results)) {
  
  cleaned_results[[i]] <- 
    interpolation_results[[i]] %>% 
    #Rename some columns in preparation for merge with other sensor data
    rename(Serial_Number = Sensor) %>% 
    select(c(Serial_Number:mab), contains("dev")) %>% 
    #Convert wide data to tall data
    gather("Date", "O2_prime", -Serial_Number, -x_pos, -y_pos, -mab) %>% 
    group_by(Serial_Number) %>% 
    nest() %>% 
    #Merge with sensor location log
    full_join(sensor_locations_log,
              by = "Serial_Number") %>% 
    #Only look at sensors in the near field
    filter(Radius %in% c(2,5,8)) %>% 
    #Add in column to define interpolation scenario
    mutate(interp_length = 30 * 2^(i - 1)) %>%  #30 min, 60 min, 120 min
    unnest() %>% 
    #Determine whether a day was pumping or not pumping
    mutate(treatment_ID_local = str_extract(string = Date, 
                                            pattern = "[CP]\\d+_PDT")) %>% 
    #Create a column to denote whether a day was a pumping or control day.
    mutate(pumping = case_when((str_detect(treatment_ID_local, "^C") == TRUE) ~ "Control",
                               TRUE ~ "Pumping")) %>% 
    #Extract whether the interpolate was for the full 3 hours, the 1st hour only, the 2nd hour only, or the 3rd hour only
    mutate(effect_window = case_when(str_detect(Date, "dev.$") == TRUE ~ "Full",
                                     str_detect(Date, "1$") == TRUE ~ "Hour 1",
                                     str_detect(Date, "2$") == TRUE ~ "Hour 2",
                                     str_detect(Date, "3$") == TRUE ~ "Hour 3")) %>% 
    #Get the date of each treatment by joining with the primary data set
    left_join(.,
              oxygenation_expt_sensor_data %>% 
                mutate(local_date = date(local_datetime)) %>% 
                select(c(treatment_ID_local, local_date)) %>% 
                distinct(),
              by = "treatment_ID_local") %>% 
    #Remove the day at the very end of the experiment (after sensors were out of the water) without a treatment ID
    filter(!is.na(treatment_ID_local))
  
}

#Merge the 30-min, 60-min, and 120-min results into a single data frame
experimental_effects <- 
  plyr::ldply(cleaned_results,
              data.frame) %>% 
  tbl_df() 

#----Summarize_results----

#Define max height of bottom of water column
bottom_cutoff <- 1.5 #m

experimental_summary <- 
  experimental_effects %>% 
  mutate(depth_group = case_when((Meters_above_bottom <= bottom_cutoff) ~ "Bottom",
                                 TRUE ~ "Surface")) %>% 
  #Summarize by depth_group at each mooring for each day for each effect for each interpolation length
  group_by(depth_group, Radius, Angle_CW_from_N, local_date, pumping, effect_window, interp_length) %>% 
  summarize(O2_prime_bar = mean(O2_prime, na.rm = TRUE)) %>% 
  ungroup()


#----Justify_choice_of_interpolation_window----

error_quantified <-
  experimental_summary %>% 
  group_by(pumping, interp_length, depth_group) %>% 
  summarize(ME = mean(O2_prime_bar, na.rm = TRUE),
            MAE = mean(abs(O2_prime_bar), na.rm = TRUE)) %>% 
  filter(depth_group == "Bottom",
         pumping == "Control")

print(error_quantified)

#Consider bottom water results for the full pumping duration based on the regressions...
#...with 120-min windows on each side of the pumping window

#Showing the regression effects for the full 3-hours pumping duration only, instead...
#... of the 1st, 2nd, or 3rd hours only because the results do not change qualitatively.
# Interested users can change the 'effect_window' filter if they are interested in seeing the results...
#... from the 1st, 2nd, or 3rd hours alone.

expt_results <- 
  experimental_summary %>% 
  filter(interp_length == 120,
         depth_group == "Bottom",
         effect_window == "Full")

#----Prepare_to_generate_plot_of_results----

#Build a character vector lookup table to return sensor depth and angle based on serial number
heading_lookup_table <-
  setNames(
    c(
      "60° (NE)",
      "180° (S)",
      "300° (NW)",
      "0° (N)",
      "120° (SE)",
      "240° (SW)"
    ),
    experimental_summary %>% 
      pull(Angle_CW_from_N) %>% 
      unique()
  )

#Build custom function for labelling radial distances           
radial_labeller <- 
  function(radius) paste0("r=",radius,"m")


#Add mooring letter to experimental results before plotting
expt_results <-
  expt_results %>% 
  group_by(Radius, Angle_CW_from_N) %>% 
  nest() %>% 
  left_join(.,
            sensor_locations_log %>% 
              select(c(Mooring, 
                       Radius,
                       Angle_CW_from_N)) %>% 
              distinct()) %>% 
  unnest() %>% 
  ungroup()

#----Generate_plot_of_experimental_results----

experimental_results_summary_figure <- 
  expt_results %>% 
  ggplot(aes(x = local_date)) +
  geom_bar(aes(y = O2_prime_bar,
               fill = pumping),
           colour = "black",
           size = 0.25,
           stat = "identity") +
  scale_fill_viridis(discrete = TRUE,
                     name = element_blank()) +
  scale_x_date(name = element_blank(),
               date_labels = "%b %d") +
  scale_y_continuous(name = expression(paste(O[2],"`",~(mg~L^{-1})))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 0.6)) +
  facet_grid(Radius ~ Angle_CW_from_N,
             labeller = labeller(Radius = as_labeller(radial_labeller),
                                 Angle_CW_from_N = heading_lookup_table)) 

#Add labels for each mooring
experimental_results_summary_figure <- 
  experimental_results_summary_figure +
  geom_text(data = expt_results,
            aes(x = ymd("2018-09-15"),
                y = 3.5,
                label = Mooring),
            family = "Helvetica",
            nudge_x = 0.5) 

#Temporarily remove color legend
temp <- 
  experimental_results_summary_figure + theme(legend.position = "none")

#Place colour legend on bottom of plot
experimental_results_summary_figure <- 
  temp %>%
  annotate_figure(bottom = get_legend(
    experimental_results_summary_figure + theme(
      legend.title.align = 0.5,
      legend.direction = "horizontal",
      legend.key.width = unit(0.75, "in"),
      legend.justification = "center"
    )
  ))

ggplot2::ggsave(filename = "figures/figure_5.png",
                plot = experimental_results_summary_figure,
                height = 6,
                width = 8,
                units = "in")

#----Hypothesis_testing_for_differences_between_control_and_pumping_days----


#Define the suite of radial/angle/effect window combinations that need to be tested
radius_angle_effect_window_combinations <- 
  expt_results %>% 
  distinct(Radius, 
           Angle_CW_from_N, 
           effect_window)

#Pre-allocate lists of results
wilcox_results <- list()
t_test_results <- list()

for (i in 1:nrow(radius_angle_effect_window_combinations)) {
  
  #Chunk the data by heading/radius
  mooring_bottom_data  <- 
    expt_results %>% 
    filter(Radius == radius_angle_effect_window_combinations[["Radius"]][i],
           Angle_CW_from_N == radius_angle_effect_window_combinations[["Angle_CW_from_N"]][i],
           effect_window == radius_angle_effect_window_combinations[["effect_window"]][i],
           O2_prime_bar != "NaN")
  
  #Perform a Wilcoxon Rank Sum test
  wilcox_results[[i]] <- 
    wilcox.test(O2_prime_bar ~ pumping,
                data = mooring_bottom_data,
                paired = FALSE,
                alternative = "less") #First group in "pumping" is control, so alternative hypothesis...
  #...is that effects on control days are less than effects on experimental days.
  
  #Perform a t-test
  t_test_results[[i]] <- 
    t.test(O2_prime_bar ~ pumping,
           data = mooring_bottom_data,
           paired = FALSE,
           alternative = "less") #First group in "pumping" is control, so alternative hypothesis...
  #...is that effects on control days are less than effects on experimental days.
  
}



test_summary <- 
  radius_angle_effect_window_combinations %>% 
  mutate(Wilcoxon_p_value = sapply(wilcox_results, function(x) x[["p.value"]]),
         Wilcoxon_significant = case_when((Wilcoxon_p_value < 0.05) ~ "significant",
                                          TRUE ~ "not significant"),
         t_test_p_value = sapply(t_test_results, function(x) x[["p.value"]]),
         t_test_significant = case_when((t_test_p_value < 0.05) ~ "significant",
                                        TRUE ~ "not significant")) %>% 
  mutate_if(is.numeric, function(x) signif(x, digits = 2)) 

print.data.frame(
  test_summary %>%
    filter(
      Wilcoxon_significant == "significant" |
        t_test_significant == "significant"
    )
)
