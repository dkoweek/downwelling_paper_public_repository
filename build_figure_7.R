# Script to build Figure 7 in "Alleviating hypoxia through induced downwelling"

# Back-of-the-envelope calculations to understand scaling from OTE results to real world systems
# Written by David Koweek on 30 April 2019
# Updated on 6 May 2019
# Key idea: Downwelling > Respiration (per m^2) / water depth (m)


#----Initialize_workspace----

#Load necessary packages
library(tidyverse)
library(viridis)

#----Set_model_parameters----

#Water body dimensions
depth <- 10 #m

log_V <- seq(6, 11, length.out = 10) #log 10 of m^3

m3_per_km3 <- 1e9

#Respiratory demand
RQ <- 1 # delta CO2/-delta O2

R_CO2 <- 300 #mmol C/m^2/d (~average estuarine respiration condition (Hopkinson and Smith 2005))

R_O2 <- R_CO2 / RQ #mmol O2/m^2/d

R_O2 <- R_O2 / 1000 / 24 #mol O2/m^2/hr

R_O2_depth_avg <- R_O2 / depth #mol O2/m^3/hr

#Downwelling conditions
OTE <- seq(10, 100, by = 10) #kg O2/kWh (~upper and lower bounds of OTE modeling exercise)

O2_molar_mass <- 32/1000 #kg O2/mol O2

OTE_molar <- OTE / O2_molar_mass #mol O2/kWh

#Cost of energy

USD_kWh <- .125 #USD/kWh
USD_MWh <- USD_kWh * 1000 #USD/MWh

hours_per_year <- 8760 #hr/year

dollars_per_million <- 1e6

# Estuarine test cases
nGOM_hypoxic_volume <- 60 #km^3 (5-year areal average of ~15000 km^2 -> 60 km^3 in Scavia et al. 2018)

chesapeake_bay_hypoxic_volume <- 10.8 #km^3

#----Calculate_model_conditions----
downwelling_scaled_df <-
  expand.grid(
    log_V,
    OTE_molar
  ) %>% 
  rename(log_V = Var1,
         OTE_molar = Var2) %>% 
  mutate(V = 10^log_V, #m^3
         unit_power = R_O2_depth_avg / OTE_molar, #kW/m^3
         unit_cost = (unit_power * m3_per_km3 / 1000) * USD_MWh, # kW/m^3 * m^3/km^3 * MW/kW * USD/MW*h = USD/km^3*h
         power = unit_power * V,  #kW
         power = power / 1000 #MW
        ) 

#----Plot_data----
downwelling_scaled_plot <- 
  downwelling_scaled_df %>% 
  ggplot(aes(x = V / m3_per_km3,
             y = power)) +
  geom_path(aes(colour = as.factor(OTE_molar))) +
  scale_x_continuous(name = expression(Water~Volume~(km^3))) +
  scale_y_continuous(name = expression(Power~(MW)),
                     sec.axis = sec_axis(~.*USD_MWh * hours_per_year / dollars_per_million,
                                         name = paste("Million USD/year (assuming $", USD_MWh, " per MWh)", sep = ""))) +
  scale_colour_viridis(option = "B",
                       discrete = TRUE,
                       begin = 0.2,
                       end = 0.8,
                       name = expression(Oxygen~Transfer~Efficiency~(kg~O[2]~kWh^{-1})),
                       breaks = OTE_molar,
                       labels = as.character(OTE)) +
  #Add annotations for hypoxic water bodies
  #Northern Gulf of Mexico
  geom_vline(xintercept = nGOM_hypoxic_volume,
             linetype = "dashed") +
  annotate("text",
           x = 75,
           y = 3.75e2,
           label = expression(atop("nGOM \"Dead Zone\"", "2014-2018 average")),
           size = 3) +
  geom_segment(aes(x = 72.5,
                   xend = 62,
                   y = 3.5e2,
                   yend = 3e2),
               arrow = arrow(angle = 15,
                             length = unit(0.03, "npc")),
               size = 0.25) +
  #Chesapeake Bay
  geom_vline(xintercept = chesapeake_bay_hypoxic_volume,
             linetype = "dashed") +
  annotate("text",
           x = 30,
           y = 2.4e2,
           label = expression(atop(atop("Chesapeake Bay", "\"Dead Zone\""), atop("2018 Maximum", "Daily Volume"))),
           size = 4) +
  geom_segment(aes(x = 25,
                   xend = 12.5,
                   y = 2e2,
                   yend = 1.5e2),
               arrow = arrow(angle = 15,
                             length = unit(0.03, "npc")),
               size = 0.25) +
  theme_bw() +
  theme(legend.direction = "horizontal", 
        legend.position = "bottom")

#----Export_plot----

ggsave(filename = "figures/figure_7.pdf",
       plot = downwelling_scaled_plot,
       device = "pdf",
       width = 8,
       height = 5,
       units = "in")


