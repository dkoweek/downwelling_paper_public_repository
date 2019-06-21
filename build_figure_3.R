# Script to build Figure 3 in "Alleviating hypoxia through induced downwelling"

# Apply OTE models to hydrographic data sets and produce figures
# Written by David Koweek on 12 Apr 2019
# Note: almost all of the code here is copied from OTE_modeling.Rmd,...
# ... and is reproduced here for production of publication-quality figure

#----Initialize_workspace----

#Source scripts to load OTE model functions
source("OTE_functions.R")

#Source scripts to load hydrographic test cases
source("scrape_hypoxia_hydrography.R")

#Load necessary packages (many already loaded in previously sourced scripts)

library(viridis)
library(cowplot)
library(ggpubr)
library(xtable)

#----Define_model_parameter_space----

alpha_pump <- seq(0.5, 1, length.out = 4)
epsilon_p <- seq(1e-7, 1e-4, length.out = 4)
alpha_air <- seq(0.5, 1, length.out = 4)
phi <- seq(30, 60, length.out = 4)
epsilon_f <- seq(0.6, 30, length.out = 4)
v <- seq(5, 10, length.out = 4)

parameter_space <- 
  expand.grid(
    alpha_pump,
    epsilon_p,
    alpha_air,
    phi,
    epsilon_f,
    v
  ) %>% 
  rename(alpha_pump = Var1,
         epsilon_p = Var2,
         alpha_air = Var3,
         phi = Var4,
         epsilon_f = Var5,
         v = Var6)

#----Generate_table_of_parameters----


parameter <-
  c(
    "$\\alpha_p$",
    "$\\epsilon_p$",
    "$\\alpha_{air}$",
    "$\\phi$",
    "$\\epsilon_f$",
    "v"
  )

parameter_description <- 
  c(
    "Water pump efficiency (\\%)",
    "Frictional losses during downwelling",
    "Air compressor efficiency (\\%)",
    "Fountain ejection angle ($^o$)",
    "Frictional losses in fountain intake",
    "Ejection velocity of fountain ($m~s^{-1}$)"
  )

parameter_range <-
  c(
    paste(min(alpha_pump),"-",max(alpha_pump)),
    paste(min(epsilon_p),"-",max(epsilon_p)),
    paste(min(alpha_air),"-",max(alpha_air)),
    paste(min(phi),"-",max(phi)),
    paste(min(epsilon_f),"-", max(epsilon_f)),
    paste(min(v),"-", max(v))
  )


table.df <-
  data.frame(
    cbind(
      parameter,
      parameter_description,
      parameter_range
    )
  )

names(table.df) <-
  c(
    "Parameter",
    "Description",
    "Range of Values Considered"
  )

parameter_table <- 
  xtable(table.df,
         caption = "OTE model parameters and range of values considered",
         sanitize.text.function = function(x) {
           x
         },
         comment = FALSE,
         include.rownames = FALSE,
         floating = TRUE,
         latex.environments = "center"
  )


#----Calculate_OTE_downwelling----


downwelling_parameters <- 
  parameter_space %>% 
  distinct(alpha_pump,
           epsilon_p)

downwelling_results_list <- list()
m <- 1 #list index

for (i in hydrocasts_df %>% distinct(Location) %>% pull()) { #For each location
  
  #Separate the data set into each location
  data_set <- 
    hydrocasts_df %>% 
    filter(Location == i)
  
  for(j in 1:nrow(downwelling_parameters)){ #For each location
    
    
    #Calculate the OTE for the water column for each set of parameter values
    downwelling_OTE <- 
      OTE_downwelling(depth = data_set[["Depth"]],
                      temperature = data_set[["Temperature"]],
                      salinity = data_set[["Salinity"]],
                      O2 = data_set[["O2"]],
                      alpha_pump = downwelling_parameters[["alpha_pump"]][j],
                      epsilon_p = downwelling_parameters[["epsilon_p"]][j])
    
    #Build a data frame of the pertinent results  
    downwelling_results_list[[m]] <- 
      data.frame(Location = i,
                 Depth = data_set[["Depth"]],
                 OTE = downwelling_OTE,
                 alpha_pump = downwelling_parameters[["alpha_pump"]][j],
                 epsilon_pump = downwelling_parameters[["epsilon_p"]][j],
                 scenario = m)
    
    m <- 
      m + 1
    
  }
}

#Combine all results into a single data frame
downwelling_results <- 
  plyr::ldply(downwelling_results_list,
              data.frame) %>% 
  tbl_df() %>% 
  mutate(technique = "Downwelling")

#----Calculate_OTE_subsurface_aeration----

subsurface_aeration_parameters <- 
  parameter_space %>% 
  distinct(alpha_air)



isothermal_results_list <- list()
adiabatic_results_list <- list()

m <- 1 #list index

for (i in hydrocasts_df %>% distinct(Location) %>% pull()) { #For each location
  
  #Separate the data set into each location
  data_set <- 
    hydrocasts_df %>% 
    filter(Location == i)
  
  for(j in 1:nrow(subsurface_aeration_parameters)){ #For each location
    
    
    #Calculate the isothermal OTE for the water column for each set of parameter values
    isothermal_OTE <- 
      OTE_isothermal(depth = data_set[["Depth"]],
                     temperature = data_set[["Temperature"]],
                     salinity = data_set[["Salinity"]],
                     alpha_air = subsurface_aeration_parameters[["alpha_air"]][j])
    
    #Build a data frame of the pertinent results from isothermal compression
    isothermal_results_list[[m]] <- 
      data.frame(Location = i,
                 Depth = data_set[["Depth"]],
                 OTE = isothermal_OTE,
                 alpha_air = subsurface_aeration_parameters[["alpha_air"]][j],
                 scenario = m)
    
    #Calculate the adiabatic OTE for the water column for each set of parameter values
    adiabatic_OTE <- 
      OTE_adiabatic(depth = data_set[["Depth"]],
                    temperature = data_set[["Temperature"]],
                    salinity = data_set[["Salinity"]],
                    alpha_air = subsurface_aeration_parameters[["alpha_air"]][j],
                    gamma = 1.4)
    
    #Build a data frame of the pertinent results from adiabatic compression
    adiabatic_results_list[[m]] <- 
      data.frame(Location = i,
                 Depth = data_set[["Depth"]],
                 OTE = adiabatic_OTE,
                 alpha_air = subsurface_aeration_parameters[["alpha_air"]][j],
                 scenario = m)
    
    
    m <- 
      m + 1
    
  }
}

#Combine all isothermal results into a single data frame
isothermal_results <- 
  plyr::ldply(isothermal_results_list,
              data.frame) %>% 
  tbl_df() %>% 
  mutate(technique = "Aeration (isothermal)")

#Combine all adiabatic results into a single data frame
adiabatic_results <- 
  plyr::ldply(adiabatic_results_list,
              data.frame) %>% 
  tbl_df() %>% 
  mutate(technique = "Aeration (adiabatic)")

#Merge isothermal and adiabatic results
subsurface_aeration_results <- 
  bind_rows(isothermal_results,
            adiabatic_results)

#----Calculate_OTE_surface_aeration----


surface_aeration_parameters <- 
  parameter_space %>% 
  distinct(alpha_pump,
           phi,
           epsilon_f,
           v)

surface_aeration_results_list <- list()
m <- 1 #list index

for (i in hydrocasts_df %>% distinct(Location) %>% pull()) { #For each location
  
  #Separate the data set into each location
  data_set <- 
    hydrocasts_df %>% 
    filter(Location == i)
  
  for(j in 1:nrow(surface_aeration_parameters)){ #For each location
    
    
    #Calculate the OTE for the water column for each set of parameter values
    fountain_OTE <- 
      OTE_fountain(depth = data_set[["Depth"]],
                   temperature = data_set[["Temperature"]],
                   salinity = data_set[["Salinity"]],
                   O2 = data_set[["O2"]],
                   alpha_pump = surface_aeration_parameters[["alpha_pump"]][j],
                   epsilon_f = surface_aeration_parameters[["epsilon_f"]][j],
                   v = surface_aeration_parameters[["v"]][j],
                   phi = surface_aeration_parameters[["phi"]][j])
    
    #Build a data frame of the pertinent results  
    surface_aeration_results_list[[m]] <- 
      data.frame(Location = i,
                 Depth = data_set[["Depth"]],
                 OTE = fountain_OTE,
                 alpha_pump = surface_aeration_parameters[["alpha_pump"]][j],
                 epsilon_f = surface_aeration_parameters[["epsilon_f"]][j],
                 v = surface_aeration_parameters[["v"]][j],
                 phi = surface_aeration_parameters[["phi"]][j],
                 scenario = m)
    
    m <- 
      m + 1
    
  }
}

#Combine all results into a single data frame
surface_aeration_results <- 
  plyr::ldply(surface_aeration_results_list,
              data.frame) %>% 
  tbl_df() %>% 
  mutate(technique = "Fountain aeration")

#----Build_individual_plots----

OTE_colour_scale <- 
  inferno(n = 4,
          begin = 0.1,
          end = 0.9)

#Angle Lake plot
angle_lake_OTE_plot <- 
  #Aggregate results from all oxygenation techniques
  bind_rows(downwelling_results,
            subsurface_aeration_results,
            surface_aeration_results) %>% 
  #Reorder techniques so that they follow presentation in the paper
  mutate(technique = factor(technique,
                            levels = c("Downwelling",
                                       "Aeration (isothermal)",
                                       "Aeration (adiabatic)",
                                       "Fountain aeration"))) %>% 
  filter(Location == "Angle Lake, WA, USA") %>% 
  #Plot for all depths except surface level, where downwelling model will give OTE of 0 b/c no O2 gradient present
  filter(Depth > min(Depth)) %>% 
  group_by(Depth, technique) %>% 
  #Calculate ranges for each technique
  dplyr::summarize(max_OTE = max(OTE),
                   min_OTE = min(OTE),
                   mean_OTE = mean(OTE)) %>% 
  ungroup() %>% 
  ggplot() +
  #Create ribbon plot
  geom_ribbon(aes(x = Depth,
                  ymin = min_OTE,
                  ymax = max_OTE,
                  fill = technique),
              alpha = 0.5) +
  #Flip axes to make depth profile
  coord_flip() +
  scale_x_reverse(name = "Depth (m)",
                  limits = c(15, 0)) +
  scale_y_log10(name = expression(Oxygen~Transfer~Efficiency~~(kg~~O[2]~~kWh^{-1})),
                limits = c(10^-6, 10^4),
                breaks = c(10^c(-6:3)),
                labels = c(expression(10^{-6}),
                           "",
                           expression(10^{-4}),
                           "",
                           expression(10^{-2}),
                           "",
                           1,
                           10,
                           100,
                           1000)) +
  scale_fill_manual(name = element_blank(),
                    values = OTE_colour_scale) +
  labs(title = "Angle Lake, WA, USA") +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 10),
        plot.title = element_text(face = "plain")) 


#Beaverdam Reservoir plot  
beaverdam_OTE_plot <- 
  #Aggregate results from all oxygenation techniques
  bind_rows(downwelling_results,
            subsurface_aeration_results,
            surface_aeration_results) %>% 
  #Reorder techniques so that they follow presentation in the paper
  mutate(technique = factor(technique,
                            levels = c("Downwelling",
                                       "Aeration (isothermal)",
                                       "Aeration (adiabatic)",
                                       "Fountain aeration"))) %>% 
  filter(Location == "Beaverdam Res., VA, USA") %>% 
  #Plot for all depths except surface level, where downwelling model will give OTE of 0 b/c no O2 gradient present
  filter(Depth > min(Depth)) %>% 
  group_by(Depth, technique) %>% 
  #Calculate ranges for each technique
  dplyr::summarize(max_OTE = max(OTE),
                   min_OTE = min(OTE),
                   mean_OTE = mean(OTE)) %>% 
  ungroup() %>% 
  #Create ribbon plot
  ggplot() +
  geom_ribbon(aes(x = Depth,
                  ymin = min_OTE,
                  ymax = max_OTE,
                  fill = technique),
              alpha = 0.5) +
  #Flip axes to make depth profile
  coord_flip() +
  scale_x_reverse(name = "Depth (m)",
                  limits = c(15, 0)) +
  scale_y_log10(name = expression(Oxygen~Transfer~Efficiency~~(kg~~O[2]~~kWh^{-1})),
                limits = c(10^-6, 10^4),
                breaks = c(10^c(-6:3)),
                labels = c(expression(10^{-6}),
                           "",
                           expression(10^{-4}),
                           "",
                           expression(10^{-2}),
                           "",
                           1,
                           10,
                           100,
                           1000)) +
  scale_fill_manual(name = element_blank(),
                    values = OTE_colour_scale) +
  labs(title = "Beaverdam Res., VA, USA")+
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 10),
        plot.title = element_text(face = "plain")) 


#Baltic Sea
baltic_sea_OTE_plot <-
  #Aggregate results from all oxygenation techniques
  bind_rows(downwelling_results,
            subsurface_aeration_results,
            surface_aeration_results) %>% 
  #Reorder techniques so that they follow presentation in the paper
  mutate(technique = factor(technique,
                            levels = c("Downwelling",
                                       "Aeration (isothermal)",
                                       "Aeration (adiabatic)",
                                       "Fountain aeration"))) %>% 
  filter(Location == "Boknis Eck, Baltic Sea") %>% 
  #Plot for all depths except surface level, where downwelling model will give OTE of 0 b/c no O2 gradient present
  filter(Depth > min(Depth)) %>% 
  group_by(Depth, technique) %>% 
  #Calculate ranges for each technique
  dplyr::summarize(max_OTE = max(OTE),
                   min_OTE = min(OTE),
                   mean_OTE = mean(OTE)) %>% 
  ungroup() %>% 
  #Create ribbon plot
  ggplot() +
  geom_ribbon(aes(x = Depth,
                  ymin = min_OTE,
                  ymax = max_OTE,
                  fill = technique),
              alpha = 0.5) +
  #Flip axes to make depth profile
  coord_flip() +
  scale_x_reverse(name = "Depth (m)",
                  limits = c(25, 0)) +
  scale_y_log10(name = expression(Oxygen~Transfer~Efficiency~~(kg~~O[2]~~kWh^{-1})),
                limits = c(10^-6, 10^4),
                breaks = c(10^c(-6:3)),
                labels = c(expression(10^{-6}),
                           "",
                           expression(10^{-4}),
                           "",
                           expression(10^{-2}),
                           "",
                           1,
                           10,
                           100,
                           1000)) +
  scale_fill_manual(name = element_blank(),
                    values = OTE_colour_scale) +
  labs(title = "Boknis Eck, Baltic Sea") +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 10),
        plot.title = element_text(face = "plain")) 

#Chesapeake Bay OTE plot
chesapeake_bay_OTE_plot <- 
  #Aggregate results from all oxygenation techniques
  bind_rows(downwelling_results,
            subsurface_aeration_results,
            surface_aeration_results) %>% 
  #Reorder techniques so that they follow presentation in the paper
  mutate(technique = factor(technique,
                            levels = c("Downwelling",
                                       "Aeration (isothermal)",
                                       "Aeration (adiabatic)",
                                       "Fountain aeration"))) %>% 
  filter(Location == "Chesapeake Bay, USA") %>% 
  #Plot for all depths except surface level, where downwelling model will give OTE of 0 b/c no O2 gradient present
  filter(Depth > min(Depth)) %>% 
  group_by(Depth, technique) %>% 
  #Calculate ranges for each technique
  dplyr::summarize(max_OTE = max(OTE),
                   min_OTE = min(OTE),
                   mean_OTE = mean(OTE)) %>% 
  ungroup() %>% 
  #Create ribbon plot
  ggplot() +
  geom_ribbon(aes(x = Depth,
                  ymin = min_OTE,
                  ymax = max_OTE,
                  fill = technique),
              alpha = 0.5) +
  #Flip axes to make depth profile
  coord_flip() +
  scale_x_reverse(name = "Depth (m)",
                  limits = c(25, 0)) +
  scale_y_log10(name = expression(Oxygen~Transfer~Efficiency~~(kg~~O[2]~~kWh^{-1})),
                limits = c(10^-6, 10^4),
                breaks = c(10^c(-6:3)),
                labels = c(expression(10^{-6}),
                           "",
                           expression(10^{-4}),
                           "",
                           expression(10^{-2}),
                           "",
                           1,
                           10,
                           100,
                           1000)) +
  scale_fill_manual(name = element_blank(),
                    values = OTE_colour_scale) +
  labs(title = "Chesapeake Bay, USA") +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 10),
        plot.title = element_text(face = "plain")) 

# northern Gulf of Mexico
nGOM_OTE_plot <- 
  #Aggregate results from all oxygenation techniques
  bind_rows(downwelling_results,
            subsurface_aeration_results,
            surface_aeration_results) %>% 
  #Reorder techniques so that they follow presentation in the paper
  mutate(technique = factor(technique,
                            levels = c("Downwelling",
                                       "Aeration (isothermal)",
                                       "Aeration (adiabatic)",
                                       "Fountain aeration"))) %>% 
  filter(Location == "Gulf of Mexico, USA") %>% 
  #Plot for all depths except surface level, where downwelling model will give OTE of 0 b/c no O2 gradient present
  filter(Depth > min(Depth)) %>% 
  group_by(Depth, technique) %>% 
  #Calculate ranges for each technique
  dplyr::summarize(max_OTE = max(OTE),
                   min_OTE = min(OTE),
                   mean_OTE = mean(OTE)) %>% 
  ungroup() %>% 
  #Create ribbon plot
  ggplot() +
  geom_ribbon(aes(x = Depth,
                  ymin = min_OTE,
                  ymax = max_OTE,
                  fill = technique),
              alpha = 0.5) +
  #Flip axes to make depth profile
  coord_flip() +
  scale_x_reverse(name = "Depth (m)",
                  limits = c(25, 0)) +
  scale_y_log10(name = expression(Oxygen~Transfer~Efficiency~~(kg~~O[2]~~kWh^{-1})),
                limits = c(10^-6, 10^4),
                breaks = c(10^c(-6:3)),
                labels = c(expression(10^{-6}),
                           "",
                           expression(10^{-4}),
                           "",
                           expression(10^{-2}),
                           "",
                           1,
                           10,
                           100,
                           1000)) +
  scale_fill_manual(name = element_blank(),
                    values = OTE_colour_scale) +
  labs(title = "Gulf of Mexico, USA") +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 10),
        plot.title = element_text(face = "plain")) 



#----Build_aggregate_figure----

#First grab legend
legend <- 
  get_legend(
    baltic_sea_OTE_plot + theme(
      legend.title.align = 0.5,
      legend.direction = "vertical",
      legend.key.width = unit(0.75, "in"),
      legend.justification = "center"
    )
  ) %>% as_ggplot()

#Combine all figure-specific sites and use common legend

all_sites_OTE_plot <-
  plot_grid(
    baltic_sea_OTE_plot + theme(legend.position = "none"),
    chesapeake_bay_OTE_plot + theme(legend.position = "none"),
    nGOM_OTE_plot + theme(legend.position = "none"),
    angle_lake_OTE_plot + theme(legend.position = "none"),
    beaverdam_OTE_plot + theme(legend.position = "none"),
    legend,
    ncol = 3,
    align = "hv",
    labels = c("A", "B", "C", "D", "E", "")
  )

#----Export_aggregate_figure----

cowplot::ggsave(filename = "figures/figure_3.pdf",
                plot = all_sites_OTE_plot,
                height = 6.5,
                width = 9.5,
                units = "in")