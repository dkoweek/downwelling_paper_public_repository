#Run OTE models

#----Initialize_workspace----

#Source scripts to load OTE model functions
source("OTE_functions.R")

#Source scripts to load hydrographic test cases
source("scrape_hypoxia_hydrography.R")

#Load necessary libraries
library(xtable)
library(tools)

#----Define_model_parameter_space----

r_to_d <- c(0.02,0.1)
v <- c(0.01, 0.1) #m/s
eta <- c(2.5, 2.84)
h_max_to_d <- c(0.25, 0.5)
K <- c(0.25, 2) #combined entrance, 90 degree bend, and contraction losses
alpha_pump <- c(0.5, 1)
alpha_air <- seq(0.5, 1, length.out = 64)
phi <- seq(30, 60, length.out = 3) #degrees from horizontal
epsilon_f <- seq(0.6, 30, length.out = 3)
v_f <- seq(5, 10, length.out = 3) #m/s

#Different numbers of parameter values to create ~equal number of model...
#...runs per oxygenation technique (ensures medians are drawn from larger samples)
#Downwelling = 2^6 = 64
#Subsurface aeration = 64
#Surface aeration = 54

parameter_space <- 
  expand.grid(
    r_to_d,
    v,
    eta,
    h_max_to_d,
    K,
    alpha_pump,
    alpha_air,
    phi,
    epsilon_f,
    v_f
  ) %>% 
  rename(r_to_d = Var1,
         v = Var2,
         eta = Var3,
         h_max_to_d = Var4,
         K = Var5,
         alpha_pump = Var6,
         alpha_air = Var7,
         phi = Var8,
         epsilon_f = Var9,
         v_f = Var10)

#----Generate_table_of_parameters----


parameter <-
  c("$\\frac{r}{d}$",
    "$v$",
    "$\\eta$",
    "$\\frac{h_{max}}{d}$",
    "$K_{90}$",
    "$\\alpha_p$",
    "$\\alpha_{air}$",
    "$\\phi$",
    "$\\epsilon_f$",
    "$v_f$"
  )

parameter_description <- 
  c("Pipe radius relative to depth of release (m m$^{-1}$)",
    "Downwelling water velocity (m s$^{-1}$)",
    "Plume entrainment coefficient",
    "Maximum plume height relative to depth of release (m m$^{-1}$)",
    "Minor losses due to pipe bend",
    "Water pump efficiency (\\%)",
    "Air compressor efficiency (\\%)",
    "Fountain ejection angle ($^o$)",
    "Frictional losses in fountain intake",
    "Ejection velocity of fountain (m s$^{-1}$)"
  )

technique <- 
  c("D",
    "D",
    "D",
    "D",
    "D",
    "D,F",
    "I,A",
    "F",
    "F",
    "F")

parameter_range <-
  c(paste(min(r_to_d),"-",max(r_to_d)," (",length(r_to_d),")"),
    paste(min(v),"-",max(v)," (",length(v),")"),
    paste(min(eta),"-",max(eta)," (",length(eta),")"),
    paste(min(h_max_to_d),"-",max(h_max_to_d)," (",length(h_max_to_d),")"),
    paste(min(K),"-",max(K)," (",length(K),")"),
    paste(min(alpha_pump),"-",max(alpha_pump)," (",length(alpha_pump),")"),
    paste(min(alpha_air),"-",max(alpha_air)," (",length(alpha_air),")"),
    paste(min(phi),"-",max(phi)," (",length(phi),")"),
    paste(min(epsilon_f),"-", max(epsilon_f)," (",length(epsilon_f),")"),
    paste(min(v_f),"-", max(v_f)," (",length(v_f),")")
  )


table.df <-
  data.frame(
    cbind(
      parameter,
      parameter_description,
      technique,
      parameter_range
    )
  )

names(table.df) <-
  c(
    "Parameter",
    "Description",
    "Technique",
    "Minimum - Maximum (No. Values)"
  )

parameter_table <- 
    xtable(table.df,
    sanitize.text.function = function(x) {
      x
    },
    comment = FALSE,
    include.rownames = FALSE,
    floating = TRUE,
    latex.environments = "center"
  )

#Export parameter table
latex_parameter_table <- 
  print.xtable(xtable(parameter_table),
               print.results = FALSE,
               sanitize.text.function = function(x) {x},
               comment = FALSE,
               include.rownames = FALSE,
               floating = TRUE,
               latex.environments = "center")

writeLines(
  c(
    "\\documentclass[12pt]{standalone}",
    "\\begin{document}",
    latex_parameter_table,
    "\\end{document}"
  ),
  "parameter_table.tex"
)

try(
  texi2pdf("parameter_table.tex",
           clean = TRUE),
  silent = TRUE
)

#Move table output to the figures folder
file.rename("parameter_table.tex", "figures/parameter_table.tex")
file.rename("parameter_table.pdf", "figures/parameter_table.pdf")

#----Calculate_OTE_downwelling----


downwelling_parameters <- 
  parameter_space %>% 
  distinct(r_to_d,
           v,
           eta,
           h_max_to_d,
           K,
           alpha_pump)

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
                      r_to_d = downwelling_parameters[["r_to_d"]][j],
                      v = downwelling_parameters[["v"]][j],
                      eta = downwelling_parameters[["eta"]][j],
                      h_max_to_d = downwelling_parameters[["h_max_to_d"]][j],
                      K = downwelling_parameters[["K"]][j],
                      alpha_pump = downwelling_parameters[["alpha_pump"]][j])
    
    #Build a data frame of the pertinent results  
    downwelling_results_list[[m]] <- 
      data.frame(Location = i,
                 Depth = data_set[["Depth"]],
                 OTE = downwelling_OTE,
                 r_to_d = downwelling_parameters[["r_to_d"]][j],
                 v = downwelling_parameters[["v"]][j],
                 eta = downwelling_parameters[["eta"]][j],
                 h_max_to_d = downwelling_parameters[["h_max_to_d"]][j],
                 K = downwelling_parameters[["K"]][j],
                 alpha_pump = downwelling_parameters[["alpha_pump"]][j],
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
           v_f)

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
                   v_f = surface_aeration_parameters[["v_f"]][j],
                   phi = surface_aeration_parameters[["phi"]][j])
    
    #Build a data frame of the pertinent results  
    surface_aeration_results_list[[m]] <- 
      data.frame(Location = i,
                 Depth = data_set[["Depth"]],
                 OTE = fountain_OTE,
                 alpha_pump = surface_aeration_parameters[["alpha_pump"]][j],
                 epsilon_f = surface_aeration_parameters[["epsilon_f"]][j],
                 v_f = surface_aeration_parameters[["v_f"]][j],
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

#----Merge_all_model_results----

OTE_model_results <- 
  #Aggregate results from all oxygenation techniques
  bind_rows(downwelling_results,
            subsurface_aeration_results,
            surface_aeration_results) %>% 
  #Reorder techniques so that they follow presentation in the paper
  mutate(technique = factor(technique,
                            levels = c("Downwelling",
                                       "Aeration (isothermal)",
                                       "Aeration (adiabatic)",
                                       "Fountain aeration")))