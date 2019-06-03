# Functions necessary to calculate oxygen transfer efficiency for downwelling, subsurface aeration, and surface aeration
# Written by David Koweek on 6 February 2019


#----Downwelling_functions----

#Delta~h (pump head height)
head_height <- function(depth, temperature, salinity, epsilon_p) {
  
  #Seawater density (kg/m^3)
  rho <-
    oce::swRho(salinity = salinity,
               temperature = temperature,
               pressure = 0)
  
  #Mean density (kg/m^3)
  rho_bar <-
    cumsum(rho) / seq_along(rho)
  
  #Surface water density (kg/m^3)
  rho_0 <- 
    data.frame(
      depth,
      rho
    ) %>% 
    arrange(depth) %>% #ensures that data frame is sorted and 1st depth value is retrieved
    top_n(-1,depth) %>% 
    pull(rho)
  
  #Calculated head height (m)
  delta_h <-
    depth * (
              (rho_bar / rho_0) -
               1 +
               epsilon_p
             )
    
  return(delta_h)
    
}

#Oxygen gradient between the surface and release point
delta_O2_pump <- function(depth, O2) {
  
  O2_0 <- 
    data.frame(
      depth,
      O2
    ) %>% 
    arrange(depth) %>% #ensures that data frame is sorted and 1st depth value is retrieved
    top_n(-1, depth) %>% 
    pull(O2)
  
  delta_O2_pump <- 
    O2_0 - O2
  
  return(delta_O2_pump)
  
}

#Oxygen transfer efficiency of downwelling
OTE_downwelling <- function(depth, temperature, salinity, O2, epsilon_p = 0, alpha_pump = 1) {
  
  #Define constant values
  g <- 9.81 #gravitational constant (m/s^2)
  
  O2_MM <- 32 #O2 molar mass (g/mol O2)
  
  J_to_kWh <- 3.6e6 #Joules per kWh
  
  dh <- #m
    head_height(depth = depth,
                temperature = temperature,
                salinity = salinity,
                epsilon_p = epsilon_p)
  
  dO2 <- #mol/kg
    delta_O2_pump(depth = depth,
                  O2 = O2)
  
  OTE <- #mol O2/J
    (alpha_pump * dO2) /
    (g * dh)
  
  OTE <- #kg O2/kWh
    OTE * 
    (O2_MM / 1000) * #mol O2 -> kg O2
    J_to_kWh #J -> kWh
      
  return(OTE)
  
}

#----Subsurface_aeration_functions----
OTE_isothermal <- function(depth, temperature, salinity, alpha_air = 1) {
  
  #Define constant values
  x_O2 <- 0.21 #mole fraction of O2 in air
  
  R <- 8.3144598 #Ideal gas constant (J/mol/K)
  
  g <- 9.81 #gravitational constant (m/s^2)
  
  Pa_per_atm <- 101325 #Pa/atm #convert Pa to atm
  
  P_surface <- 1 * Pa_per_atm
  
  J_to_kWh <- 3.6e6 #Joules per kWh
  
  #Temperature (K)
  T_K <- 
    temperature + 273.15

  #Seawater density (kg/m^3)
  rho <-
    oce::swRho(salinity = salinity,
               temperature = temperature,
               pressure = 0)
  
  #Mean density (kg/m^3)
  rho_bar <-
    cumsum(rho) / seq_along(rho)
  
  OTE <- #mol O2/J
    (alpha_air * x_O2) /
    (R * T_K *
       log(1 +
             ((rho_bar * g * depth) /
                (P_surface))))     
  
  OTE <- #kg O2/kWh
    OTE * 
    (O2_MM / 1000) * #mol O2 -> kg O2
    J_to_kWh #J -> kWh
         
  return(OTE)
  
}

OTE_adiabatic <- function(depth, temperature, salinity, alpha_air = 1, gamma = 1.4) {
  
  #Define constant values
  x_O2 <- 0.21 #mole fraction of O2 in air
  
  R <- 8.3144598 #Ideal gas constant (J/mol/K)
  
  g <- 9.81 #gravitational constant (m/s^2)
  
  Pa_per_atm <- 101325 #Pa/atm #convert Pa to atm
  
  P_surface <- 1 * Pa_per_atm
  
  J_to_kWh <- 3.6e6 #Joules per kWh
  
  #Surface Temperature
  T_0 <- #C
    data.frame(
      depth,
      temperature
    ) %>% 
    arrange(depth) %>% #ensures that data frame is sorted and 1st depth value is retrieved
    top_n(-1, depth) %>% 
    pull(temperature)
    
  T_K_0 <- #K
    T_0 + 273.15
  
  #Seawater density (kg/m^3)
  rho <-
    oce::swRho(salinity = salinity,
               temperature = temperature,
               pressure = 0)
  
  #Mean density (kg/m^3)
  rho_bar <-
    cumsum(rho) / seq_along(rho)
  
  OTE <- #mol O2/J
    ((alpha_air * (gamma - 1) * x_O2) /
       (R * T_K_0)) *
    ((1 + ((rho_bar * g * depth) / P_surface)) ^ ((gamma - 1) / gamma) - 1) ^ -1
       
  
  OTE <- #kg O2/kWh
    OTE * 
    (O2_MM / 1000) * #mol O2 -> kg O2
    J_to_kWh #J -> kWh
  
  return(OTE)
  
}
#----Surface_aeration_functions----

#Oxygen solubility (#mol/kg/atm)
S_O2 <- 
  function(temperature, salinity) {
    
    #Seawater density (kg/m^3)
    rho <-
      oce::swRho(salinity = salinity,
                 temperature = temperature,
                 pressure = 0)
    
    #Bunsen solubility coefficients for O2 (from Weiss 1970)
    T_K <-
      temperature + 273.15 #deg K
    
    A_O2 <- c(-58.3877, 85.8079, 23.8439)
    B_O2 <- c(-0.034892, 0.015568, -0.0019387)
    Bunsen_O2 <-
      exp((A_O2[1] + (A_O2[2] * (100 / T_K)) + (A_O2[3] * (log(
        T_K / 100
      )))) +
        (salinity * (B_O2[1] + (B_O2[2] * (
          T_K / 100
        )) + (B_O2[3] * ((T_K / 100) ^ 2
        )))))
    V_bar_ideal <- 22.4136 #L/mol
    S_O2_vol <- (Bunsen_O2 / V_bar_ideal) * 1e3 #mol/m^3/atm
    
    S_O2 <- 
      (S_O2_vol / rho)
    
    return(S_O2)
  }

#Oxygen deficit between observed concentration and equilibrium concentration
O2_deficit <- function(temperature, salinity, O2) {
  
  #Define constant values
  x_O2 <- 0.21 #atm O2/atm air
  
  P_surface <- 1 #atm of air
  
  
  O2_eq <- #mol O2/kg SW
    S_O2(temperature = temperature,
         salinity = salinity) *
    x_O2 *
    P_surface
  
  
  O2_prime <- #mol/kg O2
    O2 - O2_eq
  
  return(O2_prime)
  
}

#Mass transfer rate constant
k_rate <- function(v) {
  A <- 2e-4 #s/m^2
  
  k <- #s^-1
    A * v^2 
  
  return(k)
  
}

delta_O2_fountain <- function(depth, temperature, salinity, O2, v, phi) {
  
  #Define constant values
  g <- 9.81 #gravitational constant (m/s^2)
  
  phi <- #convert degrees to radians
    phi * (pi / 180)
  
  O2_prime <- #mol/kg
    O2_deficit(temperature = temperature,
               salinity = salinity,
               O2 = O2)

  k <- #s^-1
    k_rate(v)
  
  delta_O2 <- 
    O2_prime *
    (exp((-2 * k * v * sin(phi)) / g) - 1)
  
  return(delta_O2)
  
  
}

OTE_fountain <- function(depth, temperature, salinity, O2, v, phi = 60, alpha_pump = 1, epsilon_f = 6) {
  
  #Define constant values
  g <- 9.81 #gravitational constant (m/s^2)
  
  O2_MM <- 32 #O2 molar mass (g/mol O2)
  
  J_to_kWh <- 3.6e6 #Joules per kWh
  
  delta_O2 <- #mol O2/kg
    delta_O2_fountain(depth = depth,
                      temperature = temperature,
                      salinity = salinity,
                      O2 = O2,
                      v = v,
                      phi = phi)
  
  J_per_kg <- #Power per kg of water
    (depth * g) + 
      (
        (v^2 /2) *
         (epsilon_f + 1)
       )
  
  OTE <- #mol O2/J
    (alpha_pump * delta_O2) / J_per_kg
  
  OTE <- #kg O2/kWh
    OTE * 
    (O2_MM / 1000) * #mol O2 -> kg O2
    J_to_kWh #J -> kWh
  
  return(OTE)
  
}
