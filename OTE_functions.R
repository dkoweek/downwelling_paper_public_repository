# Functions necessary to calculate oxygen transfer efficiency for downwelling, subsurface aeration, and surface aeration
# Written by David Koweek on 6 February 2019
# Updated August 2019


#----Downwelling_functions----

volume_flow_per_length <- function(rho_0, delta_rho, h_max, delta_rho_h, eta) {
  
  g <- 9.81 #m/s^2
  
  q_flow <-
    (rho_0 / (g * delta_rho)) *
    (
      ((h_max * g * delta_rho_h) /
        (rho_0 * (eta ^ 2))) ^
       (3 / 2)
      )
  
  return(q_flow)
  
}

#Unit power (J/kg of seawater pumped)
unit_power <- function(depth, temperature, salinity, r_to_d, v, h_max_to_d, eta = 2.5, f_D = 0.06, K = 0.25, alpha_pump = 1) {
  
  g <- 9.81 #m/s^2
  
  #Pipe radius (m)
  r <- r_to_d * depth
  
  #Max plume height above the release point
  h_max <- h_max_to_d * depth
  
  #Seawater density (kg/m^3)
  rho <-
    oce::swRho(salinity = salinity,
               temperature = temperature,
               pressure = 0)
  
  rho_fun <-
    splinefun(x = depth,
              y = rho)
  
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
  
  #Delta rho between surface and bottom (kg/m^3)
  delta_rho <- 
    rho - rho_0
  
  #Delta rho between surface and maximum height of the plume (kg/m^3)
  delta_rho_h <- 
    rho_fun(depth - h_max) - 
    rho_0
  
  #Delta_rho_h is undefined when the plume rises shallower than the first depth observation
  delta_rho_h[which((depth - h_max < depth[1]) == TRUE)] <- NA
  
  #Flow per unit length (m^2/s)
  q_flow <- 
    volume_flow_per_length(rho_0 = rho_0,
                           delta_rho = delta_rho,
                           h_max = h_max,
                           delta_rho_h = delta_rho_h,
                           eta = eta)
  
  
  #Unit power (J/kg)
  #Break equation into terms to prevent onion function
  A <- 
    0.5 * v^2
  
  B <-
    (1 +
       (
         (f_D / 2) * 
           (
             ((pi * r * v) / q_flow) +
               (depth / r)
           )
       ) +
       K)
  
  C <-
    g * depth * ((rho_bar / rho_0) - 1)
  
  E_kg <- 
    (1 / alpha_pump) * ((A * B) + C)
  
  return(E_kg)
  
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
  
  delta_O2_pump <- #mol/kg O2
    O2_0 - O2
  
  return(delta_O2_pump)
  
}

#Oxygen transfer efficiency of downwelling
OTE_downwelling <- function(depth, temperature, salinity, O2, r_to_d, v, h_max_to_d, eta = 2.5, f_D = 0.06, K = 0.25, alpha_pump = 1) {
  
  #Define constant values
  g <- 9.81 #gravitational constant (m/s^2)
  
  O2_MM <- 32 #O2 molar mass (g/mol O2)
  
  J_to_kWh <- 3.6e6 #Joules per kWh
  
  joules_kg_seawater <- #Energy expenditure to pump seawater (J/kg)
    unit_power(depth = depth,
               temperature = temperature,
               salinity = salinity,
               r_to_d = r_to_d,
               v = v,
               h_max_to_d = h_max_to_d,
               eta = eta,
               f_D = f_D,
               K = K,
               alpha_pump = alpha_pump)
  
  dO2 <- #mol/kg
    delta_O2_pump(depth = depth,
                  O2 = O2)
  
  OTE <- #mol O2/J
    dO2 / joules_kg_seawater
  
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
k_rate <- function(v_f) {
  A <- 2e-4 #s/m^2
  
  k <- #s^-1
    A * v_f^2 
  
  return(k)
  
}

delta_O2_fountain <- function(depth, temperature, salinity, O2, v_f, phi) {
  
  #Define constant values
  g <- 9.81 #gravitational constant (m/s^2)
  
  phi <- #convert degrees to radians
    phi * (pi / 180)
  
  O2_prime <- #mol/kg
    O2_deficit(temperature = temperature,
               salinity = salinity,
               O2 = O2)

  k <- #s^-1
    k_rate(v_f)
  
  delta_O2 <- 
    O2_prime *
    (exp((-2 * k * v_f * sin(phi)) / g) - 1)
  
  return(delta_O2)
  
  
}

OTE_fountain <- function(depth, temperature, salinity, O2, v_f, phi = 60, alpha_pump = 1, epsilon_f = 6) {
  
  #Define constant values
  g <- 9.81 #gravitational constant (m/s^2)
  
  O2_MM <- 32 #O2 molar mass (g/mol O2)
  
  J_to_kWh <- 3.6e6 #Joules per kWh
  
  delta_O2 <- #mol O2/kg
    delta_O2_fountain(depth = depth,
                      temperature = temperature,
                      salinity = salinity,
                      O2 = O2,
                      v_f = v_f,
                      phi = phi)
  
  J_per_kg <- #Energy per kg of water
    (depth * g) + 
      (
        (v_f^2 /2) *
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
