# Written by David Koweek to scrape hydrographic data from hypoxic zones
# 4 February 2019

#----Initialize_workspace----

#Load necessary packages
library(tidyverse)
library(lubridate)
library(pangaear)
library(oce)

#Set system clock to avoid Mac OS X / lubridate problems with timestamps
Sys.setenv(TZ = "America/Los_Angeles")

#Clear the cache of any Pangaea data sets
pg_cache_clear(prompt = FALSE)

#----Load_data_sets----
baltic_sea_ds <- 
  pg_data(doi = "10.1594/PANGAEA.871890")

washington_lakes_ds <- 
  pg_data(doi = "10.1594/PANGAEA.884326")


# Data from Chesapeake Bay, Gulf of Mexico, and Virginia Reservoirs not easily accessible via API...
#... so downloaded onto machine and loaded from machine

chesapeake_bay_hydrocasts <- 
  read_csv(file = "data/Cast_CTD_8.csv")

virginia_reservoirs <- 
  read_csv(file = "data/CTD_Meta_13_18_final.csv",
           col_types =c(
             col_character(),
             col_character(),
             col_character(),
             col_double(),
             col_double(),
             col_double(),
             col_double(),
             col_double(),
             col_double(),
             col_double(),
             col_double(),
             col_double(),
             col_double(),
             col_double(),
             col_double(),
             col_double(),
             col_double(),
             col_double(),
             col_double(),
             col_double(),
             col_double(),
             col_double(),
             col_double(),
             col_double()
           )
  ) %>% 
  tbl_df()

nGOM_dead_zone_ds <- 
  read_delim(
    "data/PE1702_CTD.flat99.txt",
    " ",
    escape_double = FALSE,
    col_types = cols(
      cast = col_character(),
      day = col_double(),
      mon = col_double()
    ),
    trim_ws = TRUE
  )

#----Wrangle_data----

O2_MM <- 32 #molar mass of O2 (g)


# Baltic Sea
baltic_sea_latitude <- 54.529500 #Lat/Long info extracted from data web page
baltic_sea_longitude <- 10.039330

baltic_sea_df <- 
  baltic_sea_ds[[1]]$data %>% 
  mutate(datetime = parse_date_time(`Date/Time`, "%Y-%m-%d %H:%M:%S"),
         cast_date = date(datetime), #Extract date of each cast as a cast ID
         Latitude = baltic_sea_latitude,
         Longitude = baltic_sea_longitude,
         Location = "Boknis Eck, Baltic Sea") %>% 
  rename(Depth = `Depth water [m]`,
         Temperature = `Temp [°C]`,
         Salinity = Sal,
         O2 = `O2 [µmol/l]`) %>% 
  mutate(rho = swRho(salinity = Salinity,
                     temperature = Temperature,
                     pressure = 0),
         O2 = O2 / (rho / 1000) *1e-6) #umol/L -> mol/kg


# Lakes in Washington State
washington_lakes_df <- 
  washington_lakes_ds[[1]]$data %>% 
  mutate(cast_date = parse_date(`Date/Time`, "%Y-%m-%d")) %>% 
  rename(Depth = `Depth water [m]`,
         Temperature = `Temp [°C]`,
         O2 = `O2 [µmol/l]`) %>% 
  mutate(Salinity = 0, #No salinity given, assume 0 in a lake
         rho = swRho(salinity = Salinity, 
                     temperature = Temperature,
                     pressure = 0),
         O2 = O2 / (rho / 1000) * 1e-6) #umol/L -> mol/kg

# Chesapeake Bay
chesapeake_bay_df <- 
  chesapeake_bay_hydrocasts %>% 
  mutate(datetime = parse_date_time(Local_Date, "%Y%m%d"),
         cast_date = date(datetime),
         Location = "Chesapeake Bay, USA") %>% 
  rename(Depth = depth,
         Temperature = temp,
         Salinity = salinity,
         O2 = CTD_O2) %>% 
  mutate(rho = swRho(salinity = Salinity,
                     temperature = Temperature,
                     pressure = 0),
         O2 = O2 / (rho / 1000) * 1e-6) %>% #umol/L -> mol/kg
  #Depth-average the two bottles at each depth
  mutate(Depth_fac = as.factor(round(Depth, digits = 0))) %>% 
  group_by(Depth_fac) %>% 
  summarize(Depth = mean(Depth),
            Location = Location[1],
            Latitude = mean(Latitude),
            Longitude = mean(Longitude),
            cast_date = mean(cast_date),
            O2 = mean(O2),
            Temperature = mean(Temperature),
            Salinity = mean(Salinity)) %>% 
  ungroup()

# Virginia
virginia_reservoirs_df <- 
  virginia_reservoirs %>%
  mutate(cast_date = date(Date)) %>% 
  rename(Depth = Depth_m,
         Temperature = Temp_C,
         O2 = DO_mgL) %>% 
  mutate(Salinity = 0, #No salinity given, assume 0 in a lake
         rho = swRho(salinity = Salinity, 
                     temperature = Temperature,
                     pressure = 0),
         O2 = (O2 * 1e-3) / (O2_MM * (rho / 1000))) #mg/L -> mol/kg

# nGOM dead zone
nGOM_dead_zone_df <- 
  nGOM_dead_zone_ds %>% 
  filter(cast == 8) %>% 
  rename(Latitude = lat,
         Longitude = lon,
         Temperature = temp,
         Salinity = sal,
         O2 = O2_umol_kg,
         Depth = depth) %>% 
  mutate(O2 = O2 / 1e6, #umol/kg -> mol/kg
         Location = "Gulf of Mexico",
         cast_date = parse_date_time(ISO_DateTime_UTC,
                                     "%Y-%m-%d %H:%M:%S") %>% 
           date()) %>% 
  #Average upcast and downcast
  group_by(Depth) %>% 
  summarize(Location = Location[1],
            Latitude = Latitude[1],
            Longitude = Longitude[1],
            cast_date = cast_date[1],
            Temperature = mean(Temperature, na.rm = TRUE),
            Salinity = mean(Salinity, na.rm = TRUE),
            O2 = mean(O2, na.rm = TRUE)) %>% 
  ungroup()



#----Select_hydrocasts_and_merge_data_sets----

depth_interp_interval <- 0.5 #meters

#Baltic Sea
baltic_sea <- 
  baltic_sea_df %>% 
  filter(cast_date == "2013-09-05") %>% 
  select(c(Location,
           Latitude, 
           Longitude, 
           cast_date, 
           Depth, 
           Temperature, 
           Salinity, 
           O2)
  )

#Depth interpolate the Baltic Sea data set 
depth_interp <- seq(min(baltic_sea$Depth), 
                    max(baltic_sea$Depth), 
                    by = depth_interp_interval)

O2_interp <- 
  approx(x = baltic_sea$Depth,
         y = baltic_sea$O2,
         xout = depth_interp)$y

Temperature_interp <- 
  approx(x = baltic_sea$Depth,
         y = baltic_sea$Temperature,
         xout = depth_interp)$y

Salinity_interp <- 
  approx(x = baltic_sea$Depth,
         y = baltic_sea$Salinity,
         xout = depth_interp)$y

baltic_sea <- 
  data.frame(
    Depth = depth_interp,
    O2 = O2_interp,
    Temperature = Temperature_interp,
    Salinity = Salinity_interp
  ) %>% 
  mutate(Location = baltic_sea$Location[1],
         Latitude = baltic_sea$Latitude[1],
         Longitude = baltic_sea$Longitude[1],
         cast_date = baltic_sea$cast_date[1]) %>% 
  select(c(Location,
           Latitude, 
           Longitude, 
           cast_date, 
           Depth, 
           Temperature, 
           Salinity, 
           O2)
  ) 


#Angle Lake, Washington, USA
angle_lake <- 
  washington_lakes_df %>% 
  filter(Event == "Angle_Lake",
         cast_date == "2016-07-12 00:00:00") %>%
  mutate(Location = "Angle Lake, WA, USA") %>% 
  select(c(Location,
           Latitude, 
           Longitude, 
           cast_date, 
           Depth, 
           Temperature, 
           Salinity, 
           O2)
  )

#Depth interpolate the Angle Lake data set 
depth_interp <- seq(min(angle_lake$Depth), 
                    max(angle_lake$Depth), 
                    by = depth_interp_interval)

O2_interp <- 
  approx(x = angle_lake$Depth,
         y = angle_lake$O2,
         xout = depth_interp)$y

Temperature_interp <- 
  approx(x = angle_lake$Depth,
         y = angle_lake$Temperature,
         xout = depth_interp)$y

Salinity_interp <- 
  approx(x = angle_lake$Depth,
         y = angle_lake$Salinity,
         xout = depth_interp)$y

angle_lake <- 
  data.frame(
    Depth = depth_interp,
    O2 = O2_interp,
    Temperature = Temperature_interp,
    Salinity = Salinity_interp
  ) %>% 
  mutate(Location = angle_lake$Location[1],
         Latitude = angle_lake$Latitude[1],
         Longitude = angle_lake$Longitude[1],
         cast_date = angle_lake$cast_date[1]) %>% 
  select(c(Location,
           Latitude, 
           Longitude, 
           cast_date, 
           Depth, 
           Temperature, 
           Salinity, 
           O2)
  ) 


#Chesapeake Bay, USA
chesapeake_bay <- 
  chesapeake_bay_df %>% 
  select(c(Location,
           Latitude, 
           Longitude, 
           cast_date, 
           Depth, 
           Temperature, 
           Salinity, 
           O2)
  )

#Depth interpolate the Chesapeake Bay data set 
depth_interp <- seq(min(chesapeake_bay$Depth), 
                    max(chesapeake_bay$Depth), 
                    by = depth_interp_interval)

O2_interp <- 
  approx(x = chesapeake_bay$Depth,
         y = chesapeake_bay$O2,
         xout = depth_interp)$y

Temperature_interp <- 
  approx(x = chesapeake_bay$Depth,
         y = chesapeake_bay$Temperature,
         xout = depth_interp)$y

Salinity_interp <- 
  approx(x = chesapeake_bay$Depth,
         y = chesapeake_bay$Salinity,
         xout = depth_interp)$y

chesapeake_bay <- 
  data.frame(
    Depth = depth_interp,
    O2 = O2_interp,
    Temperature = Temperature_interp,
    Salinity = Salinity_interp
  ) %>% 
  mutate(Location = chesapeake_bay$Location[1],
         Latitude = chesapeake_bay$Latitude[1],
         Longitude = chesapeake_bay$Longitude[1],
         cast_date = chesapeake_bay$cast_date[1]) %>% 
  select(c(Location,
           Latitude, 
           Longitude, 
           cast_date, 
           Depth, 
           Temperature, 
           Salinity, 
           O2)
  ) 


#Beaverdam Reservoir, Virginia, USA
beaverdam_reservoir_latitude <- 37.322865
beaverdam_reservoir_longitude <- -79.824834

beaverdam_reservoir <- 
  virginia_reservoirs_df %>% 
  filter(Reservoir == "BVR",
         cast_date == "2014-06-25 00:00:00") %>% #Manually selected a summer hydrocast from Beaverdam Reservoir that shows clear hypoxia
  mutate(Location = "Beaverdam Reservoir, VA, USA",
         Latitude = beaverdam_reservoir_latitude,
         Longitude = beaverdam_reservoir_longitude ) %>% 
  select(c(Location,
           Latitude, 
           Longitude, 
           cast_date, 
           Depth, 
           Temperature, 
           Salinity, 
           O2)
  ) 

#Depth interpolate the Beaverdam Reservoir data set 
depth_interp <- seq(min(beaverdam_reservoir$Depth), 
                    max(beaverdam_reservoir$Depth), 
                    by = depth_interp_interval)

O2_interp <- 
  approx(x = beaverdam_reservoir$Depth,
         y = beaverdam_reservoir$O2,
         xout = depth_interp)$y

Temperature_interp <- 
  approx(x = beaverdam_reservoir$Depth,
         y = beaverdam_reservoir$Temperature,
         xout = depth_interp)$y

Salinity_interp <- 
  approx(x = beaverdam_reservoir$Depth,
         y = beaverdam_reservoir$Salinity,
         xout = depth_interp)$y

beaverdam_reservoir <- 
  data.frame(
    Depth = depth_interp,
    O2 = O2_interp,
    Temperature = Temperature_interp,
    Salinity = Salinity_interp
  ) %>% 
  mutate(Location = beaverdam_reservoir$Location[1],
         Latitude = beaverdam_reservoir$Latitude[1],
         Longitude = beaverdam_reservoir$Longitude[1],
         cast_date = beaverdam_reservoir$cast_date[1]) %>% 
  select(c(Location,
           Latitude, 
           Longitude, 
           cast_date, 
           Depth, 
           Temperature, 
           Salinity, 
           O2)
  ) 

#Northern Gulf of Mexico
nGOM_dead_zone <- 
  nGOM_dead_zone_df 

#Depth interpolate the Baltic Sea data set 
depth_interp <- seq(min(nGOM_dead_zone$Depth), 
                    max(nGOM_dead_zone$Depth), 
                    by = depth_interp_interval)

O2_interp <- 
  approx(x = nGOM_dead_zone$Depth,
         y = nGOM_dead_zone$O2,
         xout = depth_interp)$y

Temperature_interp <- 
  approx(x = nGOM_dead_zone$Depth,
         y = nGOM_dead_zone$Temperature,
         xout = depth_interp)$y

Salinity_interp <- 
  approx(x = nGOM_dead_zone$Depth,
         y = nGOM_dead_zone$Salinity,
         xout = depth_interp)$y

nGOM_dead_zone <- 
  data.frame(
    Depth = depth_interp,
    O2 = O2_interp,
    Temperature = Temperature_interp,
    Salinity = Salinity_interp
  ) %>% 
  mutate(Location = nGOM_dead_zone$Location[1],
         Latitude = nGOM_dead_zone$Latitude[1],
         Longitude = nGOM_dead_zone$Longitude[1],
         cast_date = nGOM_dead_zone$cast_date[1]) %>% 
  select(c(Location,
           Latitude, 
           Longitude, 
           cast_date, 
           Depth, 
           Temperature, 
           Salinity, 
           O2)
  ) 



hydrocasts_df <- 
  bind_rows(
    baltic_sea,
    angle_lake,
    chesapeake_bay,
    beaverdam_reservoir,
    nGOM_dead_zone
  ) %>% 
  group_by(Location) %>% 
  arrange(Depth,
          .by_group = TRUE) %>% 
  ungroup()