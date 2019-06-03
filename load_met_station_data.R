# Load meteorological station data
# Written by David Koweek on 20 March 2019

#Conversions
MPH_to_ms = 0.44704 


met_data <- 
  read_csv("data/meteo_station_data.txt", 
           col_types = cols(Time = col_datetime(format = "%Y-%m-%d %H:%M:%S"),
                            TemperatureF = col_double(),
                            DewpointF = col_double(),
                            PressureIn = col_double(),
                            WindDirection = col_character(),
                            WindDirectionDegrees = col_double(),
                            WindSpeedMPH = col_double(),
                            WindSpeedGustMPH = col_double(),
                            Humidity = col_double(),
                            HourlyPrecipIn = col_double(),
                            Conditions = col_character(),
                            Clouds = col_character(),
                            dailyrainin = col_double(),
                            `SolarRadiationWatts/m^2` = col_double(),
                            SoftwareType = col_double(),
                            DateUTC = col_datetime(format = "%Y-%m-%d %H:%M:%S"))) %>% 
  mutate(Time = force_tz(Time, tzone = "America/Los_Angeles"),
         Temperature = (5/9) * (TemperatureF - 32),
         Wind_Speed = WindSpeedMPH * MPH_to_ms) %>% 
  rename(local_datetime = Time,
         Solar_radiation = `SolarRadiationWatts/m^2`,
         datetime = DateUTC)
