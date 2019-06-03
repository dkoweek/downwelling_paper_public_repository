oxygenation_expt_sensor_data <- 
  read_csv("data/oxygenation_expt_sensor_data.csv",
         col_types = list(col_character(),
                          col_character(),
                          col_character(),
                          col_double(),
                          col_double(),
                          col_double(),
                          col_datetime(),
                          col_double(),
                          col_double(),
                          col_character(),
                          col_character(),
                          col_integer(),
                          col_character(),
                          col_character() 
         )
  )

oxygenation_expt_sensor_data <-
  oxygenation_expt_sensor_data %>% 
  mutate(local_datetime =  ymd_hms(local_datetime,
                                   tz = "US/Pacific"))