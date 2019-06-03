# Custom function to read in processed CastAway CTD files (.csv)

library(tidyverse)

read_castaway <- function(filename) {
  
  #Read in hydrocast metadata
  hydrocast_metadata <- 
    read_lines(filename,
               n_max = 29)
  
  #Parse timestamp (use UTC)
  timestamp <- 
    str_extract(hydrocast_metadata[3],
                "\\d+-\\d+-\\d+ \\d+:\\d+:\\d+") %>%
    ymd_hms(tz = "UTC")
  
  #Parse Longitude
  
  #Extract starting longitude
  start_longitude <- 
    str_extract(hydrocast_metadata[11],
                "\\d+.\\d+") %>%
    as.numeric()
  
  #Extract ending longitude
  end_longitude <- 
    str_extract(hydrocast_metadata[17],
                "\\d+.\\d+") %>%
    as.numeric()
  longitude <- mean(c(start_longitude, end_longitude))
  
  #Parse Latitude
  
  #Extract starting latitude
  start_latitude <- 
    str_extract(hydrocast_metadata[10],
                "\\d+.\\d+") %>%
    as.numeric ()
  
  #Extract ending latitude
  end_latitude <- 
    str_extract(hydrocast_metadata[16],
                "\\d+.\\d+") %>%
    as.numeric()
  latitude <- mean(c(start_latitude, end_latitude))
  
  
  #Read in hydrocast data
  hydrocast_data <- 
    read_csv(filename,
             comment = "%")
  
  #Output tibble with hydrocast data and metadata
  hydrocast <- 
    hydrocast_data %>%
    mutate(date_time = timestamp) %>%
    mutate(latitude = latitude) %>%
    mutate(longitude = longitude) %>%
    tbl_df()
  
  return(hydrocast)
  
}

read_ysi <- function(filename, timezone) {
  
  #Read in raw data file
  data <-
    read_csv(filename,
             comment = "\"",
             col_types = cols())
  
  #Read in data file again to extract header names    
  headers <- 
    read_csv(filename,  
             quote = "\"",
             col_types = cols())
  
    
  #Add headers to raw data
    #Creates unique header for every column
  headers[2,] <- str_replace_na(headers[2,])
  names(data) <- 
    apply(headers, 
          2, 
          function(x) str_c(x[1], 
                            x[2], 
                            sep = " "))
  
  
  #Convert date and time columns to single date-time object
  timestamp <-
    data %>% 
    select(c(starts_with("Date"), 
             starts_with("Time"))) %>% 
    unite(col = "timestring",
          sep = " ") %>% 
    pull("timestring") %>% 
    parse_date_time(orders = c("mdY HMS", "dmY HMS", "Ymd HMS"), 
                    tz = timezone) %>% 
    tbl_df()
    
  data <-
    bind_cols(timestamp, data) %>%
    dplyr::select(-c(starts_with("Date"), starts_with("Time"))) %>%
    dplyr::rename(datetime = value)
  
  
 return(data)
 
 
}

read_garmin_csv <- function(filename, skiplines) {
  
  data <- 
    read_csv(
      filename,
      skip = skiplines
    ) %>%
    tbl_df()
  
  return(data)
  
}

read_minidot <- function(filename) {
  
  #File metadata
  metadata <- 
    read_lines(filename,
               n_max = 9)
  
  #Extract serial number
  s_n <-
    str_extract(metadata[2],
                "\\d+-\\d+")
  
  #Extract salinity (entered during data download to calculate percent sat)
  salinity <- 
    str_extract(metadata[5],
                "\\d+.\\d+") %>% 
    as.numeric
  
  #Extract elevation (entered during data download to calculate percent sat)
  elevation <- 
    str_extract(metadata[6],
                "\\d+.\\d+") %>% 
    as.numeric
  
  
  #Read in data headers
  minidot_data_headers <- 
    read_csv(filename,
             skip = 7,
             col_types = cols()) %>% 
    colnames()
    
  #Read in data
  suppressWarnings(
    minidot_data_raw <-
      read_csv(filename,
               skip = 8,
               col_types = cols()) %>%
      setNames(minidot_data_headers)
  )
  
  
  #Combine metadata with MiniDOT data
  minidot_data <- 
    minidot_data_raw %>% 
    mutate(serial_number = s_n,
           salinity = salinity,
           elevation = elevation)
    
    
  return(minidot_data)  
  
}