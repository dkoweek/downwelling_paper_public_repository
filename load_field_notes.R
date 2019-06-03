# Load field notes from Searsville Lake oxygenation experiment - September 2018
# Written by David Koweek on 24 Oct 2018

#----Initialize_workspace----

# Load necessary packages
library(tidyverse)
library(readxl)

# Define field notes directory
field_notes_dir <- 
  "data/"

#----Load_treatment_log----
treatment_log <-
  read_xlsx(path = str_c(field_notes_dir,"treatment_log.xlsx"))

#----Load_study_site_activity_log----
study_site_activity <- 
  read_xlsx(path = str_c(field_notes_dir, "study_site_activity.xlsx"))

#----Load_hydrocast_log----
hydrocasts_log <- 
  read_xlsx(path = str_c(field_notes_dir, "hydrocasts_log.xlsx"))

#----Load_sensor_locations_log----
sensor_locations_log <-
  read_xlsx(path = str_c(field_notes_dir, "sensor_locations_log.xlsx"),
            col_types = c("text",
                          "text",
                          "text",
                          "numeric",
                          "numeric",
                          "numeric",
                          "date",
                          "date",
                          "text"))