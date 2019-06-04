# Script to build Figure S3 in "Alleviating hypoxia through induced downwelling"

#----Initialize_workspace----

#Load packages
library(tidyverse)
library(lubridate)
library(oce)
library(signal)
library(magrittr)
library(cowplot)
library(viridis)

#Set system clock to avoid Mac OS X / lubridate problems with timestamps
Sys.setenv(TZ = "America/Los_Angeles")

#----Load_ADP_data----

#Note: 'nf' stands for "near field" (i.e. in the vicinity of the field experiment)
adp_nf <-
  read.aquadoppProfiler(file = "data/JRBPNF02.PRF",
                        tz = "UTC",
                        orientation = "upward")

#----ADP_diagnostics----

#First, check the temperature readings
adp_temp_plot_nf <-
  plot(adp_nf,
       which="temperature")

# Ok temperatures seem envrionmentally realistic. Let's check instrument heading/pitch/roll
adp_position_nf <-
  plot(adp_nf,
       which=16:18)

#Instrument appears to remain stable and flat facing upward based on the pitch and roll time series...
#....for both instruments. It seems that there may have been a small shift in the NF pitch and roll...
#...~ 1800 on 24 September 2018. This period of time corresponds with the deployment of the 5-meter radius sensors.
# We had temporarily deployed a sensor at 5m, 180 degrees befor realizing that the sensor was near the ADP. 
# It appears that we may have actually hit the ADP during the initial deployment. Monitor processed data for discontinuities around this time period.

#Depth records show slight decline due to decline in lake level

adp_depth_plot_nf <-
  plot(adp_nf,
       which = "pressure")

#----Rotate_ADP_data_from_XYZ_to_ENU----
adp_nf_ENU <- xyzToEnuAdp(adp_nf)

#----Filter_ADP_data----

#Set up filtering coefficients
butter_coefs <- butter(2,0.1)

a <- butter_coefs$a
b <- butter_coefs$b

#Pre-allocate array for filtered data
v_filt_nf <- array(dim = dim(adp_nf_ENU[["v"]]))

#Establish depth bins
depth_bins_nf <- dim(adp_nf_ENU[["v"]])[2]

for (i in 1:3) { #east, north, and up
  for (j in 1:depth_bins_nf) { #for each depth bin
    
    velocity_vector <- adp_nf_ENU[["v"]][,j,i]
    
    v_filt_nf[,j,i] <- oce.filter(velocity_vector, 
                                  a, 
                                  b, 
                                  zero.phase = TRUE)
    
  }
}

#Copy data object
adp_nf_ENU_filtered <- 
  adp_nf_ENU

#Append filtered array
adp_nf_ENU_filtered[["v"]] <- 
  v_filt_nf

#----Diagnostic_plot_of_filtered_data----

par(mfrow=c(3,1))

#East
plot(adp_nf_ENU_filtered,
     which = 1)
lines(adp_nf_ENU_filtered[["time"]],
      adp_nf_ENU_filtered[["pressure"]])

#North
plot(adp_nf_ENU_filtered,
     which = 2)
lines(adp_nf_ENU_filtered[["time"]],
      adp_nf_ENU_filtered[["pressure"]])

#Up
plot(adp_nf_ENU_filtered,
     which = 3)
lines(adp_nf_ENU_filtered[["time"]],
      adp_nf_ENU_filtered[["pressure"]])

#----Wrangle_to_data_frames----

#Trim the data to between 0.75 mab and 2.5 mab
depth_lb <- 0.75
depth_ub <- 2.5

adp_nf_e <- adp_nf_ENU_filtered[["v"]][,,1] 
adp_nf_n <- adp_nf_ENU_filtered[["v"]][,,2]
adp_nf_up <- adp_nf_ENU_filtered[["v"]][,,3]

adp_nf_e <- as_tibble(adp_nf_e)  
adp_nf_n <- as_tibble(adp_nf_n)
adp_nf_up <- as_tibble(adp_nf_up)

names(adp_nf_e) <- adp_nf_ENU_filtered[["distance"]]
names(adp_nf_n) <- adp_nf_ENU_filtered[["distance"]]
names(adp_nf_up) <- adp_nf_ENU_filtered[["distance"]]

time_nf <- adp_nf_ENU_filtered[["time"]]
filtered_ssh_nf <- oce::oce.filter(adp_nf_ENU_filtered[["pressure"]], 
                                   a, 
                                   b, 
                                   zero.phase = TRUE)

adp_nf_e <- 
  adp_nf_e %>% 
  gather(depth,
         velocity) %>% 
  group_by(depth) %>% 
  mutate(timestamp = time_nf) %>% 
  ungroup() %>% 
  mutate(depth = as.numeric(depth)) %>% 
  mutate(direction = "East")

adp_nf_n <- 
  adp_nf_n %>% 
  gather(depth,
         velocity) %>% 
  group_by(depth) %>% 
  mutate(timestamp = time_nf) %>% 
  ungroup() %>% 
  mutate(depth = as.numeric(depth)) %>% 
  mutate(direction = "North")

adp_nf_up <- 
  adp_nf_up %>% 
  gather(depth,
         velocity) %>% 
  group_by(depth) %>% 
  mutate(timestamp = time_nf) %>% 
  ungroup() %>% 
  mutate(depth = as.numeric(depth)) %>% 
  mutate(direction = "Vertical")

adp_nf_all <- 
  bind_rows(adp_nf_e,
            adp_nf_n,
            adp_nf_up) %>%
  group_by(timestamp) %>% 
  nest() %>% 
  mutate(ssh = filtered_ssh_nf) %>% 
  unnest()

adp_nf <-
  adp_nf_all %>%
  dplyr::filter(depth >= depth_lb, #depth cut-offs
                depth <= depth_ub,
                timestamp %within% interval(ymd_hms("2018-09-02 02:30:00"),
                                            ymd_hms("2018-10-04 20:30:00")))

#----Generate_the_ADP_plot----

adp_nf_plot <-
  adp_nf %>%
  #Time of start of experiment and end of experiment
  dplyr::filter(timestamp > ymd_hms("2018-09-15 00:00:00"),
                timestamp < ymd_hms("2018-10-04 16:30:00")) %>%
  dplyr::mutate(local_datetime = with_tz(timestamp,
                                         tzone = "America/Los_Angeles")) %>% 
  ggplot(aes(x = local_datetime,
             y = depth)) +
  geom_raster(aes(fill = velocity)) +
  scale_fill_gradient2(name = expression(paste("u,v,w",~(m~~s^{-1}))),
                       limits = c(-0.05, 0.05),
                       midpoint = 0,
                       mid = "white",
                       high = "red",
                       low = "blue") +
  scale_x_datetime(name = element_blank(),
                   date_breaks = "3 days",
                   date_labels = "%b %d") +
  labs(y = "Meters Above Bottom",
       title = "Acoustic Doppler Profiler Velocity Measurements")  +
  geom_line(aes(y = ssh),
            colour = "black") +
  theme(axis.text.x = element_text(angle = -45)) +
  facet_grid(direction~.)

#Export figure
ggplot2::ggsave(filename = "figures/figure_S3.pdf",
                plot = adp_nf_plot,
                device = "pdf",
                height = 6,
                width = 8,
                units = "in")

#----Unload_signal_package----

#'signal' package causes conflicts with dplyr::filter, so unload package after completing figure
detach("package:signal", unload = TRUE)