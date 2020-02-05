# Script to build Figure 2 in "Alleviating hypoxia through induced downwelling"


#Plot hydrocasts of hypoxic data sets
#Written by David Koweek on 12 Apr 2019


#----Initiailize_workspace----

#Source data aggregation script (which sources most needed packages)
source("scrape_hypoxia_hydrography.R")

#Source OTE model functions script (contains function for O2 solubility)
source("OTE_functions.R")

#Load additional necessary packages
library(oce)
library(cowplot)
library(ggpubr)
library(viridis)

#----Build_plots_individual_variables----

#Define y-scale properties
y_limits <- c(25,0)
y_breaks = seq(max(y_limits),
               min(y_limits),
               by = -5)

#define common colour scale
palette <- 
  magma(n = 5,
        begin = 0,
        end = 0.8)

hydrocast_theme <- 
  theme(axis.text = element_text(size = 14,
                                 family = "Helvetica"),
        axis.title = element_text(size = 14,
                                  family = "Helvetica"),
        axis.text.x = element_text(angle = -45))


T_plot <- 
  hydrocasts_df %>% 
  ggplot(aes(x = Temperature,
             y = Depth)) +
  geom_path(aes(colour = Location)) +
  scale_y_continuous(name = "Depth (m)",
                     limits = y_limits,
                     breaks = y_breaks,
                     trans = "reverse") +
  scale_colour_manual(name = element_blank(),
                      values = palette) +
  scale_x_continuous(name = expression(Temperature~(degree~C))) +
  hydrocast_theme

S_plot <- 
  hydrocasts_df %>% 
  mutate(Salinity = case_when(Salinity == 0 ~ NA_real_,
                              TRUE ~ Salinity)) %>% 
  ggplot(aes(x = Salinity,
             y = Depth)) +
  geom_path(aes(colour = Location)) +
  scale_y_continuous(name = "Depth (m)",
                     limits = y_limits,
                     breaks = y_breaks,
                     trans = "reverse") +
  scale_colour_manual(name = element_blank(),
                      values = palette) +
  scale_x_continuous(name = "Salinity",
                     limits = c(10,36)) +
  hydrocast_theme

rho_plot <- 
  hydrocasts_df %>% 
  mutate(Density = swRho(salinity = Salinity,
                         temperature = Temperature,
                         pressure = Depth)) %>% 
  ggplot(aes(x = Density,
             y = Depth)) +
  geom_path(aes(colour = Location)) +
  scale_y_continuous(name = "Depth (m)",
                     limits = y_limits,
                     breaks = y_breaks,
                     trans = "reverse") +
  scale_colour_manual(name = element_blank(),
                      values = palette) +
  scale_x_continuous(name = expression(rho~(kg~m^{-3})),
                     limits = c(995, 1025),
                     breaks = seq(995, 1025, by = 10)) +
  hydrocast_theme

O2_plot <- 
  hydrocasts_df %>% 
  ggplot(aes(x = O2*1e6,
             y = Depth)) +
  geom_path(aes(colour = Location)) +
  scale_y_continuous(name = "Depth (m)",
                     limits = y_limits,
                     breaks = y_breaks,
                     trans = "reverse") +
  scale_colour_manual(name = element_blank(),
                      values = palette) +
  scale_x_continuous(name = expression(O[2]~(mu~mol~kg^{-1}))) +
  hydrocast_theme

O2_sat_plot <- 
  hydrocasts_df %>% 
  mutate(x_O2 = 0.21,
         P_surface = 1,
         sat_O2 = S_O2(temperature = Temperature, salinity = Salinity) * x_O2 * P_surface, 
         percent_sat = (O2 / sat_O2) * 100) %>% 
  ggplot(aes(x = percent_sat,
             y = Depth)) +
  geom_path(aes(colour = Location)) +
  geom_vline(xintercept = 100,
             colour = "black",
             linetype = "dotted") +
  scale_y_continuous(name = "Depth (m)",
                     limits = y_limits,
                     breaks = y_breaks,
                     trans = "reverse") +
  scale_colour_manual(name = element_blank(),
                      values = palette) +
  scale_x_continuous(name = expression(paste(O[2], " (% saturation)", sep= ""))) +
  hydrocast_theme

#----Build_map_of_hydrocasts----

locations_df <- 
  hydrocasts_df %>% 
  group_by(Location) %>% 
  summarize(Lat = mean(Latitude, na.rm = TRUE),
            Long = mean(Longitude, na.rm = TRUE))

#Plot world map
world_map <- 
  map_data("world")

map <- 
  ggplot() + 
  coord_fixed() +
  geom_polygon(data = world_map,
               aes(x = long,
                   y = lat,
                   group = group),
               fill = "grey") +
  geom_point(data = locations_df,
             aes(x = Long,
                 y = Lat,
                 colour = Location),
             size = 3) +
  scale_colour_manual(name = element_blank(),
                      values = palette,
                      labels = c(expression(atop("Angle Lake, WA, USA",
                                                 paste("(47.42750",degree,"N, -122.28710",degree,"E)"))), 
                                 expression(atop("Beaverdam Res., VA, USA",
                                                 paste("(37.32287",degree,"N, -79.82483",degree,"E)"))),
                                 expression(atop("Boknis Eck, Baltic Sea",
                                                 paste("(54.52950",degree,"N, 10.03933",degree,"E)"))),
                                 expression(atop("Chesapeake Bay, USA",
                                                 paste("(38.97567",degree,"N, -76.36850",degree,"E)"))),
                                 expression(atop("Gulf of Mexico, USA",
                                                 paste("(28.90420",degree,"N, -90.30150",degree,"E)"))))) +
  scale_x_continuous(name = expression(Longitude~(degree~East))) +
  scale_y_continuous(name = expression(Latitude~(degree~North))) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white',
                                        colour = 'white'),
        legend.key.size = unit(1.5, "cm"))




#----Combine_individual_plots----

#First grab the legend
legend <- 
  get_legend(
    O2_plot + theme(
      legend.title.align = 0.5,
      legend.direction = "vertical",
      legend.key.width = unit(0.75, "in"),
      legend.justification = "center"
    )
  ) %>% as_ggplot()


#Now combine the panels
hydrocasts_plot <-
  plot_grid(
    T_plot + theme(legend.position = "none"),
    S_plot + theme(legend.position = "none"),
    rho_plot + theme(legend.position = "none"),
    O2_plot + theme(legend.position = "none"),
    O2_sat_plot + theme(legend.position = "none"),
    legend,
    ncol = 3,
    align = "hv",
    labels = c("A", "B", "C", "D", "E", "")
  )

hydrocasts_figure <- 
  plot_grid(
    hydrocasts_plot,
    map,
    nrow = 2,
    rel_heights = c(2,1)
  )

#Export plot
cowplot::ggsave(filename = "figures/figure_2.pdf",
                plot = hydrocasts_figure,
                width = 10.5,
                height = 10,
                units = "in")