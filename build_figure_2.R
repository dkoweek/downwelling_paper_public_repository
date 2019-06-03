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
library(gridExtra)
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
  magma(n = 4,
        begin = 0,
        end = 0.8)

hydrocast_theme <- 
  theme(axis.text = element_text(size = 10,
                                 family = "Helvetica"),
        axis.title = element_text(size = 12,
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
                     limits = c(10,25)) +
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
                     limits = c(995, 1020),
                     breaks = seq(995, 1020, by = 5)) +
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

#----Combine_individual_plots----

top_row <-
  plot_grid(
  T_plot + theme(legend.position = "none"),
  S_plot + theme(legend.position = "none"),
  rho_plot + theme(legend.position = "none"),
  ncol = 3,
  align = "hv",
  labels = "AUTO"
  )

bottom_row <- 
  plot_grid(
    O2_plot + theme(legend.position = "none"),
    O2_sat_plot + theme(legend.position = "none"),
    ncol = 2,
    align = "hv",
    labels = c("D", "E")
  ) 

hydrocasts_plot <- 
  arrangeGrob(
    top_row,
    bottom_row,
    layout_matrix = rbind(c(1,1,1,1,1,1),
                          c(NA, 2, 2, 2, 2, NA))
  ) %>%
  annotate_figure(bottom = get_legend(
  O2_plot + theme(
  legend.title.align = 0.5,
  legend.direction = "vertical",
  legend.key.width = unit(0.75, "in"),
  legend.justification = "center"
  )
  ))


#Export plot
cowplot::ggsave(filename = "figures/figure_2.pdf",
                plot = hydrocasts_plot,
                width = 8,
                height = 6,
                units = "in")