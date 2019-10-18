# Plot OTE models

#----Initialize_workspace----

#Run the models
source("run_OTE_models.R")

#Load necessary packages (many already loaded in previously sourced scripts)

library(viridis)
library(cowplot)
library(ggpubr)

#----Set_common_plot_properties----

#Colour scale
OTE_colour_scale <- 
  inferno(n = 4,
          begin = 0.1,
          end = 0.9)

#Ribbon transparency
transparency <- 0.5

#Depth cutoff

  #Percent saturation cutoff
  cutoff_sat <- 0.75

depth_cutoffs <- 
  hydrocasts_df %>% 
  mutate(x_O2 = 0.21,
         P_surface = 1,
         sat_O2 = S_O2(temperature = Temperature, salinity = Salinity) * x_O2 * P_surface, 
         percent_sat = O2 / sat_O2) %>% 
  filter(percent_sat <= cutoff_sat) %>% 
  group_by(Location) %>% 
  summarize(min_depth = min(Depth))

#Depth limits
depth_limits <-
  rbind(
    c(15,0),
    c(15,0),
    c(25,0),
    c(25,0),
    c(25,0)
  )

#Test cases
sites <- 
  OTE_model_results %>% 
  mutate(Location = as.character(Location)) %>% 
  distinct(Location) %>% 
  pull()
  
#Logarithmic scale
log_scale <- 
  scale_y_log10(name = expression(Oxygen~Transfer~Efficiency~~(kg~~O[2]~~kWh^{-1})),
                limits = c(10^-6, 10^4),
                breaks = c(10^c(-6:3)),
                labels = c(expression(10^{-6}),
                           "",
                           expression(10^{-4}),
                           "",
                           expression(10^{-2}),
                           "",
                           1,
                           10,
                           100,
                           1000))


#----Build_individual_plots----

OTE_plots <- list()

for (i in 1:length(sites)) {
  
  OTE_plots[[i]] <- 
    OTE_model_results %>% 
    filter(Location == sites[i]) %>% 
    #Plot for all depths below depth cutoff
    filter(Depth > depth_cutoffs %>% filter(Location == sites[i]) %>% pull(min_depth)) %>% 
    group_by(Depth, technique) %>% 
    #Calculate ranges for each technique
    dplyr::summarize(max_OTE = max(OTE, na.rm = TRUE),
                     min_OTE = min(OTE, na.rm = TRUE),
                     median_OTE = median(OTE, na.rm = TRUE)) %>% 
    ungroup() %>% 
    ggplot() +
    #Create ribbon plot
    geom_ribbon(aes(x = Depth,
                    ymin = min_OTE,
                    ymax = max_OTE,
                    fill = technique),
                alpha = transparency) +
    geom_path(aes(x = Depth,
                  y = median_OTE,
                  group = technique,
                  colour = technique),
              show.legend = FALSE) +
    #Flip axes to make depth profile
    coord_flip() +
    scale_x_reverse(name = "Depth (m)",
                    limits = depth_limits[i,]) +
    scale_y_continuous(name = expression(Oxygen~Transfer~Efficiency~~(kg~~O[2]~~kWh^{-1}))) +
    scale_fill_manual(name = element_blank(),
                      values = OTE_colour_scale) +
    scale_colour_manual(values = OTE_colour_scale) +
    labs(title = sites[i]) +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(size = 10),
          plot.title = element_text(face = "plain"))
  
}


#----Build_aggregate_figures----

#First grab legend
legend <- 
  get_legend(
    OTE_plots[[1]] + theme(
      legend.title.align = 0.5,
      legend.direction = "vertical",
      legend.key.width = unit(0.75, "in"),
      legend.justification = "center"
    )
  ) %>% as_ggplot()

#Combine all figure-specific sites and use common legend

all_sites_OTE_plot <-
  plot_grid(
    OTE_plots[[3]] + theme(legend.position = "none"),
    OTE_plots[[4]] + theme(legend.position = "none"),
    OTE_plots[[5]] + theme(legend.position = "none"),
    OTE_plots[[1]] + theme(legend.position = "none"),
    OTE_plots[[2]] + theme(legend.position = "none"),
    legend,
    ncol = 3,
    align = "hv",
    labels = c("A", "B", "C", "D", "E", "")
  )

#Build log-scale version of the figure
all_sites_OTE_plot_log_scale <-
  plot_grid(
    OTE_plots[[3]] + theme(legend.position = "none") + log_scale,
    OTE_plots[[4]] + theme(legend.position = "none") + log_scale,
    OTE_plots[[5]] + theme(legend.position = "none") + log_scale,
    OTE_plots[[1]] + theme(legend.position = "none") + log_scale,
    OTE_plots[[2]] + theme(legend.position = "none") + log_scale,
    legend,
    ncol = 3,
    align = "hv",
    labels = c("A", "B", "C", "D", "E", "")
  )

#----Export_aggregate_figures----

cowplot::ggsave(filename = "figures/figure_3.pdf",
                plot = all_sites_OTE_plot,
                height = 6.5,
                width = 9.5,
                units = "in")

cowplot::ggsave(filename = "figures/figure_S1.pdf",
                plot = all_sites_OTE_plot_log_scale,
                height = 6.5,
                width = 9.5,
                units = "in")