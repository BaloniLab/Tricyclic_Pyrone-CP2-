# monte carlo simulation for CP2 project 
library(dplyr)
library(tidyr)
library(tidyverse)
library(tidyplots)
library(ggbreak)

###### working with flux sampling results 
# read the flux sampling output 
sampling_list <- "/Users/egabal/Library/CloudStorage/Box-Box/CP2_project/Flux_sampling/"
sampling_files <- list.files(path = sampling_list, pattern = "\\.csv$", full.names = TRUE)
sampling <- setNames(
  lapply(sampling_files, read.csv),
  basename(tools::file_path_sans_ext(sampling_files))
)

# combine all flux values per sample for each reactions
cp2_samples <- grep("^iMAT_CP2-APPPS1-[1-5]_sampling$", names(sampling), value = TRUE)
common_cp2 <- Reduce(intersect, lapply(sampling[cp2_samples], colnames))
CP2 <- do.call(rbind, lapply(sampling[cp2_samples], function(df) df[, common_cp2, drop = FALSE]))
rownames(CP2) <- NULL

appps1_samples <- grep("^iMAT_Veh-APPPS1-[1-5]_sampling$", names(sampling), value = TRUE)
common_appps1 <- Reduce(intersect, lapply(sampling[appps1_samples], colnames))
APPPS1 <- do.call(rbind, lapply(sampling[appps1_samples], function(df) df[, common_appps1, drop = FALSE]))
row.names(APPPS1) <- NULL

ctrl_samples <- grep("^iMAT_Veh-NTG[1-5]_sampling$", names(sampling), value = TRUE)
common_ctrl <- Reduce(intersect, lapply(sampling[ctrl_samples], colnames))
control <- do.call(rbind, lapply(sampling[ctrl_samples], function(df) df[, common_ctrl, drop = FALSE]))
row.names(control) <- NULL

# extract the reaction of interest 
reaction_flux <- data.frame(
  APPPS1_CP2 = CP2[['ACRNtm']])



# First, pivot the data into a long format
reaction_flux_longer <- reaction_flux %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "FluxValue")

reaction_flux_longer <- reaction_flux_longer %>%
  mutate(
    FluxValue = as.numeric(unlist(FluxValue)),  # Unlist and convert FluxValue to numeric
    Sample = factor(Sample, levels = c("NTG_Veh", "APPPS1_Veh", "APPPS1_CP2"))                 # Ensure Sample is a factor
  ) %>%
  drop_na(FluxValue)  # Remove rows with NA in 'FluxValue'

colors <- new_color_scheme(c("#458B00","#E2D269", "#ECA669"))

p <- reaction_flux_longer |>      
  tidyplot(x = Sample, y = FluxValue, color = Sample) |>    
  add_violin(alpha = 0.7, color = "black", draw_quantiles = 0.3) |>
  add_boxplot(alpha = 0.7, color = "black", fill = "white") |>  # removed width
  adjust_colors(new_colors = colors) |>
  add_test_asterisks(method = "wilcox_test", p.adjust.method = "BH")

# Apply custom ggplot2 theme modifications
p <- p + 
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, color = "black"),
    axis.title.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.text.y = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 14, face = "bold")) + labs(title = "MDHm", y = "Flux value (mmol/gDCW/h)") +
  scale_x_discrete(labels = c("NTG-Veh", "APP/PS1-Veh", "APP/PS1-CP2"))



#ggsave("~/Desktop/Flux_sampling_CP2_FIGURES/cholesterol_metabolism/MDHm.png", plot = p, bg = "transparent", limitsize = F, dpi = 600)



####### HISTOGRAM

# Define your custom color palette as a named vector

reaction_flux <- data.frame(
  NTG_Veh = control[['ACRNtm']],
  APPPS1_Veh = APPPS1[['ACRNtm']],
  APPPS1_CP2 = CP2[['ACRNtm']])

custom_colors <- c(
  "NTG_Veh" = "dodgerblue",
  "APPPS1_Veh" = "grey",
  "APPPS1_CP2" = "goldenrod1"
)

# Pivot longer
reaction_flux_longer_hist <- reaction_flux %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "FluxValue") %>%
  mutate(
    FluxValue = as.numeric(FluxValue),
    Sample = factor(Sample, levels = names(custom_colors))
  ) %>%
  drop_na(FluxValue)

# Plot
reaction_hist <- ggplot(reaction_flux_longer_hist, aes(x = FluxValue, fill = Sample)) +
  geom_histogram(alpha = 0.7, position = "identity", bins = 30, color = "black") +
  geom_freqpoly(aes(color = Sample), bins = 30, linewidth = 1.2) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  labs(x = "Flux Value", y = "Frequency") +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 18, face = "bold", color = "black"),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold", color = "black"),
    axis.title.y = element_text(size = 15, face = "bold", color = "black"),
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5)
  ) +
  labs(title = "MDHm", x = "Flux value (mmol/gDCW/h)") 


ggsave("/Users/egabal/Library/CloudStorage/Box-Box/Dr_Trushina_project/selected_data_my_paper/Flux_sampling_CP2_FIGURES/energy_metabolism/MDHm.png", plot = reaction_hist, limitsize = F, dpi = 600, bg = "white")


