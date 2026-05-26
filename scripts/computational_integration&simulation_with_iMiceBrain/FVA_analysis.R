# flux variability analysis (FVA) for CP2 study
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(tidyplots)
library(pheatmap)
library(VennDiagram)
library(grid)
library(futile.logger)
library(ggVennDiagram)
library(UpSetR)
library(gridExtra)
library(igraph)
library(ggraph)
library(plotly)
library(forcats)


path <- "/Users/egabal/Library/CloudStorage/Box-Box/CP2_project/FVA_output"
files_list <- list.files(path = path, pattern = ".csv", full.names = TRUE)

FVA_data <- list()

# loop through files and assign a zero for reactions with vmin and vmax zero and other reactions a value 1
for (file in files_list) {
  file_name <- tools::file_path_sans_ext(basename(file))
  df <- read.csv(file)
  # Assign 0 if both vmin and vmax are 0, otherwise 1
  df[[file_name]] <- ifelse(df[,"minimum"] == 0 & df[,"maximum"] == 0, 0, 1)
  df <- df[, !(names(df) %in% c("minimum", "maximum"))] # exclude after the pinarization
  FVA_data[[file_name]] <- df
}

FVA_combined <- reduce(FVA_data, full_join, by = c("reaction_id", "subsystem", "genes", "reaction_name"))
colnames(FVA_combined) <- gsub("_FVA", "", colnames(FVA_combined))
colnames(FVA_combined) <- gsub("iMAT_", "", colnames(FVA_combined))

FVA_combined_new <- FVA_combined[,-c(2,4)]
head(FVA_combined_new,5)

sample_cols <- names(FVA_combined_new)[-(1:2)]

# Build the group map correctly
group_map <- data.frame(
  SampleID = sample_cols,
  Group  = c(rep("APP/PS1-CP2", 5),  
             rep("APP/PS1-Veh", 5),  
             rep("NTG-Veh", 5)))

# convert to longer format
fva_long <- FVA_combined_new %>%
  pivot_longer(
    cols = -c(reaction_id, subsystem),
    names_to = "SampleID",
    values_to = "active"
  ) %>%
  filter(active == 1)

fva_long <- left_join(fva_long, group_map, by = "SampleID")

fva_summary <- fva_long %>%
  group_by(subsystem, Group) %>%
  summarise(n_active = n(), .groups = "drop")

fva_summary <- fva_summary[-c(1:3),]

# i will exclude the extracellular transport 
fva_summary <- fva_summary %>% filter(!subsystem %in% c("Transport, extracellular", "Exchange/demand reaction", "Miscellaneous"))

# get the top 10 subsystems
top_subsystems <- fva_summary %>%
  group_by(subsystem) %>%
  summarise(total_active = sum(n_active)) %>%
  arrange(desc(total_active)) %>%
  slice_head(n = 10) %>%
  pull(subsystem)

# Filter fva_summary to only keep those top 10
fva_summary_top10 <- fva_summary %>%
  filter(subsystem %in% top_subsystems)

# plot the stackplot top 10 subsystems
fva_summary_top10 <- fva_summary_top10 %>%
  mutate(subsystem = fct_reorder(subsystem, n_active))

my_colors <- new_color_scheme(c("#CD0000", "#FFD700", "#0000FF"), name = "new_color")

p <- fva_summary_top10 |>
  tidyplot(x = n_active, y = subsystem, color = Group) |>
  add_barstack_absolute() |>
  adjust_x_axis_title("Frequency of active reactions") |>
  adjust_y_axis_title("iMiceBrain Subsystem") |>
  adjust_colors(new_colors = my_colors)

ggsave("~/Desktop/our_study_FVA_active_reactions_subsystems.png", plot = p, dpi = 600, bg = "white", limitsize = F)

