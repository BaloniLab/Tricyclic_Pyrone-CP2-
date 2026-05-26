# flux balance analysis (FBA/FBA) CP2 metabolic simulation
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

################# looping through the files and extract the second sheet only 
path <- "/Users/eso1993/Library/CloudStorage/Box-Box/CP2_project/FBA_output/"
files_list <- list.files(path, pattern = ".xlsx", full.names = TRUE)
data_list <- list()

# loop through the files and read it then exclude the reactions with fluxes = zero which represent inactive reactions
for (file in files_list) {
  file_name <- tools::file_path_sans_ext(basename(file))
  df = readxl::read_excel(file, sheet = 2)
  df = df %>% filter(!Flux == 0)
  data_list[[file_name]] <- df
}

# export these files as csv files
for (name in names(data_list)) {
df <- data_list[[name]]  # extract the data frame
write.csv(df, file = paste0("/Users/eso1993/Library/CloudStorage/Box-Box/CP2_project/FBA_output/files_with_active_reactions/active_fluxes_", name, ".csv"), row.names = FALSE)}


################# identify the reactions which was consumed in objective function under all conditions
shared_ntg <- Reduce(intersect, data_list[c("iMAT_Veh-NTG1_FBA", 
                                            "iMAT_Veh-NTG2_FBA", 
                                            "iMAT_Veh-NTG3_FBA", 
                                            "iMAT_Veh-NTG4_FBA", 
                                            "iMAT_Veh-NTG5_FBA")])

shared_ntg <- shared_ntg %>% filter(Subsystem != "Exchange/demand reaction")

shared_apps1 <- Reduce(intersect, data_list[c("iMAT_Veh-APPPS1-1_FBA",
                                              "iMAT_Veh-APPPS1-2_FBA",
                                              "iMAT_Veh-APPPS1-3_FBA",
                                              "iMAT_Veh-APPPS1-4_FBA",
                                              "iMAT_Veh-APPPS1-5_FBA")])
shared_apps1 <- shared_apps1 %>% filter(Subsystem != "Exchange/demand reaction")


shared_cp2 <- Reduce(intersect, data_list[c("iMAT_CP2-APPPS1-1_FBA",
                                              "iMAT_CP2-APPPS1-2_FBA",
                                              "iMAT_CP2-APPPS1-3_FBA",
                                              "iMAT_CP2-APPPS1-4_FBA",
                                              "iMAT_CP2-APPPS1-5_FBA")])
shared_cp2 <- shared_cp2 %>% filter(Subsystem != "Exchange/demand reaction")


rxns_list <- list(
  `NTG-Veh` = trimws(shared_ntg$`Reaction ID`),
  `APP/PS1-Veh` = trimws(shared_apps1$`Reaction ID`),
  `APP/PS1-CP2` = trimws(shared_cp2$`Reaction ID`)
)

rxns_list <- lapply(rxns_list, function(x) x[!is.na(x) & x != ""]) # remove any NA

venn_rxn <- venn.diagram(
  x = rxns_list,
  category.names = c("NTG-Veh", "APP/PS1-Veh", "APP/PS1-CP2"),
  filename = "~/Desktop/Venn_FBA_CP2.png",
  
  # Venn Diagram settings
  imagetype = "png",
  height = 3000,
  width = 3000,
  resolution = 600,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#31B57BFF", "#4DACD6", "#fb8500"),
  
  # Numbers
  cex = 2,
  fontface = "bold",
  fontfamily = "sans",
  
  # Labels
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.03, 0.01, 0.03),
  cat.fontfamily = "sans",
  cat.col = "navy")




