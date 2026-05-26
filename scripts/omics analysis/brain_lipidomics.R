# lipidomics analysis for CP2 project
#library(BiocManager)#BiocManager::install("lipidr")
#BiocManager::install(c("ComplexHeatmap", "circlize"))

library(lipidr)
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(limma)
library(EnhancedVolcano)
library(tidyplots)
library(ComplexHeatmap)
library(circlize)
library(readr)

#################### Read data and prepare for lipidr  
lipid_data <- read.csv("/Users/egabal/Library/CloudStorage/Box-Box/CP2_project/brain_lipidomics/raw_lipidomics.csv")

# process to start with lipidr
lipid_object <- as_lipidomics_experiment(lipid_data) 

# add the meta data 
lipid_object <- add_sample_annotation(lipid_object, "/Users/egabal/Library/CloudStorage/Box-Box/CP2_project/brain_lipidomics/lipidomics_meta.csv")
colData(lipid_object) # inspect

################### Raw data quality check
lipid_object$Group <- factor(lipid_object$Group, levels = c("NTG-Veh", "APPPS1-Veh", "APPPS1-CP2"))
quality_plot <- plot_samples(lipid_object, type = "boxplot", log = TRUE, color = "Group")

quality_plot <- quality_plot + theme_minimal() +
  theme(
  axis.text.x = element_text(size = 16, face = "bold", angle = 360),
  axis.title.y = element_text(size = 17),
  axis.title.x = element_text(size = 17),
  axis.text.y = element_text(size = 14)) + labs(title = "Quality check")

#ggsave("~/Desktop/quality_check_lipidomics_cp2.png", plot = quality_plot,bg = "white", dpi = 600, limitsize = F)


################### lipid class distribution
lipid_Clsss_plot <- plot_lipidclass(lipid_object, "boxplot")
lipid_Clsss_plot <- lipid_Clsss_plot +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, size = 15, color = "black", margin = margin(t = 6, r = 5)),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 17, face = "bold"),
        axis.title.y = element_text(size = 17))
#ggsave("~/Desktop/abundance_lipidomics_cp2.png", plot = lipid_Clsss_plot,bg = "white", dpi = 600, limitsize = F)

####################### Normalization 
# Probabilistic Quotient Normalization (PQN)
# data already have single value per lipid molecule in the dataset (already summarized)

d_normalized = normalize_pqn(lipid_object, measure = "Area", exclude = "blank", log = TRUE)
norm_plot <- plot_samples(d_normalized, "boxplot", color = "Group") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 360, size = 17, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))

#ggsave("~/Desktop/normalized_lipidomics_cp2.png", plot = norm_plot,bg = "white", dpi = 600, limitsize = F)

######################## multivariate analysis PCA (unsupervised)
colData(d_normalized)$Group <- factor(colData(d_normalized)$Group, levels = c("NTG-Veh", "APPPS1-Veh", "APPPS1-CP2"))
mvaresults = mva(d_normalized, measure="Area", method="PCA")
pca <- plot_mva(mvaresults, color_by="Group", ellipse = F)
pca <- pca + geom_point(size = 9) + theme_minimal() +
  scale_color_manual(values = c(
    "APPPS1-CP2" = "#D55E00",
    "NTG-Veh"    = "#009E73",
    "APPPS1-Veh" = "#0072B2")) +
  theme(axis.text.x = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 30, color = "black"),
        axis.title.x = element_text(size = 35, color = "navy", face = "bold"),
        axis.title.y = element_text(size = 35, color = "navy", face = "bold"))

#ggsave("~/Desktop/PCA_lipidomics_cp2.png", plot = pca,bg = "white", dpi = 600, limitsize = F)

######################## multivariate analysis oPLSDA (supervised)
# 1- APP/PS1-Veh vs NTG-Veh 
mvaresults_1 <- mva(d_normalized,method = "OPLS-DA", group_col = "Group", groups = c("APPPS1-Veh", "NTG-Veh"))
OPLS_1 <- plot_mva(mvaresults_1, color_by = "Group", ellipse = F)
OPLS_1 <- OPLS_1 + geom_point(size = 6.5) + 
  stat_ellipse(aes(color = Group), linewidth = 1.2) +
  scale_color_manual(values = c(
    "NTG-Veh" = "#009E73",
    "APPPS1-Veh"    = "#0072B2")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 30, color = "black"),
        axis.title.x = element_text(size = 35, color = "navy", face = "bold"),
        axis.title.y = element_text(size = 35, color = "navy", face = "bold"))

#ggsave("~/Desktop/OPLS-DA-APPPS1_vs_NTG_lipidomics_cp2.png", plot = OPLS_1, bg = "white", dpi = 600, limitsize = F)

# loadings or lipidomics contributed to that varaiance
loading1 <- plot_mva_loadings(mvaresults_1, color_by = "Class", top.n = 40)
loading1 <- loading1 + theme_minimal() +
  theme(axis.text.x = element_text(size = 22, color = "black"),
        axis.text.y = element_text(size = 22, color = "black"),
        axis.title.x = element_text(size = 25, color = "navy", face = "bold"),
        axis.title.y = element_text(size = 25, color = "navy", face = "bold"))

# change the size for labels
idx <- which(sapply(loading1$layers, function(x)
  inherits(x$geom, "GeomLabelRepel")))
loading1$layers[[idx]]$aes_params$size <- 5.3  

#ggsave("~/Desktop/top40_lipids_OPLS-DA-APPPS1_vs_NTG_lipidomics_cp2.png", plot = loading1, bg = "white", dpi = 600, limitsize = F)

# get the information about top 40 lipid species 
apps1_top40 <- top_lipids(mvaresults_1, top.n = 40)
#write.csv(apps1_top40, "~/Desktop/top40_lipids_APPS1_vs_NTG.csv", row.names = F)

# 2- APP/PS1-CP2 vs APP/PS1-Veh 
mvaresults_2 <- mva(d_normalized,method = "OPLS-DA", group_col = "Group", groups = c("APPPS1-CP2", "APPPS1-Veh"))
OPLS_2 <- plot_mva(mvaresults_2, color_by = "Group", ellipse = F)
OPLS_2 <- OPLS_2 + geom_point(size = 6.5) + 
  stat_ellipse(aes(color = Group), linewidth = 1.2) +
  scale_color_manual(values = c(
    "APPPS1-CP2" = "#D55E00",
    "APPPS1-Veh"    = "#0072B2")) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 30, color = "black"),
        axis.title.x = element_text(size = 35, color = "navy", face = "bold"),
        axis.title.y = element_text(size = 35, color = "navy", face = "bold"))

#ggsave("~/Desktop/OPLS-DA-CP2_vs_APPPS1_lipidomics_cp2.png", plot = OPLS_2, bg = "white", dpi = 600, limitsize = F)

# loadings or lipidomics contributed to that varaiance
loading2 <- plot_mva_loadings(mvaresults_2, color_by = "Class", top.n = 40)
loading2 <- loading2 + theme_minimal() +
  theme(axis.text.x = element_text(size = 22, color = "black"),
        axis.text.y = element_text(size = 22, color = "black"),
        axis.title.x = element_text(size = 25, color = "navy", face = "bold"),
        axis.title.y = element_text(size = 25, color = "navy", face = "bold"))

# change the size for labels
idx <- which(sapply(loading2$layers, function(x)
  inherits(x$geom, "GeomLabelRepel")))
loading2$layers[[idx]]$aes_params$size <- 5.3  

#ggsave("~/Desktop/top40_lipids_OPLS-DA-cp2_vs_APPPS1_lipidomics_cp2.png", plot = loading2, bg = "white", dpi = 600, limitsize = F)

# get the information about top 40 lipid species 
cp2_top40 <- top_lipids(mvaresults_2, top.n = 40)
#write.csv(cp2_top40, "~/Desktop/top40_lipids_CP2_vs_APPPS1.csv", row.names = F)

# 3- APP/PS1-CP2 vs NTG-Veh 
mvaresults_3 <- mva(d_normalized,method = "OPLS-DA", group_col = "Group", groups = c("APPPS1-CP2", "NTG-Veh"))
OPLS_3 <- plot_mva(mvaresults_3, color_by = "Group", ellipse = F)
OPLS_3 <- OPLS_3 + geom_point(size = 6.5) + 
  stat_ellipse(aes(color = Group), linewidth = 1.2) +
  scale_color_manual(values = c(
    "APPPS1-CP2" = "#D55E00",
    "NTG-Veh"    = "#009E73")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 30, color = "black"),
        axis.title.x = element_text(size = 35, color = "navy", face = "bold"),
        axis.title.y = element_text(size = 35, color = "navy", face = "bold"))

#ggsave("~/Desktop/OPLS-DA-CP2_vs_NTG_lipidomics_cp2.png", plot = OPLS_3, bg = "white", dpi = 600, limitsize = F)

# loadings or lipidomics contributed to that varaiance
loading3 <- plot_mva_loadings(mvaresults_3, color_by = "Class", top.n = 40)
loading3 <- loading3 + theme_minimal() +
  theme(axis.text.x = element_text(size = 22, color = "black"),
        axis.text.y = element_text(size = 22, color = "black"),
        axis.title.x = element_text(size = 25, color = "navy", face = "bold"),
        axis.title.y = element_text(size = 25, color = "navy", face = "bold"))

# change the size for labels
idx <- which(sapply(loading3$layers, function(x)
  inherits(x$geom, "GeomLabelRepel")))
loading3$layers[[idx]]$aes_params$size <- 5.3  

#ggsave("~/Desktop/top40_lipids_OPLS-DA-cp2_vs_NTG_lipidomics_cp2.png", plot = loading3, bg = "white", dpi = 600, limitsize = F)

# get the information about top 40 lipid species 
cp2_2_top40 <- top_lipids(mvaresults_3, top.n = 40)
#write.csv(cp2_2_top40, "~/Desktop/top40_lipids_CP2_vs_NTG.csv", row.names = F)

######################## Differential lipids 
table(colData(d_normalized)$Group)
# replace - by _ in the group names 
colData(d_normalized)$Group <- gsub("-", "_", colData(d_normalized)$Group)

# run the differential analysis analysis 
# 1- comparison APPPS1-Veh vs NTG-Veh
de_results_1 = de_analysis(data= d_normalized, APPPS1_Veh - NTG_Veh, group_col = "Group",
                         measure = "Area")

# cut off p-value 0.05, lfc of 1.0 
sig_lipids_1 <- de_results_1 %>% filter(P.Value <= 0.05, abs(logFC) >= 1.0) %>%
  select(everything())

sig_lipids_1 # 27 lipid species differ
#write.csv(sig_lipids_1, "~/Desktop/sig_lipids_APPPS1-Veh_vs_NTG-Veh.csv", row.names = F)

# plot the volcano plotting
maxY <- max(-log10(de_results_1$P.Value), na.rm = TRUE)

volcano_1 <- EnhancedVolcano(
  de_results_1,
  lab = de_results_1$Molecule,
  x   = 'logFC',
  y   = 'P.Value',
  pCutoff  = 0.05, 
  FCcutoff = 1.0,
  title    = 'APPPS1-Veh vs NTG-Veh',
  subtitle = NULL,
  labSize = 5,
  labCol = 'blue',
  labFace = 'bold',
  boxedLabels = TRUE,
  parseLabels = FALSE,
  col = c('black', 'pink', 'green3', 'red3'),
  colAlpha = 4/5,
  legendPosition = 'bottom',
  legendLabSize = 10,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.3,
  colConnectors = 'black', max.overlaps = Inf, pointSize = 2) + coord_cartesian(ylim = c(0, maxY + 1))
#ggsave("~/Desktop/volcano_APPS1-Veh_vs_NTG-Veh_lipidomics.png", plot = volcano_1, dpi = 600, limitsize = F, bg ="white")

# 2- APPPS1-CP2 vs APPPS1-Veh
de_results_2 = de_analysis(data= d_normalized, APPPS1_CP2 - APPPS1_Veh, group_col = "Group",
                           measure = "Area")

# cut off raw p-value 0.05, lfc of 1.0
sig_lipids_2 <- de_results_2 %>% filter(P.Value <= 0.05, abs(logFC) >= 1.0) %>%
  select(everything())

sig_lipids_2 # 2 lipid species differ
#write.csv(sig_lipids_2, "~/Desktop/sig_lipids_APPPS1-CP2_vs_APPS1-Veh.csv", row.names = F)

# plot the volcano plotting
maxY <- max(-log10(de_results_2$P.Value), na.rm = TRUE)

volcano_2 <- EnhancedVolcano(
  de_results_2,
  lab = de_results_2$Molecule,
  x   = 'logFC',
  y   = 'P.Value',
  pCutoff  = 0.05, 
  FCcutoff = 1.0,
  title    = 'APPPS1-CP2 vs APPPS1-Veh',
  subtitle = NULL,
  labSize = 4,
  labCol = 'blue',
  labFace = 'bold',
  boxedLabels = TRUE,
  parseLabels = FALSE,
  col = c('black', 'pink', 'green3', 'red3'),
  colAlpha = 4/5,
  legendPosition = 'bottom',
  legendLabSize = 10,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.3,
  colConnectors = 'black', max.overlaps = Inf, pointSize = 2) + coord_cartesian(ylim = c(0, maxY + 1))
#ggsave("~/Desktop/volcano_APPS1-CP2_vs_APPS1-Veh_lipidomics.png", plot = volcano_2, dpi = 600, limitsize = F, bg ="white")

# 3- APPPS1-CP2 vs NTG-Veh
de_results_3 = de_analysis(data= d_normalized, APPPS1_CP2 - NTG_Veh, group_col = "Group",
                           measure = "Area")

# cut off raw p-value 0.05, lfc of 1.0
sig_lipids_3 <- de_results_3 %>% filter(P.Value <= 0.05, abs(logFC) >= 1.0) %>%
  select(everything())

sig_lipids_3 # 25 lipid species differ
#write.csv(sig_lipids_3, "~/Desktop/sig_lipids_APPPS1-CP2_vs_NTG-Veh.csv", row.names = F)

# plot the volcano plotting
maxY <- max(-log10(de_results_3$P.Value), na.rm = TRUE)

volcano_3 <- EnhancedVolcano(
  de_results_3,
  lab = de_results_3$Molecule,
  x   = 'logFC',
  y   = 'P.Value',
  pCutoff  = 0.05, 
  FCcutoff = 1.0,
  title    = 'APPPS1-CP2 vs NTG-Veh',
  subtitle = NULL,
  labSize = 4,
  labCol = 'blue',
  labFace = 'bold',
  boxedLabels = TRUE,
  parseLabels = FALSE,
  col = c('black', 'pink', 'green3', 'red3'),
  colAlpha = 4/5,
  legendPosition = 'bottom',
  legendLabSize = 10,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.3,
  colConnectors = 'black', max.overlaps = Inf, pointSize = 2) + coord_cartesian(ylim = c(0, maxY + 1))
#ggsave("~/Desktop/volcano_APPS1-CP2_vs_NTG-Veh_lipidomics.png", plot = volcano_3, dpi = 600, limitsize = F, bg ="white")

############################### lipid class analysis
# visualizing the total lipid class levels across treatment group
# extract the normalized values 
assay_mat <- assay(d_normalized, "Area")

# lipid classes
lipid_info <- as.data.frame(rowData(d_normalized)) %>%
  mutate(Feature = rownames(.)) %>% 
  select(Feature, Molecule, Class)

# sample annotation
sample_info <- as.data.frame(colData(d_normalized)) %>%
  mutate(Sample = rownames(.))

# convert to longer format 
df_long <- assay_mat %>%
  as.data.frame() %>%
  mutate(Feature = rownames(.)) %>%
  pivot_longer(
    cols = -Feature,
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  left_join(lipid_info,  by = "Feature") %>%
  left_join(sample_info, by = "Sample")

# summarise per class and drop NA
class_sample <- df_long %>%
  filter(!is.na(Class)) %>%               
  group_by(Class, Sample, Group) %>%
  summarise(
    ClassAbundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop")

class_sample$Group <- gsub("_", "-", class_sample$Group)
class_sample$Group <- factor(class_sample$Group, levels = c("NTG-Veh", "APPPS1-Veh", "APPPS1-CP2"))

# plotting
my_colors <- new_color_scheme(c("#31B57BFF", "#4DACD6", "#fb8500"), name = "new_color")

p <- class_sample |> 
  tidyplot(x = Group, y = ClassAbundance, color = Group) |>
  add_boxplot(alpha = 1, fill = "white", saturation = 1, linewidth = 0.85) |> 
  add_data_points_beeswarm(color = "black", size = 2) |>
  add_test_asterisks(method = "tukey_hsd", hide_info = TRUE) |> 
  adjust_colors(new_colors = my_colors) |>
  remove_x_axis_title() |>
  adjust_x_axis(rotate_labels = 45) |>
  adjust_y_axis_title("Lipid class abundance") |>
  adjust_font(fontsize = 14) |>
  adjust_title(fontsize = 16, face = "bold")|>
  split_plot(by = Class, nrow = 5, ncol = 5) 

#ggsave("~/Desktop/lipids_classes.png", plot = p, dpi = 600, width = 30, height = 30, bg = "white")

################################## lipids heatmaps
# prepare the annotations
lipid_info <- lipid_info[match(rownames(assay_mat), lipid_info$Feature), ]
sample_info <- sample_info[match(colnames(assay_mat), sample_info$Sample), ]

# transpose the normalized matrix
mat_t <- t(assay_mat)   # now: rows = samples, cols = Lipids

# remove the lipids of no class idnetified 
keep_idx <- !is.na(lipid_info$Class)
mat_t    <- mat_t[, keep_idx]
feat_info_sub <- lipid_info[keep_idx, ]

# Z-score each lipid across samples (so colors are comparable)
mat_scaled <- scale(mat_t) 

# Group as factor, ordered
sample_info$Group <- gsub("_", "-", sample_info$Group)
sample_info$Group <- factor(
  sample_info$Group,
  levels = c("NTG-Veh", "APPPS1-Veh", "APPPS1-CP2"))

# row annotation: sample groups
row_ha <- rowAnnotation(
  Group = sample_info$Group,
  col = list(
    Group = c(
      "NTG-Veh"    = "#009E73",
      "APPPS1-Veh" = "#0072B2",
      "APPPS1-CP2" = "#D55E00"
    )
  ),
  show_annotation_name = TRUE
)

# column annotation: lipid class
classes <- sort(unique(feat_info_sub$Class))

class_colors <- c(
  AcCa    = "#8DD3C7",
  LPE     = "#FB8072",
  LPC     = "#00FF7F",  
  SPH     = "#984EA3",
  PE      = "#0000FF",
  FA      = "#7FFFD4",
  LPI     = "#A65628",
  TG      = "#FFBBFF",  
  PG      = "#7FFF00",  
  PI      = "#FF4040",
  PC      = "#999999",  
  PS      = "#A2CD5A",
  SM      = "#EEDC82", 
  Hex1Cer = "#CC79A7",  
  Cer     = "#D55E00",  
  Hex2Cer = "#1B9E77",
  DG      = "#56B4E9", 
  DGO     = "#7570B3",
  TGO     = "#E7298A",
  CL      = "#66A61E",
  ChE     = "#FF7F00"   
)

col_ha <- HeatmapAnnotation(
  Class = feat_info_sub$Class,
  col = list(
    Class = class_colors
  ),
  show_annotation_name = TRUE
)

# color scale for the heatmap 
col_fun <- colorRamp2(
  c(-2, 0, 2),                        
  c("#000000", "#8968CD", "#FFFF66")   
)

heatmap_lipid <- Heatmap(
  mat_scaled,
  name = "z-score",
  col  = col_fun,
  top_annotation    = col_ha,
  left_annotation   = row_ha,
  cluster_rows      = TRUE,
  cluster_row_slices = TRUE,
  row_split         = sample_info$Group,
  cluster_columns   = TRUE,
  cluster_column_slices = TRUE,
  column_split      = feat_info_sub$Class,
  show_row_names    = TRUE,           
  show_column_names = FALSE,           
  column_title      = "Lipid species (grouped by class)",
  row_title         = "Treatment groups",
  column_title_side = "bottom",
  heatmap_legend_param = list(
    title = "z-score\n(normalized abundance)"))

png(
  filename = "~/Desktop/lipidomics_normalized_heatmap.png",
  width  = 17.2,
  height = 5.79,
  units  = "in",
  res    = 600,
  bg     = "white"
)

# draw the heatmap object
draw(
  heatmap_lipid,
  heatmap_legend_side    = "right",
  annotation_legend_side = "right"
)

dev.off()

#################### heatmap for cholesterol species only 
cholest_id <- lipid_info %>% filter(Class == "ChE")
cholest_species <- rownames(cholest_id)

matrix_cholest <- assay_mat[rownames(assay_mat) %in% cholest_species, , drop = FALSE]
matrix_cholest <- t(matrix_cholest)
matrix_cholest <- scale(matrix_cholest)

row_ha_cholest <- rowAnnotation(
  Group = factor(sample_info$Group, levels = c("NTG-Veh", "APPPS1-Veh", "APPPS1-CP2")),
  col = list(
    Group = c(
      "NTG-Veh"    = "#009E73",
      "APPPS1-Veh" = "#0072B2",
      "APPPS1-CP2" = "#D55E00"
    )
  ),
  show_annotation_name = TRUE, annotation_name_rot = 45
)

choles_map <- Heatmap(
  matrix_cholest,
  name = "z-score",
  col  = col_fun,
  cluster_rows      = TRUE,
  cluster_row_slices = FALSE,
  row_split         = sample_info$Group,
  left_annotation   = row_ha_cholest,
  cluster_columns   = TRUE,
  cluster_column_slices = TRUE,
  show_row_names    = TRUE,           
  show_column_names = TRUE,           
  column_title      = "Cholesterol species",
  row_title         = "Treatment groups",
  column_title_side = "bottom",
  heatmap_legend_param = list(
    title = "z-score\n(normalized abundance)"), column_names_rot = 45)


png(
  filename = "~/Desktop/cholesterol_normalized_heatmap.png",
  width  = 17.2,
  height = 5.79,
  units  = "in",
  res    = 600,
  bg     = "white"
)

# draw the heatmap object
draw(
  choles_map,
  heatmap_legend_side    = "right",
  annotation_legend_side = "right")
dev.off()

###################### erichment analysis [output of LION]
apps1_veh <- read_csv("/Users/egabal/Library/CloudStorage/Box-Box/Papers preparations/CP2 paper/supplmentary files/lipidomics/lion_enrichment_APPPS1-Veh_lipidomics.csv")
apps1_cp1 <- read_csv("/Users/egabal/Library/CloudStorage/Box-Box/Papers preparations/CP2 paper/supplmentary files/lipidomics/lion_enrichment_APPPS1-CP2_vs_NTG_lipidomics.csv")

apps1_cp1$`-log FDR (Q-value)` <- -log10(apps1_cp1$`FDR q-value`)
apps1_cp1$Discription <- fct_reorder(apps1_cp1$Discription,
                                     apps1_cp1$`-log FDR (Q-value)`)
plot_cp2 <- ggplot(apps1_cp1,
       aes(x = Discription,
           y = `-log FDR (Q-value)`,
           fill = `-log FDR (Q-value)`)) +
  geom_col(color = "black") +
  coord_flip() +
  scale_fill_gradient(low = "#FDBE85", high = "#D7301F") +
  scale_size_continuous(
    range = c(7, 35),            
    name = "-log10(FDR q-value)") + 
  labs(x = NULL,
       y = "-log10(FDR q-value)",
       title = "APP/PS1-CP2 vs NTG-Veh ") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.x = element_text(size = 20, color = "black"), axis.text.y = element_text(size = 20, color = "black"), axis.title.x = element_text(size = 18))

#ggsave("~/Desktop/apps1_cp2_enrcihment.png", plot = plot_cp2, dpi = 600, limitsize = F, bg = "white")


apps1_veh$`-log FDR (Q-value)` <- -log10(apps1_veh$`FDR q-value`)
apps1_veh$Discription <- fct_reorder(apps1_veh$Discription,
                                     apps1_veh$`-log FDR (Q-value)`)
plot_veh <- ggplot(apps1_veh,
                   aes(x = Discription,
                       y = `-log FDR (Q-value)`,
                       fill = `-log FDR (Q-value)`)) +
  geom_col(color = "black") +
  coord_flip() +
  scale_fill_gradient(low = "#FDBE85", high = "#D7301F") +
  scale_size_continuous(
    range = c(7, 35),            
    name = "-log10(FDR q-value)") + 
  labs(x = NULL,
       y = "-log10(FDR q-value)",
       title = "APP/PS1-Veh vs NTG-Veh ") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.x = element_text(size = 20, color = "black"), axis.text.y = element_text(size = 20, color = "black"), axis.title.x = element_text(size = 18))

#ggsave("~/Desktop/apps1_veh_enrcihment.png", plot = plot_veh, dpi = 600, limitsize = F, bg = "white")
