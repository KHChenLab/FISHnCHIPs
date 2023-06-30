setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
pacman::p_load(dplyr, plyr, Seurat, ggplot2, dittoSeq, ComplexHeatmap, tidyverse, ggpubr, ggrepel)

### Loading data and setup output directory ------------------------------------
seurat <- readRDS("./input_files/Fig3-MouseBrain_FISHnCHIPs_seurat_renamed.rds")
output.dir <- "./output_files/"
dir.create(output.dir, recursive = T, showWarnings = F)


## Assign clusternames 
cell.labels <- seurat@meta.data %>% select(seurat_clusters)
colnames(cell.labels) <- "cluster_id"
from.clusters <- levels(cell.labels[, "cluster_id"])
print(from.clusters)
new.cluster.names <- c("Glutamatergic Neurons",
                       "Glutamatergic Neurons",
                       "NA", 
                       "NA",
                       "Astrocytes",
                       "Oligodendrocytes",
                       "Glutamatergic Neurons",
                       "Astrocytes", 
                       "GABAergic Neurons",
                       "Oligodendrocytes",
                       "Endothelial",
                       "Macrophage",
                       "Perivascular cells",
                       "Vascular lepotomeningeal cells")
print(length(new.cluster.names))
show.clusters <- unique(new.cluster.names[new.cluster.names!="NA"])
print(show.clusters)

### ---------------------------------------------------------------------------
seurat@meta.data$Celltype <- mapvalues(cell.labels[,"cluster_id"], from = from.clusters, to = new.cluster.names)
seurat <- subset(seurat, Celltype %in% show.clusters)
seurat@meta.data$Celltype <- factor(seurat@meta.data$Celltype, 
                                    levels = c("Glutamatergic Neurons", 
                                               "GABAergic Neurons",
                                               "Astrocytes",
                                               "Oligodendrocytes",
                                               "Endothelial",
                                               "Macrophage",
                                               "Perivascular cells",
                                               "Vascular lepotomeningeal cells"
                                    ))
Idents(seurat) <- "Celltype"
#write.csv(seurat@meta.data, file = paste0(output.dir, paste0("MouseBrain_cells_meta_Anno.csv")), row.names = T)

celltypeAbundance <- as.data.frame(table(seurat@meta.data$Celltype))
colnames(celltypeAbundance) <- c("Celltype", "ncells")
#write.csv(celltypeAbundance, file = paste0(output.dir, paste0("MouseBrain_cellsAbundance_Anno.csv")), row.names = F)

### Plottting ---------------------------------------------------------------
cluster.colors <- c("#738595", "#bc13fe", "#0804f9", "#ff000d",
                             "#26f7fd", "#ff9408", "#21fc0d", "#ffff14")
                             
# Clustermap
cells.info <- seurat@meta.data
cells.info %>%
  ggplot(aes(y, x,  colour = Celltype)) +
  geom_point(size = 0.5) +
  facet_wrap(cells.info$Celltype, ncol = 4) +
  scale_color_manual(values= cluster.colors) +
  theme(panel.background = element_rect(fill = "black", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "black"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "black"),
        axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
ggsave(filename = "Fig3-Clustermap.pdf", path = output.dir, height = 4, width = 9)


### Heatmap
bits.order <- c("M12", "M13", "M7", "M18", "M2", "M1", 
                "M6", "M8", "M16", "M17", "M4", "M5", 
                "M3", "M10", "M11", "M9", "M14", "M15")
dh <- dittoHeatmap(seurat, genes = bits.order, slot = "scale.data",
                   breaks = seq(-3, 3, length.out = 256), annot.colors = cluster.colors,
                   heatmap.colors = colorRampPalette(c("magenta", "black", "yellow"))(256),
                   cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = TRUE,
                   scale = "none", scaled.to.max = FALSE, annot.by = "Celltype", 
                   main = "", name = "Expression \n (z-scaled)", complex = T)
filename <- paste0(output.dir, "Fig3-Heatmap.pdf")
pdf(file = filename, width = 8, height = 4)
draw(dh)
dev.off()


### Plotting cluster size correlation between MERFISH and FISHnCHIPs 
data <- read.csv("./input_files/Lm67_ClusterSize_Comparison_0421_subtract.csv")
rownames(data) <- data$Celltype
# data$Celltype <- NULL
# View(data)
row_sub = apply(data, 1, function(row) all(row !=0 ))


data[row_sub,] %>% ggplot(aes(x = capFISH, y = MERFISH)) + geom_point() +
  geom_point() + 
  # scale_x_continuous(limits=c(0, 0.6)) +
  # scale_y_continuous(imits=c(0, 0.6)) +
  xlim(0, 0.6) +
  xlim(0, 0.6) +
  geom_smooth(method='lm', formula= y~x, se = TRUE) +
  # stat_regline_equation(label.x=2.3, label.y=0.6) +
  # stat_cor(aes(label=..rr.label..), label.x=2.3, label.y=0.5)+
  # geom_text(label = names(row_sub), cex= 2, hjust = 0.4, vjust =-1, 
  #           check_overlap = TRUE, nudge_x = 0.01) +
  geom_text_repel(aes(label = rownames(data)), cex = 2.5) +
  labs(x = "Fraction of total cells in FISHnCHIPs", 
       y = "Fraction of total cells in MERFISH") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = "./corr_capFISH_pct_MERFISH_pct_point_removeEquations.pdf", path = output.dir,
       width = 6, height = 6)



