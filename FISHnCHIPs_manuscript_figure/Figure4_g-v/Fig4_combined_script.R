setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
pacman::p_load(Seurat, dplyr, plyr, capFISHImage, ggplot2, dittoSeq, ComplexHeatmap)

### Loading data --------------------------------------------------
seurat <- readRDS("./input_files/seurat_anno_renamed.rds")
output.dir <- "./output_files/"
dir.create(output.dir, recursive = T, showWarnings = F)
unique(seurat@meta.data$cluster)

celltype.order <- c("L2/3", "L3/4", "L4/5", "L5p1", "L5/6", "L6p1", 
                    "IntPv", "IntSst", "IntNpy/CckVip", "Hip", "Sub")
n.exc <- 6
length(celltype.order)

seurat@meta.data$cluster <- factor(seurat@meta.data$cluster, 
                                   levels = celltype.order)
Idents(seurat) <- "cluster"
levels(seurat)
crop.name <- "halfImgHip"
x_lims <- c(-1625, 23486)
y_lims <- c(0, 17000)

h <- x_lims[2]-x_lims[1]
w <- 17000

cells.info <- seurat@meta.data

### # update colors  --------------------------------------------
cluster.colors <- c("#0804f9", "#bc13fe", "#ff9408",  
                             "#738595", "#ff000d", "#CCCFCB",  
                             "#26f7fd", "#21fc0d", "#ffff14",
                             "#047495", "#009E73")
                             # cluster.colors <- dittoColors()

exc.colors <- cluster.colors[1:6]
int.colors <- cluster.colors[7:11]
#----------------------------------------------------------------


### Plottting ---------------------------------------------------------------

### Heatmap for all bits -------------------------------------------------
cbm.allbits <- read.csv("./input_files/image_cbm_norm_20Bits_renamed.csv") 
rownames(cbm.allbits) <- cbm.allbits$X
dim(cbm.allbits)

cbm.allbits <- cbm.allbits %>% filter(X %in% rownames(cells.info))
cbm.allbits$X <- NULL
dim(cbm.allbits)
dim(seurat)
cells.info.new <- cells.info[rownames(cbm.allbits),]

seurat.allbits <- CreateSeuratObject(counts = t(cbm.allbits),
                                     project = "Mouse Cortex 20 Bits",
                                     meta.data = cells.info.new,
                                     min.cells = 0,
                                     min.features = 0)
seurat.allbits <- ScaleData(seurat.allbits, features = rownames(seurat.allbits))

slot = "scale.data"
bits.order <- c("ExcL2", "ExcL3", "ExcL4",
                "ExcL5p1", "ExcL5p2", "ExcL5p3",
                "ExcL6p1", "ExcL6p2", "IntPv",
                "IntSst", "IntNpy", "IntCckVip")
activity.bits <- c("Hip", "Sub", 
                   "Erp", "LrpD", 
                   "LrpS", "NS", 
                   "Other", "Syn")
bits.order.all <- append(bits.order, activity.bits)

dh <- dittoHeatmap(seurat.allbits,
                   genes = bits.order.all,
                   slot = slot, breaks = seq(-3, 3, length.out = 256),
                   annot.colors = cluster.colors,
                   heatmap.colors = colorRampPalette(c("magenta", "black", "yellow"))(256),
                   show_rownames = TRUE, cluster_cols = FALSE, cluster_rows = FALSE,
                   scale = "none", scaled.to.max = FALSE, annot.by = "cluster",
                   main = "", name = "Expression \n (z-scaled)", complex = T)
pdf(file = paste0(output.dir, "Fig4-Heatmap_20bits.pdf"), width = 8, height = 4)
draw(dh)
dev.off()


# Clustermap

# Clustermap for Exc and Inv
Exc.cells <- cells.info %>% filter(cluster %in% celltype.order[1:n.exc]) 
Exc.cells %>% 
  ggplot(aes(y, x,  colour = cluster)) +
  geom_point(size = 0.5) +
  facet_wrap(Exc.cells$cluster, nrow = 1) +
  scale_color_manual(values= exc.colors) +
  theme(panel.background = element_rect(fill = "black", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "black"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "black"),
        legend.position = "none",
        axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
ggsave(filename = paste0(output.dir, "Clustermap-Exc.pdf"),  height = 4*h/w, width = 4*n.exc)

other.cells <- cells.info %>% filter(!cluster %in% celltype.order[1:n.exc])
other.cells %>% 
  ggplot(aes(y, x,  colour = cluster)) +
  geom_point(size = 0.5) +
  facet_wrap(other.cells$cluster, nrow = 1) +
  scale_color_manual(values= int.colors) +
  theme(panel.background = element_rect(fill = "black", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
        axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
ggsave(filename = paste0(output.dir, "Clustermap-others.pdf"), height = 4*h/w, width = 4*(length(celltype.order)-n.exc))

#in the end arc was done with lots ofgrouping and ungrouping and aligning in Illustrator



# Fused clustermap 
cells.info %>%
  ggplot(aes(y, x,  colour = cluster)) +
  geom_point(size = 1) +
  scale_color_manual(values=cluster.colors) +
  theme(panel.background = element_rect(fill = "black", colour = "black", size = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
        axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
ggsave(filename = paste0(output.dir, "Fig4-Clustermap_fused_all.pdf"), height = 4*h/w, width = 4)




# UMAP
DimPlot(seurat, reduction = "umap", cols = cluster.colors,
        label = T, pt.size=1.5, repel = T, label.size = 5) 
ggsave(paste0(output.dir,"UMAP_all.pdf"), width = 8, height = 6)


### ------------------------- Gratident Analyisis ---------------------------------------------------------------
### Define edges and plot distance distribution 
gg_circle <- function(r, xc, yc, color="black", fill=NA, ...) {
  x <- xc + r*cos(seq(0, pi, length.out=100))
  ymax <- yc + r *sin(seq(0, pi, length.out=100))
  ymin <- yc + r *sin(seq(0, -pi, length.out=100))
  annotate("ribbon", x=x, ymin=ymin, ymax=ymax, color=color, fill=fill, ...)
}

fig.fuse <- cells.info %>% ggplot(aes(y, x)) +
  geom_point(aes(colour = cluster), size = 0.1) +
  scale_color_manual(values= cluster.colors) +
  theme(panel.background = element_rect(fill = "black", colour = "black", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_blank(), legend.position = "none",
        axis.ticks = element_blank(), axis.text = element_blank()
  )
print(fig.fuse)


# ### Definition (circles with the same R) --------------------------------- #used this for paper!!!
xc = 25500
yc = 10000
R = 25500
span = 10000
fig.fuse +
  gg_circle(r=R, xc= xc, yc= yc, color="white", alpha=0.2, size= 1) +
  gg_circle(r=R, xc =xc+span, yc= yc, color="white", alpha=0.2, size=1)


cells.info$r1 = sqrt((cells.info$y - xc)**2 + (cells.info$x - yc)**2)
cells.info$r2 = sqrt((cells.info$y - (xc+span) )**2 + (cells.info$x - yc)**2)
cells.info$d = R -cells.info$r1 # distance to the outer edge
cells.info$d.norm = cells.info$d/span

cells.in.between <- cells.info %>% filter((r1-R)*(r2-R)<0)


#Plot cell distance distribution
cells.in.between %>%
  select(cluster, d.norm) %>%
  filter(cluster %in% celltype.order[1:n.exc]) %>%
  ggplot(aes(x = d.norm, group = cluster)) +
  stat_density(aes(x=d.norm, colour=cluster), lwd = 2,
               geom="line",position="identity") +
  xlab("Normalised distance to the edge") + ylab("density") +
  scale_color_manual(values= exc.colors) +
  # xlim(0, 1) +
  theme_bw(base_size = 16) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(paste0(output.dir,"Fig4-cells_distance_distribution_ExcLCells_sameR.pdf"), width = 8, height = 4)


cells.in.between %>%
  select(cluster, d.norm) %>%
  filter(cluster %in% celltype.order[(n.exc+1): (n.exc+3)]) %>%
  ggplot(aes(x = d.norm, group = cluster)) +
  stat_density(aes(x=d.norm, colour=cluster), lwd = 2,
               geom="line",position="identity") +
  xlab("Normalised distance to the edge") + ylab("density") +
  scale_color_manual(values= int.colors) +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(paste0(output.dir,"Fig4-cells_distance_distribution_IntCells_sameR.pdf"), width = 8, height = 4)


