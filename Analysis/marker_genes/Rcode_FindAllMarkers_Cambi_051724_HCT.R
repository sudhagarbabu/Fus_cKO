library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
options(future.globals.maxSize = 3e+09)

setwd("F:/Colab/Cambi_lab/Data/rds")
seurat_obj <- readRDS("Cambi_051724_final_decontX.Rds")
seurat_obj #19059 features across 115555 samples

seurat_obj <- JoinLayers(seurat_obj)
seurat_obj <- NormalizeData(seurat_obj)
Idents(seurat_obj) <- "HCT"

HCT_order <- c("Astro-TE", "Astro-NT", "Ependymal", "CHOR",
               "IT-ET Glut", "NP-CT-L6b Glut", 
               "DG Glut", "DG-PIR Ex IMN", "OB-STR-CTX Inh IMN", "HPF CR Glut",
               "OPC", "COP", "NFOL", "MFOL",
               "MOL_1", "MOL_2", "MOL_3", "MOL_4", "MOL_5", "MOL_6", "MOL_7",
               "Endo", "VLMC", "Peri", "SMC", "ABC",
               "Microglia_1", "Microglia_2", "Microglia_3", "Microglia_4", "Microglia_5", "Microglia_6", 
               "BAM", "Monocytes", "DC", "T cells", "B cells")

setdiff(unique(seurat_obj$HCT), HCT_order)
levels(seurat_obj) <- HCT_order

#####
#####

## Micro
Idents(seurat_obj) <- "MCT"
pbmc <- subset(seurat_obj, idents = "Microglia")
Idents(pbmc) <- "HCT"

HCT_order <- c("Microglia_1", "Microglia_2", "Microglia_3", "Microglia_4", "Microglia_5", "Microglia_6")

setdiff(unique(pbmc$HCT), HCT_order)
levels(pbmc) <- HCT_order

Microglia_markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, max.cells.per.ident = 500)
write.csv(Microglia_markers, file = "F:/Colab/Cambi_lab/Analysis/marker_genes/Cambi_051724_Microglia_markers.csv")

Microglia_markers <- read.csv("F:/Colab/Cambi_lab/Analysis/marker_genes/Cambi_051724_Microglia_markers.csv")
Microglia_markers <- subset(Microglia_markers, Microglia_markers$p_val_adj < 0.05)

Microglia_markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top100
write.csv(top100, file = "F:/Colab/Cambi_lab/Analysis/marker_genes/Cambi_051724_Microglia_top_100_markers.csv")

Microglia_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
write.csv(top10, file = "F:/Colab/Cambi_lab/Analysis/marker_genes/Cambi_051724_Microglia_top_10_markers.csv")

#####
#####

## MOL
Idents(seurat_obj) <- "MCT"
pbmc <- subset(seurat_obj, idents = "MOL")
Idents(pbmc) <- "HCT"

HCT_order <- c("MOL_1", "MOL_2", "MOL_3", "MOL_4", "MOL_5", "MOL_6", "MOL_7")

setdiff(unique(pbmc$HCT), HCT_order)
levels(pbmc) <- HCT_order

MOL_markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, max.cells.per.ident = 500)
write.csv(MOL_markers, file = "F:/Colab/Cambi_lab/Analysis/marker_genes/Cambi_051724_MOL_markers.csv")

MOL_markers <- read.csv("F:/Colab/Cambi_lab/Analysis/marker_genes/Cambi_051724_MOL_markers.csv")
MOL_markers <- subset(MOL_markers, MOL_markers$p_val_adj < 0.05)

MOL_markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top100
write.csv(top100, file = "F:/Colab/Cambi_lab/Analysis/marker_genes/Cambi_051724_MOL_top_100_markers.csv")

MOL_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
write.csv(top10, file = "F:/Colab/Cambi_lab/Analysis/marker_genes/Cambi_051724_MOL_top_10_markers.csv")

#####
#####

## Microglia plots
Idents(seurat_obj) <- "MCT"
obj <- subset(seurat_obj, idents = "Microglia")
obj <- ScaleData(obj)
Idents(obj) <- "HCT"

top10 <- Microglia_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10) %>%
  arrange(cluster, desc(avg_log2FC))

HCT_colors_df <- read.csv("F:/Colab/Cambi_lab/Figures/colors/HCT_colors.csv")
HCT_colors <- setNames(HCT_colors_df$color, HCT_colors_df$HCT)

hm <- DoHeatmap(subset(obj, downsample = 500), 
                features = unique(top10$gene), 
                group.by = "HCT", 
                slot = "scale.data",
                group.colors = HCT_colors, 
                size = 4) 
hm <- hm +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(
    axis.text.x = element_text(size = 10, face = "plain"),
    axis.text.y = element_text(size = 10, face = "plain"),
    legend.text = element_text(size = 10, face = "plain"),
    legend.title = element_text(size = 10, face = "plain"),
    strip.text = element_text(size = 10, face = "plain")
  )

pdf("F:/Colab/Cambi_lab/Figures/Heatmap_Microglia_markers.pdf", width = 10, height = 10)
hm
dev.off()

#####
#####

## MOL plots
Idents(seurat_obj) <- "MCT"
obj <- subset(seurat_obj, idents = "MOL")
obj <- ScaleData(obj)
Idents(obj) <- "HCT"

top10 <- MOL_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10) %>%
  arrange(cluster, desc(avg_log2FC))

HCT_colors_df <- read.csv("F:/Colab/Cambi_lab/Figures/colors/HCT_colors.csv")
HCT_colors <- setNames(HCT_colors_df$color, HCT_colors_df$HCT)

hm <- DoHeatmap(subset(obj, downsample = 500), 
                features = unique(top10$gene), 
                group.by = "HCT", 
                slot = "scale.data",
                group.colors = HCT_colors, 
                size = 4) 
hm <- hm +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(
    axis.text.x = element_text(size = 10, face = "plain"),
    axis.text.y = element_text(size = 10, face = "plain"),
    legend.text = element_text(size = 10, face = "plain"),
    legend.title = element_text(size = 10, face = "plain"),
    strip.text = element_text(size = 10, face = "plain")
  )

pdf(file = "F:/Colab/Cambi_lab/Figures/Heatmap_MOL_markers.pdf", width = 10, height = 10)
hm
dev.off()

#####
#####
