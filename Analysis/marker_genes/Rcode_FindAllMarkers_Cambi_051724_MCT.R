library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(scales)
library(RColorBrewer)
options(future.globals.maxSize = 3e+09)

setwd("F:/Colab/Cambi_lab/Data/rds")
seurat_obj <- readRDS("Cambi_051724_final.Rds")
seurat_obj #19059 features across 115555 samples

seurat_obj <- JoinLayers(seurat_obj)
seurat_obj <- NormalizeData(seurat_obj)
Idents(seurat_obj) <- "MCT"

MCT_order <- c("Astro-TE", "Astro-NT", "Ependymal", "CHOR",
               "IT-ET Glut", "NP-CT-L6b Glut", 
               "DG Glut", "DG-PIR Ex IMN", "OB-STR-CTX Inh IMN", "HPF CR Glut",
               "OPC", "COP", "NFOL", "MFOL", "MOL",
               "Endo", "VLMC", "Peri", "SMC", "ABC",
               "Microglia", "BAM", "Monocytes", "DC", "T cells", "B cells")

setdiff(unique(seurat_obj$MCT), MCT_order)
levels(seurat_obj) <- MCT_order

## FindAllMarkers
obj_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, 
                              max.cells.per.ident = 500)

write.csv(obj_markers, file = "F:/Colab/Cambi_lab/Analysis/marker_genes/Cambi_051724_MCT_markers.csv")
obj_markers <- read.csv("F:/Colab/Cambi_lab/Analysis/marker_genes/Cambi_051724_MCT_markers.csv")
obj_markers <- subset(obj_markers, obj_markers$p_val_adj < 0.05)

obj_markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top100
write.csv(top100, file = "F:/Colab/Cambi_lab/Analysis/marker_genes/Cambi_051724_MCT_top_100_markers.csv")

obj_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
write.csv(top10, file = "F:/Colab/Cambi_lab/Analysis/marker_genes/Cambi_051724_MCT_top_10_markers.csv")

select_genes <- subset(top100, top100$pct.1 > 0.8 & top100$pct.2 < 0.05)

#####
#####

MCT_colors_df <- read.csv("F:/Colab/Cambi_lab/Figures/colors/MCT_colors.csv")

features <- c(
  "Aqp4",
  "Agt",
  "Tmem212",
  "Clic6",
  "Cck",
  "Hs3st4",
  "Fam163b",
  "Igfbpl1",
  "Insm1",
  "Reln",
  "Vcan",
  "Gpr17",
  "Rassf10",
  "Ctps",
  "Opalin",
  "Cldn5",
  "Slc6a13",
  "Higd1b",
  "Tpm2",
  "Slc47a1",
  "Tmem119",
  "F13a1",
  "Cxcr2",
  "Cd74",
  "Ms4a4b",
  "Cd79a")

feature_to_MCT <- c(
  "Aqp4" = "Astro-TE",
  "Agt" = "Astro-NT",
  "Tmem212" = "Ependymal",
  "Clic6" = "CHOR",
  "Cck" = "IT-ET Glut",
  "Hs3st4" = "NP-CT-L6b Glut",
  "Fam163b" = "DG Glut",
  "Igfbpl1" = "DG-PIR Ex IMN",
  "Insm1" = "OB-STR-CTX Inh IMN",
  "Reln" = "HPF CR Glut",
  "Vcan" = "OPC",
  "Gpr17" = "COP",
  "Rassf10" = "NFOL",
  "Ctps" = "MFOL",
  "Opalin" = "MOL",
  "Cldn5" = "Endo",
  "Slc6a13" = "VLMC",
  "Higd1b" = "Peri",
  "Tpm2" = "SMC",
  "Slc47a1" = "ABC",
  "Tmem119" = "Microglia",
  "F13a1" = "BAM",
  "Cxcr2" = "Monocytes",
  "Cd74" = "DC",
  "Ms4a4b" = "T cells",
  "Cd79a" = "B cells")

mapped_MCTs <- feature_to_MCT[features]
colors_to_use <- MCT_colors_df$color[match(mapped_MCTs, MCT_colors_df$MCT)]
names(colors_to_use) <- features

## VlnPlot
vln <- VlnPlot(seurat_obj, features = features, layer = "data",
              stack = TRUE, sort = FALSE, flip = TRUE, cols = colors_to_use) +
  geom_violin(color = NA) +  
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 12)
  )

ggsave("F:/Colab/Cambi_lab/Figures/Vlnplot_MCT_markers.pdf", 
       plot = vln, width = 8, height = 8, dpi = 300)
