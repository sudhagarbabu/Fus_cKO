library(Seurat)
library(SingleCellExperiment)
library(muscat)
library(limma)
library(purrr)
library(dplyr)
options(future.globals.maxSize = 3e+09)

setwd("F:/Colab/Cambi_lab/Data/rds")
pbmc <- LoadSeuratRds("Cambi_051724_final_decontX.Rds")

ncells_by_geno <- data.frame(table(pbmc$genotype, pbmc$MCT))
ncells_by_geno <- subset(ncells_by_geno, ncells_by_geno$Freq > 0)
ncells_by_geno <- subset(ncells_by_geno, ncells_by_geno$Freq < 20)
celltypes <- as.character(unique(ncells_by_geno$Var2))

Idents(pbmc)= "MCT"
pbmc$msex <- ifelse(grepl("female", pbmc$sex), 0, 1)
pbmc <- subset(pbmc, idents= celltypes, invert=T)
pbmc

pbmc$genotype <- as.character(pbmc$genotype)
table(pbmc$genotype)

Idents(pbmc)= "genotype"
pbmc <- subset(pbmc, idents= c("AD/cKO", "AD"))
pbmc$genotype_num <- ifelse(grepl("AD/cKO", pbmc$genotype), 1, 0)
table(pbmc$genotype_num)

## celltypes
ncells_per_celltype <- data.frame(table(pbmc$MCT))
ncells_per_celltype <- ncells_per_celltype[order(ncells_per_celltype$Freq, decreasing = T), ]
celltypes <- as.character(ncells_per_celltype$Var1)

for (i in 1:length(celltypes)){
  
  celltype <- celltypes[[i]]
  Idents(pbmc) <- "MCT"
  pbmc_subset <- subset(pbmc, idents = celltype)
  pbmc_subset <- JoinLayers(pbmc_subset)
  
  metadata <- pbmc_subset@meta.data
  pbmc.sce <- SingleCellExperiment(list(counts= pbmc_subset[["RNA"]]$counts), colData= metadata)
  rm(pbmc_subset)
  
  pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 0) > 0, ]
  pbmc.sce <- pbmc.sce[rowSums(counts(pbmc.sce) > 1) >= 10, ]
  dim(pbmc.sce)
  
  sce <- prepSCE(pbmc.sce, kid = "MCT", gid = "genotype", sid = "sample", drop = FALSE)
  nk <- length(kids <- levels(sce$cluster_id))
  ns <- length(sids <- levels(sce$sample_id))
  names(kids) <- kids; names(sids) <- sids
  
  pb <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))
  assayNames(pb)
  
  rm(sce)
  rm(pbmc.sce)
  
  formula <- as.formula(paste0("~ genotype_num + age_max + msex"))
  cd <- as.data.frame(colData(pb))
  cd2 <- cd[, c("genotype_num", "age_max", "msex")]
  design <- model.matrix(formula, cd2)
  pb2 <- pb[, rownames(design)]
  
  ## muscat
  res <- pbDS(pb2, design = design, coef = as.list(2:4), min_cells = 5)
  
  tbl <- res$table[[1]]
  DEG_tbl <- tbl[[1]]
  DEG_tbl <- DEG_tbl[order(DEG_tbl$p_adj.loc), ]
  
  outdir <- paste0("F:/Colab/Cambi_lab/Analysis/muscat/AD_cKO_vs_AD/MCT")
  dir.create(outdir, recursive = T)
  setwd(outdir)
  celltype <- gsub(" ","_",celltype)
  filename <- paste0(celltype,"_AD_cKO_vs_AD.csv")
  write.csv(DEG_tbl, file= filename)
}
