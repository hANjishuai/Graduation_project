# 00_rawdatamerge.R

# 从命令行参数中读取输入和输出
projectname <- snakemake@config[["projectname"]]  # 第一个参数：projectname
rawdatadir <- snakemake@config[["rawdatadir"]]  # 第二个参数：rawdatadir
tissues <- snakemake@config[["tissues"]]     # 第三个参数：tissues
omics <- snakemake@config[["omics"]]        # 第四个参数：omics
integrated_pdf <- as.character(snakemake@config[["output"]]["figure"])               # 第五个参数：integrated_pdf
seurat_rdata <- as.character(snakemake@config[["output"]]["result"])                   # 第六个参数：seurat_rdata

if(T){
# 加载库
.libPaths("~/R/4.4.1/library/")
library(Seurat)
library(SeuratData)
library(DoubletFinder)
library(patchwork)
library(stringr)
library(ggplot2)
library(harmony)

# 创建输出目录
if (!dir.exists("figure_out/01.QC/")) dir.create("figure_out/01.QC/",recursive = TRUE)
if (!dir.exists("result_out/01.QC/")) dir.create("result_out/01.QC/",recursive = TRUE)


# 依次读入样本
samples <- dir(rawdatadir)

sceList <- list()
for (pro in samples) {
  folder <- file.path(rawdatadir,pro)
  pj_name <- str_split(pro,pattern = "_",simplify = T)[1]
  mobject <- CreateSeuratObject(counts =Read10X(folder), project = pj_name)
  sceList[[pj_name]] <- mobject
}

# 合并样本
cell_ids <- names(sceList)
seuratdata <- merge(
  x = sceList[[1]], 
  y = sceList[-1], 
  add.cell.ids = cell_ids, 
  project = projectname
)

# 不整合
seuratdata <- NormalizeData(seuratdata)
seuratdata <- FindVariableFeatures(seuratdata)
seuratdata <- ScaleData(seuratdata)
seuratdata <- RunPCA(seuratdata)
seuratdata <- RunHarmony(seuratdata, "orig.ident")
seuratdata <- RunUMAP(seuratdata, dims = 1:30, reduction = "harmony", reduction.name = "umap.unintegrated")
umap_ui <- DimPlot(seuratdata, reduction = "umap.unintegrated", group.by = c("orig.ident"))

# 执行整合
seuratdata <- IntegrateLayers(object = seuratdata, method = CCAIntegration,
                              orig.reduction = "harmony", new.reduction = "integrated.cca",
                              verbose = FALSE)
seuratdata[["RNA"]] <- JoinLayers(seuratdata[["RNA"]])
seuratdata <- FindNeighbors(seuratdata, reduction = "integrated.cca", dims = 1:30)
seuratdata <- FindClusters(seuratdata, resolution = 0.5)
seuratdata <- RunUMAP(seuratdata, dims = 1:30, reduction = "integrated.cca")
umap_i <- DimPlot(seuratdata, reduction = "umap", group.by = c("orig.ident"))

# 保存整合图
p <- umap_ui + umap_i
ggsave(integrated_pdf, plot = p, width = 12, height = 5)

# 去除双细胞
sweep.res.list <- paramSweep(seuratdata, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
annotations <- seuratdata@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075 * nrow(seuratdata@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
seuratdata <- doubletFinder(seuratdata, PCs = 1:30, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
colnames(seuratdata@meta.data)[7] <- "DF.classifications"
seuratdata <- subset(seuratdata, subset = `DF.classifications` == "Singlet")
seuratdata$condition <- gsub("\\d$","",  seuratdata$orig.ident)
# 保存 Seurat 对象
save(seuratdata, file = seurat_rdata)
}