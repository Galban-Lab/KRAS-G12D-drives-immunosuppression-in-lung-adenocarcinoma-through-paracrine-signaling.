#This script recreates Figures 7A, 7C-F, S7A
#This data was downloaded from GEO (GSE136246)

#Maroni G, Bassal MA, Krishnan I, Fhu CW et al. Identification of a targetable 
#KRAS-mutant epithelial population in non-small cell lung cancer. Commun Biol 
#2021 Apr 14;4(1):370. PMID: 33854168

#Data was processed in line with the Seurat workflow:
#Website: https://satijalab.org/seurat/index.html

#Reference: Hao et al., Integrated analysis of multimodal single-cell data. 
#Cell. 2021 Jun 24;184(13):3573-3587.e29. doi: 10.1016/j.cell.2021.04.048. 
#PMID: 34062119; PMCID: PMC8238499. 

#Seurat Version 4.3.0.1
#SeuratObject Version 4.1.3
#"R version 4.2.3 (2023-03-15)" -- "Shortstop Beagle"

#Load packages
library(Seurat)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(scater)
library(pheatmap)


#### Pre-processing ----------------------------------------------------------------------------------------------------####
#Load Raw Data:
NSC010.t1.counts <- read.table("~/GSE136246_RAW/GSM4043243_NSC010.t1.counts.tsv", header = T, row.names = 1)
NSC010.t1.counts <- t(NSC010.t1.counts)
NSC010.t2.counts <- read.table("~/GSE136246_RAW/GSM4043244_NSC010.t2.counts.tsv", header = T, row.names = 1)
NSC010.t2.counts <- t(NSC010.t2.counts)

NSC016.t1.counts <- read.table("~/GSE136246_RAW/GSE136246_RAW 2/GSM4043245_NSC016.t1.counts.tsv", header = T, row.names = 1)
NSC016.t1.counts <- t(NSC016.t1.counts)
NSC016.t2.counts <- read.table("~/GSE136246_RAW/GSE136246_RAW 2/GSM4043246_NSC016.t2.counts.tsv", header = T, row.names = 1)
NSC016.t2.counts <- t(NSC016.t2.counts)
NSC016.t3.counts <- read.table("~/GSE136246_RAW/GSE136246_RAW 2/GSM4043247_NSC016.t3.counts.tsv", header = T, row.names = 1)
NSC016.t3.counts <- t(NSC016.t3.counts)

NSC020.t1.counts <- read.table("~/GSE136246_RAW/GSM4043253_NSC020.t1.counts.tsv", header = T, row.names = 1)
NSC020.t1.counts <- t(NSC020.t1.counts)
NSC020.t2.counts <- read.table("~/GSE136246_RAW/GSM4043254_NSC020.t2.counts.tsv", header = T, row.names = 1)
NSC020.t2.counts <- t(NSC020.t2.counts)

NSC035.counts <- read.table("~/GSE136246_RAW/GSM4043257_NSC035.counts.tsv", header = T, row.names = 1)
NSC035.counts <- t(NSC035.counts)

NSC037.counts <- read.table("~/GSE136246_RAW/GSM4043259_NSC037.counts.tsv", header = T, row.names = 1)
NSC037.counts <- t(NSC037.counts)



#Make Seurat Objects (include features/genes expressed by at least 3 cells and include cells with at least 100 features/genes)
NSC010.t1<- CreateSeuratObject(NSC010.t1.counts, min.cells = 3, min.features = 100)
NSC010.t2<- CreateSeuratObject(NSC010.t2.counts, min.cells = 3, min.features = 100)

NSC016.t1<- CreateSeuratObject(NSC016.t1.counts, min.cells = 3, min.features = 100)
NSC016.t2<- CreateSeuratObject(NSC016.t2.counts, min.cells = 3, min.features = 100)
NSC016.t3<- CreateSeuratObject(NSC016.t3.counts, min.cells = 3, min.features = 100)

NSC020.t1<- CreateSeuratObject(NSC020.t1.counts, min.cells = 3, min.features = 100)
NSC020.t2<- CreateSeuratObject(NSC020.t2.counts, min.cells = 3, min.features = 100)

NSC035<- CreateSeuratObject(NSC035.counts, min.cells = 3, min.features = 100)

NSC037<- CreateSeuratObject(NSC037.counts, min.cells = 3, min.features = 100)

#Add Metadata:
#Sample
NSC010.t1[["Sample"]] <- "NSC010.t1"
NSC010.t2[["Sample"]] <- "NSC010.t2"

NSC016.t1[["Sample"]] <- "NSC016.t1"
NSC016.t2[["Sample"]] <- "NSC016.t2"
NSC016.t3[["Sample"]] <- "NSC016.t3"

NSC020.t1[["Sample"]] <- "NSC020.t1"
NSC020.t2[["Sample"]] <- "NSC020.t2"

NSC035[["Sample"]] <- "NSC035"

NSC037[["Sample"]] <- "NSC037"

#Patient
NSC010.t1[["Patient"]] <- "NSC010"
NSC010.t2[["Patient"]] <- "NSC010"

NSC016.t1[["Patient"]] <- "NSC016"
NSC016.t2[["Patient"]] <- "NSC016"
NSC016.t3[["Patient"]] <- "NSC016"

NSC020.t1[["Patient"]] <- "NSC020"
NSC020.t2[["Patient"]] <- "NSC020"

NSC035[["Patient"]] <- "NSC035"

NSC037[["Patient"]] <- "NSC0237"

#Group
NSC010.t1[["Group"]] <- "KRAS-wt"
NSC010.t2[["Group"]] <- "KRAS-wt"

NSC016.t1[["Group"]] <- "KRAS-G12D"
NSC016.t2[["Group"]] <- "KRAS-G12D"
NSC016.t3[["Group"]] <- "KRAS-G12D"

NSC020.t1[["Group"]] <- "KRAS-G12D"
NSC020.t2[["Group"]] <- "KRAS-G12D"

NSC035[["Group"]] <- "KRAS-G12D"

NSC037[["Group"]] <- "KRAS-wt"

#Smoke
NSC010.t1[["Smoke"]] <- "Former Smoker"
NSC010.t2[["Smoke"]] <- "Former Smoker"

NSC016.t1[["Smoke"]] <- "Never Smoker"
NSC016.t2[["Smoke"]] <- "Never Smoker"
NSC016.t3[["Smoke"]] <- "Never Smoker"

NSC020.t1[["Smoke"]] <- "Former Smoker"
NSC020.t2[["Smoke"]] <- "Former Smoker"

NSC035[["Smoke"]] <- "Never Smoker"

NSC037[["Smoke"]] <- "Never Smoker"

#Merge Samples:
HumanLung <- merge(x = NSC010.t1, y = c(NSC010.t2,NSC016.t1,NSC016.t2,NSC016.t3,NSC020.t1,NSC020.t2,NSC035,NSC037), add.cell.ids = (c("a","b","c","d","e","f","g","h","i")))
#Normalize
HumanLung <- NormalizeData(object = HumanLung, normalization.method = "LogNormalize", scale.factor = 10000)

#Aply Unbiased QC Cutoffs
HumanLung[["percent.mt"]] <- PercentageFeatureSet(object = HumanLung, pattern = "^MT.")

HumanLung <- subset(x = HumanLung, subset = nCount_RNA > 800 & nCount_RNA < 100000 & percent.mt < 15)

#Integrate Datasets by Patient
Idents(object = HumanLung) <- "Patient"
HumanLung.list <- SplitObject(HumanLung, split.by = "Patient")

HumanLung.list <- lapply(X = HumanLung.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

HumanLung.anchors <- FindIntegrationAnchors(object.list = HumanLung.list)
HumanLung <- IntegrateData(anchorset = HumanLung.anchors)

DefaultAssay(HumanLung) <- "integrated"
HumanLung <- FindVariableFeatures(HumanLung, selection.method = "vst", nfeatures = 2000)

#Scale Data
all.genes <- rownames(HumanLung)
HumanLung <- ScaleData(HumanLung, verbose = T, features = all.genes)
HumanLung <- RunPCA(HumanLung, npcs = 30, verbose = FALSE)

#Calculate 90% Variance:
st_dev <- HumanLung@reductions$pca@stdev
var <- st_dev^2
sum(var[1:16])/ sum(var)

HumanLung <- FindNeighbors(object = HumanLung, dims = 1:16)
HumanLung <- FindClusters(object = HumanLung, resolution = 2)
HumanLung <- RunUMAP(HumanLung, reduction = "pca", dims = 1:16, verbose = F)
DimPlot(HumanLung, reduction = "umap", label = T, repel = T)
DimPlot(HumanLung, reduction = "umap", label = T, split.by = "Sample")
DimPlot(HumanLung, reduction = "umap", label = T, split.by = "Group")

#Determine Clusters:
DefaultAssay(HumanLung) <- "RNA"
DotPlot(HumanLung, features = c("PTPRC","CD79A","CD19","MS4A1","CD3E","CD4","CCR7","IL2RA","FOXP3","CTLA4","NKG7","KLRB1","CCR5","CD8A","ITGAM","CD14","LY6G6C","LY6H","S100A8","ACOD1","PDGFRA","IL6","CLEC3B","IL33","COL1A1","COL1A2","PDPN","CD68","ADGRE1","C1QA","APOE",
                                "MRC1","FCGR3A","HLA-DRB1","SIGLEC5","ITGAX","MERTK","SPP1","NR2F6","IGHA1","SDC1","ITGAE","FLT3","BATF3","XCR1","TAGLN","ACTA2","PDE5A","CLU","PECAM1","CDH5","FLT1","KIT","SNCA","SCGB1A1","LAMP3","CXCL8","AQP5","SFTPD","LYVE1","KRT8","KRT18","KRT19","SEC14L3",
                                "FOXJ1","MSLN","WT1","MYL7","TNNI3","TNNT2","G0S2","CMA1","HDC","GATA2","MS4A2","CCNA2","MKI67","HBB"), cols = "RdYlBu", dot.scale = 6, assay = "RNA") + RotatedAxis()

#Label Manual Clusters:
Idents(HumanLung) <- "seurat_clusters"
HumanLung <- RenameIdents(HumanLung,
                          "0" = "B Cell", 
                          "1" = "T Cell", 
                          "2" = "B Cell", 
                          "3" = "T Cell", 
                          "4" = "T Cell", 
                          "5" = "T Cell", 
                          "6" = "Macrophage", 
                          "7" = "T Cell", 
                          "8" = "Cycling",
                          "9" = "Macrophage", 
                          "10" = "Macrophage", 
                          "11" = "Macrophage", 
                          "12" = "Plasma Cell", 
                          "13" = "T Cell", 
                          "14" = "T Cell", 
                          "15" = "Epithelial", 
                          "16" = "Epithelial", 
                          "17" = "Monocyte",
                          "18" = "Cycling", 
                          "19" = "DC", 
                          "20" = "NK Cell", 
                          "21" = "Granulocyte",
                          "22" = "Plasma Cell",
                         "23" = "Plasma Cell", 
                         "24" = "T Cell", 
                         "25" = "Fibroblast", 
                         "26" = "Endothelial", 
                         "27" = "DC",
                         "28" = "Cycling", 
                         "29" = "Cycling", 
                         "30" = "Epithelial", 
                         "31" = "DC",
                         "32" = "B Cell")
HumanLung[["manual_clusters"]] <- HumanLung@active.ident

new_order <- c("Epithelial",
               "Fibroblast",
               "Macrophage",
               "Monocyte",
               "Granulocyte",
               "T Cell",
               "NK Cell",
               "DC",
               "Endothelial",
               "B Cell",
               "Plasma Cell",
               "Cycling")

HumanLung@active.ident <- factor(HumanLung@active.ident, levels = new_order)
HumanLung[["manual_clusters"]] <- HumanLung@active.ident

Idents(HumanLung) <- "Group"
new_order <- c("KRAS-G12D","KRAS-wt")
HumanLung@active.ident <- factor(HumanLung@active.ident, levels = new_order)

#### Pre-processing Epithelial Cells ----------------------------------------------------------------------------------------------------####
Idents(HumanLung) <- "manual_clusters"
HumanLungEpithelialthelial <- subset(HumanLung, idents = "Epithelial")
HumanLungEpithelial <- FindVariableFeatures(HumanLungEpithelial, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(HumanLungEpithelial)
HumanLungEpithelial <- ScaleData(HumanLungEpithelial, verbose = T, features = all.genes)
HumanLungEpithelial <- RunPCA(HumanLungEpithelial, npcs = 30, verbose = FALSE)

#Calculate 90% Variance:
st_dev <- HumanLungEpithelial@reductions$pca@stdev
var <- st_dev^2
sum(var[1:22])/ sum(var)

HumanLungEpithelial <- FindNeighbors(object = HumanLungEpithelial, dims = 1:22)
HumanLungEpithelial <- FindClusters(object = HumanLungEpithelial, resolution = 1)
HumanLungEpithelial <- RunUMAP(HumanLungEpithelial, reduction = "pca", dims = 1:22, verbose = F)
DimPlot(HumanLungEpithelial, reduction = "umap", label = T)

DefaultAssay(HumanLungEpithelial) <- "RNA"
DotPlot(HumanLungEpithelial, features = c("PTPRC","CD79A","CD19","MS4A1","CD3E","CD4","CCR7","IL2RA","FOXP3","CTLA4","NKG7","KLRB1","CCR5","CD8A","ITGAM","CD14","LY6G6C","LY6H","S100A8","ACOD1","PDGFRA","IL6","CLEC3B","IL33","COL1A1","COL1A2","PDPN","CD68","ADGRE1","C1QA","APOE",
                                   "MRC1","FCGR3A","HLA-DRB1","SIGLEC5","ITGAX","MERTK","SPP1","NR2F6","IGHA1","SDC1","ITGAE","FLT3","BATF3","XCR1","TAGLN","ACTA2","PDE5A","CLU","PECAM1","CDH5","FLT1","KIT","SNCA","SCGB1A1","LAMP3","CXCL8","AQP5","SFTPD","LYVE1","KRT8","KRT18","KRT19","SEC14L3",
                                   "FOXJ1","MSLN","WT1","MYL7","TNNI3","TNNT2","G0S2","CMA1","HDC","GATA2","MS4A2","CCNA2","MKI67","HBB"), cols = "RdYlBu", dot.scale = 6, assay = "RNA") + RotatedAxis()

#Remove non-epithelial cells
HumanLungEpithelial <- subset(HumanLungEpithelial, idents = c("7","8"), invert = T) #remove contamination

HumanLungEpithelial <- FindVariableFeatures(HumanLungEpithelial, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(HumanLungEpithelial)
HumanLungEpithelial <- ScaleData(HumanLungEpithelial, verbose = T, features = all.genes)
HumanLungEpithelial <- RunPCA(HumanLungEpithelial, npcs = 30, verbose = FALSE)

#Calculate 90% Variance:
st_dev <- HumanLungEpithelial@reductions$pca@stdev
var <- st_dev^2
sum(var[1:22])/ sum(var)

HumanLungEpithelial <- FindNeighbors(object = HumanLungEpithelial, dims = 1:22)
HumanLungEpithelial <- FindClusters(object = HumanLungEpithelial, resolution = 1)
HumanLungEpithelial <- RunUMAP(HumanLungEpithelial, reduction = "pca", dims = 1:22, verbose = F)
DimPlot(HumanLungEpithelial, reduction = "umap", label = T)
DimPlot(HumanLungEpithelial, reduction = "umap", label = T, split.by = "Group")
DimPlot(HumanLungEpithelial, reduction = "umap", label = T, split.by = "Sample")

DefaultAssay(HumanLungEpithelial) <- "RNA"

DotPlot(HumanLungEpithelial, features = c("PTPRC","CD79A","CD19","MS4A1","CD3E","CD4","CCR7","IL2RA","FOXP3","CTLA4","NKG7","KLRB1","CCR5","CD8A","ITGAM","CD14","LY6G6C","LY6H","S100A8","ACOD1","PDGFRA","IL6","CLEC3B","IL33","COL1A1","COL1A2","PDPN","CD68","ADGRE1","C1QA","APOE",
                                   "MRC1","FCGR3A","HLA-DRB1","SIGLEC5","ITGAX","MERTK","SPP1","NR2F6","IGHA1","SDC1","ITGAE","FLT3","BATF3","XCR1","TAGLN","ACTA2","PDE5A","CLU","PECAM1","CDH5","FLT1","KIT","SNCA","SCGB1A1","LAMP3","CXCL8","AQP5","SFTPD","LYVE1","KRT8","KRT18","KRT19","SEC14L3",
                                   "FOXJ1","MSLN","WT1","MYL7","TNNI3","TNNT2","G0S2","CMA1","HDC","GATA2","MS4A2","CCNA2","MKI67","HBB"), cols = "RdYlBu", dot.scale = 6, assay = "RNA") + RotatedAxis()

#### Pre-processing Fibroblasts ----------------------------------------------------------------------------------------------------####
Idents(HumanLung) <- "manual_clusters"
HumanLungFB <- subset(HumanLung, idents = "Fibroblast")
HumanLungFB <- FindVariableFeatures(HumanLungFB, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(HumanLungFB)
HumanLungFB <- ScaleData(HumanLungFB, verbose = T, features = all.genes)
HumanLungFB <- RunPCA(HumanLungFB, npcs = 30, verbose = FALSE)

#Calculate 90% Variance:
st_dev <- HumanLungFB@reductions$pca@stdev
var <- st_dev^2
sum(var[1:26])/ sum(var)

HumanLungFB <- FindNeighbors(object = HumanLungFB, dims = 1:26)
HumanLungFB <- FindClusters(object = HumanLungFB, resolution = 1)
HumanLungFB <- RunUMAP(HumanLungFB, reduction = "pca", dims = 1:26, verbose = F)
DimPlot(HumanLungFB, reduction = "umap", label = T)
DimPlot(HumanLungFB, reduction = "umap", label = T, split.by = "Group")

DefaultAssay(HumanLungFB) <- "RNA"
DotPlot(HumanLungFB, features = c("PTPRC","CD79A","CD19","MS4A1","CD3E","CD4","CCR7","IL2RA","FOXP3","CTLA4","NKG7","KLRB1","CCR5","CD8A","ITGAM","CD14","LY6G6C","LY6H","S100A8","ACOD1","PDGFRA","IL6","CLEC3B","IL33","COL1A1","COL1A2","PDPN","CD68","ADGRE1","C1QA","APOE",
                                  "MRC1","FCGR3A","HLA-DRB1","SIGLEC5","ITGAX","MERTK","SPP1","NR2F6","IGHA1","SDC1","ITGAE","FLT3","BATF3","XCR1","TAGLN","ACTA2","PDE5A","CLU","PECAM1","CDH5","FLT1","KIT","SNCA","SCGB1A1","LAMP3","CXCL8","AQP5","SFTPD","LYVE1","KRT8","KRT18","KRT19","SEC14L3",
                                  "FOXJ1","MSLN","WT1","MYL7","TNNI3","TNNT2","G0S2","CMA1","HDC","GATA2","MS4A2","CCNA2","MKI67","HBB"), cols = "RdYlBu", dot.scale = 6, assay = "RNA") + RotatedAxis()

#Remove non-fibroblasts 
HumanLungFB <- subset(HumanLungFB, idents = c("4"), invert = T) #remove smooth muscle

HumanLungFB <- FindVariableFeatures(HumanLungFB, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(HumanLungFB)
HumanLungFB <- ScaleData(HumanLungFB, verbose = T, features = all.genes)
HumanLungFB <- RunPCA(HumanLungFB, npcs = 30, verbose = FALSE)

#Calculate 90% Variance:
st_dev <- HumanLungFB@reductions$pca@stdev
var <- st_dev^2
sum(var[1:26])/ sum(var)

HumanLungFB <- FindNeighbors(object = HumanLungFB, dims = 1:26)
HumanLungFB <- FindClusters(object = HumanLungFB, resolution = 1)
HumanLungFB <- RunUMAP(HumanLungFB, reduction = "pca", dims = 1:26, verbose = F)
DimPlot(HumanLungFB, reduction = "umap", label = T)
DimPlot(HumanLungFB, reduction = "umap", label = T, split.by = "Group")

DefaultAssay(HumanLungFB) <- "RNA"
DotPlot(HumanLungFB, features = c("PTPRC","CD79A","CD19","MS4A1","CD3E","CD4","CCR7","IL2RA","FOXP3","CTLA4","NKG7","KLRB1","CCR5","CD8A","ITGAM","CD14","LY6G6C","LY6H","S100A8","ACOD1","PDGFRA","IL6","CLEC3B","IL33","COL1A1","COL1A2","PDPN","CD68","ADGRE1","C1QA","APOE",
                                  "MRC1","FCGR3A","HLA-DRB1","SIGLEC5","ITGAX","MERTK","SPP1","NR2F6","IGHA1","SDC1","ITGAE","FLT3","BATF3","XCR1","TAGLN","ACTA2","PDE5A","CLU","PECAM1","CDH5","FLT1","KIT","SNCA","SCGB1A1","LAMP3","CXCL8","AQP5","SFTPD","LYVE1","KRT8","KRT18","KRT19","SEC14L3",
                                  "FOXJ1","MSLN","WT1","MYL7","TNNI3","TNNT2","G0S2","CMA1","HDC","GATA2","MS4A2","CCNA2","MKI67","HBB"), cols = "RdYlBu", dot.scale = 6, assay = "RNA") + RotatedAxis()

#### Pre-processing Myeloid Cells ----------------------------------------------------------------------------------------------------####
Idents(HumanLung) <- "manual_clusters"
HumanLungMyeloid <- subset(HumanLung, idents = c("Macrophage","Monocyte","DC"))

HumanLungMyeloid <- FindVariableFeatures(HumanLungMyeloid, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(HumanLungMyeloid)
HumanLungMyeloid <- ScaleData(HumanLungMyeloid, verbose = T, features = all.genes)
HumanLungMyeloid <- RunPCA(HumanLungMyeloid, npcs = 30, verbose = FALSE)

#Calculate 90% Variance:
st_dev <- HumanLungMyeloid@reductions$pca@stdev
var <- st_dev^2
sum(var[1:22])/ sum(var)

HumanLungMyeloid <- FindNeighbors(object = HumanLungMyeloid, dims = 1:22)
HumanLungMyeloid <- FindClusters(object = HumanLungMyeloid, resolution = 1)
HumanLungMyeloid <- RunUMAP(HumanLungMyeloid, reduction = "pca", dims = 1:22, verbose = F)
DimPlot(HumanLungMyeloid, reduction = "umap", label = T)

DefaultAssay(HumanLungMyeloid) <- "RNA"

DotPlot(HumanLungMyeloid, features = c("PTPRC","CD79A","CD19","MS4A1","CD3E","CD4","CCR7","IL2RA","FOXP3","CTLA4","NKG7","KLRB1","CCR5","CD8A","ITGAM","CD14","LY6G6C","LY6H","S100A8","ACOD1","PDGFRA","IL6","CLEC3B","IL33","COL1A1","COL1A2","PDPN","CD68","ADGRE1","C1QA","APOE",
                                       "MRC1","FCGR3A","HLA-DRB1","SIGLEC5","ITGAX","MERTK","SPP1","NR2F6","IGHA1","SDC1","ITGAE","FLT3","BATF3","XCR1","TAGLN","ACTA2","PDE5A","CLU","PECAM1","CDH5","FLT1","KIT","SNCA","SCGB1A1","LAMP3","CXCL8","AQP5","SFTPD","LYVE1","KRT8","KRT18","KRT19","SEC14L3",
                                       "FOXJ1","MSLN","WT1","MYL7","TNNI3","TNNT2","G0S2","CMA1","HDC","GATA2","MS4A2","CCNA2","MKI67","HBB"), cols = "RdYlBu", dot.scale = 6, assay = "RNA") + RotatedAxis()

#remove non-myeloid cells
HumanLungMyeloid <- subset(HumanLungMyeloid, idents = c("14","17"), invert = T) #remove contamination

HumanLungMyeloid <- FindVariableFeatures(HumanLungMyeloid, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(HumanLungMyeloid)
HumanLungMyeloid <- ScaleData(HumanLungMyeloid, verbose = T, features = all.genes)
HumanLungMyeloid <- RunPCA(HumanLungMyeloid, npcs = 30, verbose = FALSE)

#Calculate 90% Variance:
st_dev <- HumanLungMyeloid@reductions$pca@stdev
var <- st_dev^2
sum(var[1:22])/ sum(var)

HumanLungMyeloid <- FindNeighbors(object = HumanLungMyeloid, dims = 1:22)
HumanLungMyeloid <- FindClusters(object = HumanLungMyeloid, resolution = 1)
HumanLungMyeloid <- RunUMAP(HumanLungMyeloid, reduction = "pca", dims = 1:22, verbose = F)
DimPlot(HumanLungMyeloid, reduction = "umap", label = T)

DefaultAssay(HumanLungMyeloid) <- "RNA"
DotPlot(HumanLungMyeloid, features = c("PTPRC","CD79A","CD19","MS4A1","CD3E","CD4","CCR7","IL2RA","FOXP3","CTLA4","NKG7","KLRB1","CCR5","CD8A","ITGAM","CD14","LY6G6C","LY6H","S100A8","ACOD1","PDGFRA","IL6","CLEC3B","IL33","COL1A1","COL1A2","PDPN","CD68","ADGRE1","C1QA","APOE",
                                       "MRC1","FCGR3A","HLA-DRB1","SIGLEC5","ITGAX","MERTK","SPP1","NR2F6","IGHA1","SDC1","ITGAE","FLT3","BATF3","XCR1","TAGLN","ACTA2","PDE5A","CLU","PECAM1","CDH5","FLT1","KIT","SNCA","SCGB1A1","LAMP3","CXCL8","AQP5","SFTPD","LYVE1","KRT8","KRT18","KRT19","SEC14L3",
                                       "FOXJ1","MSLN","WT1","MYL7","TNNI3","TNNT2","G0S2","CMA1","HDC","GATA2","MS4A2","CCNA2","MKI67","HBB"), cols = "RdYlBu", dot.scale = 6, assay = "RNA") + RotatedAxis()

#remove non-myeloid cells
HumanLungMyeloid <- subset(HumanLungMyeloid, idents = c("12","13"), invert = T) #remove contamination

HumanLungMyeloid <- FindVariableFeatures(HumanLungMyeloid, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(HumanLungMyeloid)
HumanLungMyeloid <- ScaleData(HumanLungMyeloid, verbose = T, features = all.genes)
HumanLungMyeloid <- RunPCA(HumanLungMyeloid, npcs = 30, verbose = FALSE)

#Calculate 90% Variance:
st_dev <- HumanLungMyeloid@reductions$pca@stdev
var <- st_dev^2
sum(var[1:22])/ sum(var)

HumanLungMyeloid <- FindNeighbors(object = HumanLungMyeloid, dims = 1:22)
HumanLungMyeloid <- FindClusters(object = HumanLungMyeloid, resolution = 1)
HumanLungMyeloid <- RunUMAP(HumanLungMyeloid, reduction = "pca", dims = 1:22, verbose = F)
DimPlot(HumanLungMyeloid, reduction = "umap", label = T)

DefaultAssay(HumanLungMyeloid) <- "RNA"
DotPlot(HumanLungMyeloid, features = c("PTPRC","CD79A","CD19","MS4A1","CD3E","CD4","CCR7","IL2RA","FOXP3","CTLA4","NKG7","KLRB1","CCR5","CD8A","ITGAM","CD14","LY6G6C","LY6H","S100A8","ACOD1","PDGFRA","IL6","CLEC3B","IL33","COL1A1","COL1A2","PDPN","CD68","ADGRE1","C1QA","APOE",
                                       "MRC1","FCGR3A","HLA-DRB1","SIGLEC5","ITGAX","MERTK","SPP1","NR2F6","IGHA1","SDC1","ITGAE","FLT3","BATF3","XCR1","TAGLN","ACTA2","PDE5A","CLU","PECAM1","CDH5","FLT1","KIT","SNCA","SCGB1A1","LAMP3","CXCL8","AQP5","SFTPD","LYVE1","KRT8","KRT18","KRT19","SEC14L3",
                                       "FOXJ1","MSLN","WT1","MYL7","TNNI3","TNNT2","G0S2","CMA1","HDC","GATA2","MS4A2","CCNA2","MKI67","HBB"), cols = "RdYlBu", dot.scale = 6, assay = "RNA") + RotatedAxis()

#cDC1 markers
DotPlot(HumanLungMyeloid, features = c("BATF3", "IRF8", "CLEC9A", "ITGAE", "ITGAX", "CADM1","CD8A", "THBD", "XCR1"), cols = c("RdYlBu"), dot.scale = 8, assay = "RNA") + RotatedAxis()
#active DC markers
DotPlot(HumanLungMyeloid, features = c("CCL22", "CCL17", "CCR7"), cols = c("RdYlBu"), dot.scale = 8, assay = "RNA") + RotatedAxis()
#Macrophages
DotPlot(HumanLungMyeloid, features = c("PTPRC","ITGAM","CD14","LY6G6C","LY6H","S100A8","ACOD1","CD68","ADGRE1","C1QA","APOE","MRC1","FCGR3A","HLA-DRB1",
                                       "SIGLEC5","ITGAX","MERTK","SPP1","G0S2","CCNA2","MKI67"), cols = "RdYlBu", dot.scale = 6, assay = "RNA") + RotatedAxis()
#monocyte
DotPlot(HumanLungMyeloid, features = c("LY6H", "CCR2",  "ITGAM",  "CX3CR1", "CD14", "FCGR3A"), cols = c("RdYlBu"), dot.scale = 8, assay = "RNA") + RotatedAxis()
#classical monocyte (FCGR3A aka Cd16 negative)
DotPlot(HumanLungMyeloid, features = c("CD14", "FCGR3A", "CCR2",  "CCR5",  "SELL"), cols = c("RdYlBu"), dot.scale = 8, assay = "RNA") + RotatedAxis() 

Idents(HumanLungMyeloid) <- "seurat_clusters"
HumanLungMyeloid <- RenameIdents(HumanLungMyeloid,
                                 "0" = "Macrophage", 
                                 "1" = "Macrophage", 
                                 "2" = "Macrophage", 
                                 "3" = "Macrophage", 
                                 "4" = "Monocyte", 
                                 "5" = "Classical Monocyte", 
                                 "6" = "Macrophage", 
                                 "7" = "Classical Monocyte", 
                                 "8" = "Macrophage",
                                 "9" = "Macrophage", 
                                 "10" = "Macrophage", 
                                 "11" = "Macrophage", 
                                 "12" = "Macrophage", 
                                 "13" = "cDC1", 
                                 "14" = "Active DC",
                                 "15" = "Macrophage")
HumanLungMyeloid[["manual_clusters"]] <- HumanLungMyeloid@active.ident

DimPlot(HumanLungMyeloid, label = T, repel = T)+NoLegend()

HumanLungMacrophage <- subset(HumanLungMyeloid, idents = "Macrophage")


#### Figures ----------------------------------------------------------------------------------------------------####


### Figure 7A ###
Idents(HumanLung) <- "manual_clusters"

DimPlot(HumanLung,cols = c("#DF4945","#D7A8FF","#52E7ED","#528AFF",
                           "#FFE18F","#73D279","#FFB300","#7FA4CA",
                           "#3A9255","#009193","#FF9360","#C1D274"))+NoLegend()

### Figure 7C ###
#G12D vs wt
HumanLungEpithelial_DE <- FindMarkers(HumanLungEpithelial, ident.1 = "KRAS-G12D", ident.2 = "KRAS-wt", group.by = "Group")

genes <- c("IL37","IL4R","IL6ST","IL32","IL13RA1","IL10","IL18","IL17RE","CXCL3",
                "CXCL8","CXCL17","CXCL16","CXCL14","CXCL2","CXCL1","CCL28","CCL20",
                "TGFBR1","TGFA","IGF2R","EGFR","HDGF","VEGFA","FGFR1","HBEGF")

data <- FetchData(HumanLungEpithelial, vars = c(genes, "Group"))

df_avg <- data.frame()

for (id in levels(factor(data$Group))) {
  data_subset <- data %>% filter(Group == id)
  data_subset_avg <- apply(data_subset[,1:(ncol(data_subset)-1)], 2, mean)
  df_avg <- rbind(df_avg, data_subset_avg)
}

off_score <- c()
on_score <- c()

for (i in 1:(ncol(df_avg))) {
  if (df_avg[1,i] > df_avg[2,i]){
    ratio_1 <- df_avg[2,i]/df_avg[1,i]
    off_score <- append(off_score, ratio_1)
    on_score <- append(on_score, 1)
  } else {
    ratio_2 <- df_avg[1,i]/df_avg[2,i]
    on_score <- append(on_score, ratio_2)
    off_score <- append(off_score, 1)
  }
  df_avg_ratios <- rbind(df_avg, on_score, off_score)
  df_avg_ratios <- df_avg_ratios[3:4,]
}

colnames(df_avg_ratios) <- colnames(data)[1:(ncol(data)-1)]
rownames(df_avg_ratios) <- levels(factor(data$Group))
pheatmap(as.matrix(df_avg_ratios), fontsize = 14, color = brewer.pal(9, "OrRd"), scale = 'none', cluster_rows = F, cluster_cols = T)

### Figure 7D ###

VlnPlot(HumanLungEpithelial, features = "CXCL1", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(HumanLungEpithelial, features = "CXCL2", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(HumanLungEpithelial, features = "HDGF", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

### Figure 7E ###
#G12D vs wt
HumanLungMacrophage_DE <- FindMarkers(HumanLungMacrophage, ident.1 = "KRAS-G12D", ident.2 = "KRAS-wt", group.by = "Group")

genes <- c("IL17RB","IL6","IL3RA","IL6R","IL32","IL10RA","IL4I1","IL1RN","CXCL9","CXCL8",
           "CXCL2","CXCL1","CXCL3","CXCL10","CCL18","CCL3","CCL4","CCL20","CCL2","S100A8",
           "S100A6","S100A4","MRC1","APOE","FN1","STAT1","C1QB","C1QC","C1QA","CCR1","CCR7",
           "CD83","TLR2","CSF1R")

data <- FetchData(HumanLungMacrophage, vars = c(genes, "Group"))

df_avg <- data.frame()

for (id in levels(factor(data$Group))) {
  data_subset <- data %>% filter(Group == id)
  data_subset_avg <- apply(data_subset[,1:(ncol(data_subset)-1)], 2, mean)
  df_avg <- rbind(df_avg, data_subset_avg)
}

off_score <- c()
on_score <- c()

for (i in 1:(ncol(df_avg))) {
  if (df_avg[1,i] > df_avg[2,i]){
    ratio_1 <- df_avg[2,i]/df_avg[1,i]
    off_score <- append(off_score, ratio_1)
    on_score <- append(on_score, 1)
  } else {
    ratio_2 <- df_avg[1,i]/df_avg[2,i]
    on_score <- append(on_score, ratio_2)
    off_score <- append(off_score, 1)
  }
  df_avg_ratios <- rbind(df_avg, on_score, off_score)
  df_avg_ratios <- df_avg_ratios[3:4,]
}

colnames(df_avg_ratios) <- colnames(data)[1:(ncol(data)-1)]
rownames(df_avg_ratios) <- levels(factor(data$Group))
pheatmap(as.matrix(df_avg_ratios), fontsize = 14, color = brewer.pal(9, "YlGnBu"), scale = 'none', cluster_rows = F, cluster_cols = T)

### Figure 7F ###
Idents(HumanLungMacrophage) <- "Group"

VlnPlot(HumanLungMacrophage, features = "C1QB", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(HumanLungMacrophage, features = "C1QC", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

### Supplemental Figure 7A ###

#G12D vs wt
HumanLungFB_DE <- FindMarkers(HumanLungFB, ident.1 = "KRAS-G12D", ident.2 = "KRAS-wt", group.by = "Group")

genes <- c("IL13RA1","IL17RA","IL32","IL6","IL1B","CXCL2","CXCL12","CXCL1","CXCL14",
                "CCL2","CCL28","CCL21","CCL19","CXCR4","TGFBI","TGFBR1","TGFBR2","TGFBR3",
                "TGFB1I1","CTGF","FGFR1","PDGFRB","IGF1","HGF","PDGFRA","FGF13","VEGFA",
                "PDGFA","MMP14","MMP11","TIMP2","TIMP3","TIMP1","TNFAIP6","IFNGR1","IFNAR2",
                "FBLN5","FBLN1","COL6A1","COL4A1","COL8A1","JAK3","STAT1","STAT3")

data <- FetchData(HumanLungFB, vars = c(genes, "Group"))

df_avg <- data.frame()

for (id in levels(factor(data$Group))) {
  data_subset <- data %>% filter(Group == id)
  data_subset_avg <- apply(data_subset[,1:(ncol(data_subset)-1)], 2, mean)
  df_avg <- rbind(df_avg, data_subset_avg)
}

off_score <- c()
on_score <- c()

for (i in 1:(ncol(df_avg))) {
  if (df_avg[1,i] > df_avg[2,i]){
    ratio_1 <- df_avg[2,i]/df_avg[1,i]
    off_score <- append(off_score, ratio_1)
    on_score <- append(on_score, 1)
  } else {
    ratio_2 <- df_avg[1,i]/df_avg[2,i]
    on_score <- append(on_score, ratio_2)
    off_score <- append(off_score, 1)
  }
  df_avg_ratios <- rbind(df_avg, on_score, off_score)
  df_avg_ratios <- df_avg_ratios[3:4,]
}

colnames(df_avg_ratios) <- colnames(data)[1:(ncol(data)-1)]
rownames(df_avg_ratios) <- levels(factor(data$Group))
pheatmap(as.matrix(df_avg_ratios), fontsize = 14, color = brewer.pal(9, "BuPu"), scale = 'none', cluster_rows = F, cluster_cols = T)
