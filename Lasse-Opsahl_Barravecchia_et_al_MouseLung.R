#This script recreates Figures 2A-B, 2D-G, S2A-C, 3D-F, S3B-E, 4C-G, S4E.
#Raw data files for the novel dataset generated in this manuscript are available through the NIH Gene Expression Omnibus (GEO), accession number GSE281744.
#L-iKras 21 weeks ON (GSM8627458)
#L-iKras 20 weeks ON, 1 week OFF (GSM8627459)

#Data was processed in line with the Seurat workflow:
#Website: https://satijalab.org/seurat/index.html

#Reference: Hao et al., Integrated analysis of multimodal single-cell data. 
#Cell. 2021 Jun 24;184(13):3573-3587.e29. doi: 10.1016/j.cell.2021.04.048. 
#PMID: 34062119; PMCID: PMC8238499. 

#Seurat Version 4.0.2
#SeuratObject Version 4.0.1
#R version 4.2.3 (2023-03-15) --- Shortstop Beagle

#Load packages
library(Seurat)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(scater)
library(pheatmap)
library(enrichplot)
library(dplyr)
library(msigdbr)
library(fgsea)

#### Pre-processing ----------------------------------------------------------------------------------------------------####

#Load Raw Data:
iKras_ON_data <- Read10X_h5("~/Lasse-Opsahl_Barravecchia_et_al_GSM8627458/filtered_feature_bc_matrix.h5")
iKras_OFF_data <- Read10X_h5("~/Lasse-Opsahl_Barravecchia_et_al_GSM8627459/filtered_feature_bc_matrix.h5")

#Make Seurat Objects (include features/genes expressed by at least 3 cells and include cells with at least 100 features/genes)
iKras_ON <- CreateSeuratObject(iKras_ON_data, min.cells = 3, min.features = 100)
iKras_OFF <- CreateSeuratObject(iKras_OFF_data, min.cells = 3, min.features = 100)

#Add Metadata:
iKras_ON[["Group"]] <- "Kras ON"
iKras_OFF[["Group"]] <- "Kras OFF"

iKras_ON[["Run_Date"]] <- "APR_2021"
iKras_OFF[["Run_Date"]] <- "APR_2021"

iKras_ON[["Time"]] <- "KrasON_21wk"
iKras_OFF[["Time"]] <- "KrasON_20wk_OFF_1wk"

iKras_ON[["Sample"]] <- "KrasON_21wk_lung"
iKras_OFF[["Sample"]] <- "KrasON_20wk_OFF_1wk_lung"

#Merge Samples:
MouseLung <- merge(x = iKras_ON, y = c(iKras_OFF), add.cell.ids = (c("ON","OFF")))

#Log Normalize:
MouseLung <- NormalizeData(object = MouseLung, normalization.method = "LogNormalize", scale.factor = 10000)

#Check percent mitochondrial genes to check for doublets and poor cell quality
MouseLung[["percent.mt"]] <- PercentageFeatureSet(object = MouseLung, pattern = "^mt-")

#Keep cells with between 800 and 100,000 RNA reads and less than 15% mitochondrial genes
MouseLung <- subset(x = MouseLung, subset = nCount_RNA > 800 & nCount_RNA < 100000 & percent.mt < 15)

#Identify variable genes
MouseLung <- FindVariableFeatures(MouseLung, selection.method = "vst", nfeatures = 2000)

#Scale data and run PCA dimensionality reduction
all.genes <- rownames(MouseLung)
MouseLung <- ScaleData(MouseLung, verbose = T, features = all.genes)
MouseLung <- RunPCA(MouseLung, npcs = 30, verbose = FALSE)

#Calculate 90% Variance:
st_dev <- MouseLung@reductions$pca@stdev
var <- st_dev^2
sum(var[1:15])/ sum(var) 

#Find Neighbors and cluster cells
MouseLung <- FindNeighbors(object = MouseLung, dims = 1:15)
MouseLung <- FindClusters(object = MouseLung, resolution = 1)
MouseLung <- RunUMAP(MouseLung, reduction = "pca", dims = 1:15, verbose = F)
DimPlot(MouseLung, reduction = "umap", label = T)

#Identify clusters based on marker expression
DotPlot(MouseLung, features=c("Ptprc", "Cd79a", "Cd19", "Ms4a1","Cd3e", "Cd4","Ccr7","Il2ra", "Foxp3","Ctla4", "Nkg7","Klrb1c","Ccr5", "Cd8a", "Itgam","Cd14","Ly6g", "Ly6c2", "S100a8","Acod1","Pdgfra","Il6", "Clec3b", "Il33", "Col1a1","Col1a2","Pdpn","Cd68","Adgre1","C1qa","Apoe", 
                           "Mrc1","Fcgr3","H2-Eb1","Siglecf","Itgax","Mertk","Spp1","Ear2","Igha", "Sdc1","Itgae","Flt3", "Batf3","Xcr1","Tagln","Acta2","Pde5a","Clu","Pecam1", "Cdh5","Flt1","Kit", "Snca","Scgb1a1","Lamp3","Cxcl15","Aqp5", "Sftpd","Lyve1","Krt8", "Krt18", "Krt19","Sec14l3",
                           "Foxj1","Msln", "Wt1", "Myl7", "Tnni3","Tnnt2","G0s2","Mcpt4","Mcpt8","Gata2", "Ms4a2","Ccna2", "Mki67","Hbb-bt"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

#Label Manual Clusters:
Idents(MouseLung) <- "seurat_clusters"
MouseLung <- RenameIdents(MouseLung,
                       "0" = "T Cell", 
                       "1" = "T Cell", 
                       "2" = "RBC", 
                       "3" = "RBC", 
                       "4" = "B Cell", 
                       "5" = "RBC", 
                       "6" = "B Cell", 
                       "7" = "NK Cell", 
                       "8" = "Neutrophil",
                       "9" = "Fibroblast", 
                       "10" = "Macrophage", 
                       "11" = "Macrophage", 
                       "12" = "Eosinophil", 
                       "13" = "Fibroblast", 
                       "14" = "Dendritic Cell", 
                       "15" = "T Cell", 
                       "16" = "Pericyte", 
                       "17" = "Endothelial",
                       "18" = "Plasma Cell", 
                       "19" = "Epithelial", 
                       "20" = "Smooth Muscle", 
                       "21" = "Endothelial",
                       "22" = "Epithelial", 
                       "23" = "Cycling T Cell", 
                       "24" = "Mesothelial", 
                       "25" = "Cardiomyocyte", 
                       "26" = "Endothelial")
MouseLung[["manual_clusters"]] <- MouseLung@active.ident
Idents(MouseLung) <- "manual_clusters"

new_order <- c("Epithelial",
               "Fibroblast",
               "Neutrophil",
               "Macrophage",
               "T Cell",
               "Cycling T Cell",
               "NK Cell",
               "Eosinophil",
               "Mesothelial",
               "Pericyte",
               "Dendritic Cell",
               "Smooth Muscle",
               "Endothelial",
               "B Cell",
               "Plasma Cell",
               "Cardiomyocyte",
               "RBC")

MouseLung@active.ident <- factor(MouseLung@active.ident, levels = new_order)
MouseLung[["manual_clusters"]] <- MouseLung@active.ident

Idents(MouseLung) <- "manual_clusters"
MouseLung$Group <- factor(MouseLung$Group, levels = c("Kras ON","Kras OFF"))

#### Pre-processing Epithelial Cells ----------------------------------------------------------------------------------------------------####

Idents(MouseLung) <- "manual_clusters"
MouseEpithelial <- subset(MouseLung, idents = c("Epithelial"))

MouseEpithelial <- FindVariableFeatures(MouseEpithelial, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(MouseEpithelial)
MouseEpithelial <- ScaleData(MouseEpithelial, verbose = T, features = all.genes)
MouseEpithelial <- RunPCA(MouseEpithelial, npcs = 30, verbose = FALSE)

#Calculate 90% Variance:
st_dev <- MouseEpithelial@reductions$pca@stdev
var <- st_dev^2
sum(var[1:21])/ sum(var)

MouseEpithelial <- FindNeighbors(object = MouseEpithelial, dims = 1:21)
MouseEpithelial <- FindClusters(object = MouseEpithelial, resolution = 2)
MouseEpithelial <- RunUMAP(MouseEpithelial, reduction = "pca", dims = 1:21, verbose = F)
DimPlot(MouseEpithelial, reduction = "umap", label = T)

DotPlot(MouseEpithelial, features=c("Ptprc", "Cd79a", "Cd19", "Ms4a1","Cd3e", "Cd4","Ccr7","Il2ra", "Foxp3","Ctla4", "Nkg7","Klrb1c","Ccr5", "Cd8a", "Itgam","Cd14","Ly6g", "Ly6c2", "S100a8","Acod1","Pdgfra","Il6", "Clec3b", "Il33", "Col1a1","Col1a2","Pdpn","Cd68","Adgre1","C1qa","Apoe", 
                               "Mrc1","Fcgr3","H2-Eb1","Siglecf","Itgax","Mertk","Spp1","Ear2","Igha", "Sdc1","Itgae","Flt3", "Batf3","Xcr1","Tagln","Acta2","Pde5a","Clu","Pecam1", "Cdh5","Flt1","Kit", "Snca","Scgb1a1","Lamp3","Cxcl15","Aqp5", "Sftpd","Lyve1","Krt8", "Krt18", "Krt19","Sec14l3",
                               "Foxj1","Msln", "Wt1", "Myl7", "Tnni3","Tnnt2","G0s2","Mcpt4","Mcpt8","Gata2", "Ms4a2","Ccna2", "Mki67","Hbb-bt"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

#Remove non-epithelial cells
MouseEpithelial <- subset(MouseEpithelial, idents = c("6","9","11"), invert = T) #remove contamination

MouseEpithelial <- FindVariableFeatures(MouseEpithelial, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(MouseEpithelial)
MouseEpithelial <- ScaleData(MouseEpithelial, verbose = T, features = all.genes)
MouseEpithelial <- RunPCA(MouseEpithelial, npcs = 30, verbose = FALSE)

#Calculate 90% Variance:
st_dev <- MouseEpithelial@reductions$pca@stdev
var <- st_dev^2
sum(var[1:21])/ sum(var)

MouseEpithelial <- FindNeighbors(object = MouseEpithelial, dims = 1:21)
MouseEpithelial <- FindClusters(object = MouseEpithelial, resolution = 1)
MouseEpithelial <- RunUMAP(MouseEpithelial, reduction = "pca", dims = 1:21, verbose = F)
DimPlot(MouseEpithelial, reduction = "umap", label = T)

DotPlot(MouseEpithelial, features=c("Ptprc", "Cd79a", "Cd19", "Ms4a1","Cd3e", "Cd4","Ccr7","Il2ra", "Foxp3","Ctla4", "Nkg7","Klrb1c","Ccr5", "Cd8a", "Itgam","Cd14","Ly6g", "Ly6c2", "S100a8","Acod1","Pdgfra","Il6", "Clec3b", "Il33", "Col1a1","Col1a2","Pdpn","Cd68","Adgre1","C1qa","Apoe", 
                               "Mrc1","Fcgr3","H2-Eb1","Siglecf","Itgax","Mertk","Spp1","Ear2","Igha", "Sdc1","Itgae","Flt3", "Batf3","Xcr1","Tagln","Acta2","Pde5a","Clu","Pecam1", "Cdh5","Flt1","Kit", "Snca","Scgb1a1","Lamp3","Cxcl15","Aqp5", "Sftpd","Lyve1","Krt8", "Krt18", "Krt19","Sec14l3",
                               "Foxj1","Msln", "Wt1", "Myl7", "Tnni3","Tnnt2","G0s2","Mcpt4","Mcpt8","Gata2", "Ms4a2","Ccna2", "Mki67","Hbb-bt"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

#Remove non-epithelial cells
MouseEpithelial <- subset(MouseEpithelial, idents = c("6"), invert = T) #remove contamination

MouseEpithelial <- FindVariableFeatures(MouseEpithelial, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(MouseEpithelial)
MouseEpithelial <- ScaleData(MouseEpithelial, verbose = T, features = all.genes)
MouseEpithelial <- RunPCA(MouseEpithelial, npcs = 30, verbose = FALSE)

#Calculate 90% Variance:
st_dev <- MouseEpithelial@reductions$pca@stdev
var <- st_dev^2
sum(var[1:20])/ sum(var)

MouseEpithelial <- FindNeighbors(object = MouseEpithelial, dims = 1:20)
MouseEpithelial <- FindClusters(object = MouseEpithelial, resolution = 1)
MouseEpithelial <- RunUMAP(MouseEpithelial, reduction = "pca", dims = 1:20, verbose = F)
DimPlot(MouseEpithelial, reduction = "umap", label = T)

# bronchial epithelium
DotPlot(MouseEpithelial, features=c("Krt5", "Tp63", "Muc5b", "Muc5ac", "Foxj1", "Actub", "Ccsp", "Scgb1a1", "Scgb3a2", "Pou2f3", "Syp", "Chga","Pgp9.5","Robo2","Eno2","Foxi1","Cftr", "Sec14l3"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()
#Club Cell
DotPlot(MouseEpithelial, features= c("Scgb1a1", "Scgb3a1","Scgb3a2","Tuba1a"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()
#AT1
DotPlot(MouseEpithelial, features= c("Cldn18", "Gprc5a", "Lmo7", "Ager"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()
#AT2
DotPlot(MouseEpithelial, features= c("Sftpc", "Sftpb", "Ctsh","Dram1", "Lamp3"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

#Label Manual Clusters:
Idents(MouseEpithelial) <- "seurat_clusters"
MouseEpithelial <- RenameIdents(MouseEpithelial,
                           "0" = "Bronchial Epithelial", 
                           "1" = "AT2", 
                           "2" = "AT2", 
                           "3" = "Bronchial Epithelial",
                           "4" = "Clara Cell",
                           "5" = "AT1")
MouseEpithelial[["manual_clusters"]] <- MouseEpithelial@active.ident

#### Pre-processing Fibroblasts and Pericytes ----------------------------------------------------------------------------------------------------####

Idents(MouseLung) <- "manual_clusters"
MouseLungFB <- subset(MouseLung, idents = c("Fibroblast","Pericyte"))

MouseLungFB <- FindVariableFeatures(MouseLungFB, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(MouseLungFB)
MouseLungFB <- ScaleData(MouseLungFB, verbose = T, features = all.genes)
MouseLungFB <- RunPCA(MouseLungFB, npcs = 30, verbose = FALSE)

#Calculate 90% Variance:
st_dev <- MouseLungFB@reductions$pca@stdev
var <- st_dev^2
sum(var[1:22])/ sum(var)

MouseLungFB <- FindNeighbors(object = MouseLungFB, dims = 1:22)
MouseLungFB <- FindClusters(object = MouseLungFB, resolution = 1)
MouseLungFB <- RunUMAP(MouseLungFB, reduction = "pca", dims = 1:22, verbose = F)
DimPlot(MouseLungFB, reduction = "umap", label = T)

DotPlot(MouseLungFB, features=c("Ptprc", "Cd79a", "Cd19", "Ms4a1","Cd3e", "Cd4","Ccr7","Il2ra", "Foxp3","Ctla4", "Nkg7","Klrb1c","Ccr5", "Cd8a", "Itgam","Cd14","Ly6g", "Ly6c2", "S100a8","Acod1","Pdgfra","Il6", "Clec3b", "Il33", "Col1a1","Col1a2","Pdpn","Cd68","Adgre1","C1qa","Apoe", 
                                "Mrc1","Fcgr3","H2-Eb1","Siglecf","Itgax","Mertk","Spp1","Ear2","Igha", "Sdc1","Itgae","Flt3", "Batf3","Xcr1","Tagln","Acta2","Pde5a","Clu","Pecam1", "Cdh5","Flt1","Kit", "Snca","Scgb1a1","Lamp3","Cxcl15","Aqp5", "Sftpd","Lyve1","Krt8", "Krt18", "Krt19","Sec14l3",
                                "Foxj1","Msln", "Wt1", "Myl7", "Tnni3","Tnnt2","G0s2","Mcpt4","Mcpt8","Gata2", "Ms4a2","Ccna2", "Mki67","Hbb-bt"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

#remove cells that aren't fibroblasts or pericytes
MouseLungFB <- subset(MouseLungFB, idents = c("6","8","10","11"), invert = T) #remove contamination

MouseLungFB <- FindVariableFeatures(MouseLungFB, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(MouseLungFB)
MouseLungFB <- ScaleData(MouseLungFB, verbose = T, features = all.genes)
MouseLungFB <- RunPCA(MouseLungFB, npcs = 30, verbose = FALSE)

#Calculate 90% Variance:
st_dev <- MouseLungFB@reductions$pca@stdev
var <- st_dev^2
sum(var[1:22])/ sum(var)

MouseLungFB <- FindNeighbors(object = MouseLungFB, dims = 1:22)
MouseLungFB <- FindClusters(object = MouseLungFB, resolution = 1)
MouseLungFB <- RunUMAP(MouseLungFB, reduction = "pca", dims = 1:22, verbose = F)
DimPlot(MouseLungFB, reduction = "umap", label = T)
DimPlot(MouseLungFB, reduction = "umap", label = T, split.by = "Group")

#remove cells that aren't fibroblasts or pericytes
MouseLungFB <- subset(MouseLungFB, idents = c("3","8","11"), invert = T) #remove contamination

MouseLungFB <- FindVariableFeatures(MouseLungFB, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(MouseLungFB)
MouseLungFB <- ScaleData(MouseLungFB, verbose = T, features = all.genes)
MouseLungFB <- RunPCA(MouseLungFB, npcs = 30, verbose = FALSE)

#Calculate 90% Variance:
st_dev <- MouseLungFB@reductions$pca@stdev
var <- st_dev^2
sum(var[1:22])/ sum(var)

MouseLungFB <- FindNeighbors(object = MouseLungFB, dims = 1:22)
MouseLungFB <- FindClusters(object = MouseLungFB, resolution = 1)
MouseLungFB <- RunUMAP(MouseLungFB, reduction = "pca", dims = 1:22, verbose = F)
DimPlot(MouseLungFB, reduction = "umap", label = T)

DotPlot(MouseLungFB, features=c("Ptprc", "Cd79a", "Cd19", "Ms4a1","Cd3e", "Cd4","Ccr7","Il2ra", "Foxp3","Ctla4", "Nkg7","Klrb1c","Ccr5", "Cd8a", "Itgam","Cd14","Ly6g", "Ly6c2", "S100a8","Acod1","Pdgfra","Il6", "Clec3b", "Il33", "Col1a1","Col1a2","Pdpn","Cd68","Adgre1","C1qa","Apoe", 
                                "Mrc1","Fcgr3","H2-Eb1","Siglecf","Itgax","Mertk","Spp1","Ear2","Igha", "Sdc1","Itgae","Flt3", "Batf3","Xcr1","Tagln","Acta2","Pde5a","Clu","Pecam1", "Cdh5","Flt1","Kit", "Snca","Scgb1a1","Lamp3","Cxcl15","Aqp5", "Sftpd","Lyve1","Krt8", "Krt18", "Krt19","Sec14l3",
                                "Foxj1","Msln", "Wt1", "Myl7", "Tnni3","Tnnt2","G0s2","Mcpt4","Mcpt8","Gata2", "Ms4a2","Ccna2", "Mki67","Hbb-bt"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

#Lipofibroblast markers
DotPlot(MouseLungFB, features = c("Tcf21","Plin2","Fgf10","G0s2","Lpl","Ear1", "Ear2", "Mrc1", "Clec4n"))
#Myofibroblasts
DotPlot(MouseLungFB, features = c("Tagln","Fap","Mmp11","Myl9","Hopx","Postn","Tpm1","Tpm2","Thy1","Col12a1","Thbs2","Tnc","Tgfb1","Pdgfrb",
                                  "Hhip","Aspn","Mustn1","Enpp2","Igfbp3","Grem2","Lum","Bmp5","Acta2","Gja1","Cyp2e1","Des"), cols = "RdYlBu")+RotatedAxis()
#Col13a1 Matrix FB
DotPlot(MouseLungFB, features = c("Col13a1","Itga8","Cxcl14","Npnt","Hsd11b1","Tcf21","Mfap4","Spon1","Limch1","Cdh11","Slc7a10","Gm14964","Scn3a","Nebl","Mrc2","Bmp3"), cols = "RdYlBu")+RotatedAxis()
#Col14a1 Matrix FB
DotPlot(MouseLungFB, features = c("Pi16","Mmp3","Clec3b","Cygb","Dcn","Rbp4","Gsn","Col14a1","Dpep1","Fbln1","Htra3","Igfbp4","Cxcl12","Gas7"), cols = "RdYlBu")+RotatedAxis()
#FB-like
DotPlot(MouseLungFB, features = c("Cyp1b1","Apod","Pdgfrb"), cols = "RdYlBu")
#Pericytes
DotPlot(MouseLungFB, features = c("Rgs5","Pdgfrb"), cols = "RdYlBu")

#Label Manual Clusters:
Idents(MouseLungFB) <- "seurat_clusters"
MouseLungFB <- RenameIdents(MouseLungFB,
                            "0" = "COL14A1 + Matrix FB", 
                            "1" = "COL14A1 + Matrix FB", 
                            "2" = "COL14A1 + Matrix FB", 
                            "3" = "Pericyte", 
                            "4" = "Lipofibroblast", 
                            "5" = "COL13A1 + Matrix FB", 
                            "6" = "COL13A1 + Matrix FB", 
                            "7" = "COL14A1 + Matrix FB", 
                            "8" = "FB-like",
                            "9" = "Myofibroblast")
MouseLungFB[["manual_clusters"]] <- MouseLungFB@active.ident

new_order <- c("Myofibroblast",
               "Lipofibroblast",
               "COL13A1 + Matrix FB",
               "COL14A1 + Matrix FB",
               "FB-like",
               "Pericyte")

#for dotplot order
new_order <- c("Pericyte",
               "FB-like",
               "COL14A1 + Matrix FB",
               "COL13A1 + Matrix FB",
               "Lipofibroblast",
               "Myofibroblast")

MouseLungFB@active.ident <- factor(MouseLungFB@active.ident, levels = new_order)
MouseLungFB[["manual_clusters"]] <- MouseLungFB@active.ident

MouseLungFB_noPericytes <- subset(MouseLungFB, idents = "Pericyte", invert = T)

#### Pre-processing Myeloid Cells ----------------------------------------------------------------------------------------------------####
Idents(MouseLung) <- "manual_clusters"
MouseLungMyeloid <- subset(MouseLung, idents = c("Macrophage", "Neutrophil","Eosinophil","Dendritic Cell"))

MouseLungMyeloid <- FindVariableFeatures(MouseLungMyeloid, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(MouseLungMyeloid)
MouseLungMyeloid <- ScaleData(MouseLungMyeloid, verbose = T, features = all.genes)
MouseLungMyeloid <- RunPCA(MouseLungMyeloid, npcs = 30, verbose = FALSE)

#Calculate 90% Variance:
st_dev <- MouseLungMyeloid@reductions$pca@stdev
var <- st_dev^2
sum(var[1:20])/ sum(var)

MouseLungMyeloid <- FindNeighbors(object = MouseLungMyeloid, dims = 1:20)
MouseLungMyeloid <- FindClusters(object = MouseLungMyeloid, resolution = 4)
MouseLungMyeloid <- RunUMAP(MouseLungMyeloid, reduction = "pca", dims = 1:20, verbose = F)
DimPlot(MouseLungMyeloid, reduction = "umap", label = T)
DimPlot(MouseLungMyeloid, reduction = "umap", label = T, group.by = "manual_clusters", split.by = "Group")

DotPlot(MouseLungMyeloid, features=c("Ptprc", "Cd79a", "Cd19", "Ms4a1","Cd3e", "Cd4","Ccr7","Il2ra", "Foxp3","Ctla4", "Nkg7","Klrb1c","Ccr5", "Cd8a", "Itgam","Cd14","Ly6g", "Ly6c2", "S100a8","Acod1","Pdgfra","Il6", "Clec3b", "Il33", "Col1a1","Col1a2","Pdpn","Cd68","Adgre1","C1qa","Apoe", 
                                "Mrc1","Fcgr3","H2-Eb1","Siglecf","Itgax","Mertk","Spp1","Ear2","Igha", "Sdc1","Itgae","Flt3", "Batf3","Xcr1","Tagln","Acta2","Pde5a","Clu","Pecam1", "Cdh5","Flt1","Kit", "Snca","Scgb1a1","Lamp3","Cxcl15","Aqp5", "Sftpd","Lyve1","Krt8", "Krt18", "Krt19","Sec14l3",
                                "Foxj1","Msln", "Wt1", "Myl7", "Tnni3","Tnnt2","G0s2","Mcpt4","Mcpt8","Gata2", "Ms4a2","Ccna2", "Mki67","Hbb-bt"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

#Remove non-myeloid cells
MouseLungMyeloid <- subset(MouseLungMyeloid, idents = c("2","7","12","14","15","16","18","20","21","22","26","27","28","31","32","33"), invert = T) #remove contamination

MouseLungMyeloid <- FindVariableFeatures(MouseLungMyeloid, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(MouseLungMyeloid)
MouseLungMyeloid <- ScaleData(MouseLungMyeloid, verbose = T, features = all.genes)
MouseLungMyeloid <- RunPCA(MouseLungMyeloid, npcs = 30, verbose = FALSE)

#Calculate 90% Variance:
st_dev <- MouseLungMyeloid@reductions$pca@stdev
var <- st_dev^2
sum(var[1:19])/ sum(var)

MouseLungMyeloid <- FindNeighbors(object = MouseLungMyeloid, dims = 1:19)
MouseLungMyeloid <- FindClusters(object = MouseLungMyeloid, resolution = 2)
MouseLungMyeloid <- RunUMAP(MouseLungMyeloid, reduction = "pca", dims = 1:19, verbose = F)
DimPlot(MouseLungMyeloid, reduction = "umap", label = T)
DimPlot(MouseLungMyeloid, reduction = "umap", label = T, group.by = "manual_clusters", split.by = "Group")

DotPlot(MouseLungMyeloid, features=c("Ptprc", "Cd79a", "Cd19", "Ms4a1","Cd3e", "Cd4","Ccr7","Il2ra", "Foxp3","Ctla4", "Nkg7","Klrb1c","Ccr5", "Cd8a", "Itgam","Cd14","Ly6g", "Ly6c2", "S100a8","Acod1","Pdgfra","Il6", "Clec3b", "Il33", "Col1a1","Col1a2","Pdpn","Cd68","Adgre1","C1qa","Apoe", 
                                     "Mrc1","Fcgr3","H2-Eb1","Siglecf","Itgax","Mertk","Spp1","Ear2","Igha", "Sdc1","Itgae","Flt3", "Batf3","Xcr1","Tagln","Acta2","Pde5a","Clu","Pecam1", "Cdh5","Flt1","Kit", "Snca","Scgb1a1","Lamp3","Cxcl15","Aqp5", "Sftpd","Lyve1","Krt8", "Krt18", "Krt19","Sec14l3",
                                     "Foxj1","Msln", "Wt1", "Myl7", "Tnni3","Tnnt2","G0s2","Mcpt4","Mcpt8","Gata2", "Ms4a2","Ccna2", "Mki67","Hbb-bt"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

#Remove non-myeloid cells
MouseLungMyeloid <- subset(MouseLungMyeloid, idents = c("7","12","16","17"), invert = T) #remove contamination

MouseLungMyeloid <- FindVariableFeatures(MouseLungMyeloid, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(MouseLungMyeloid)
MouseLungMyeloid <- ScaleData(MouseLungMyeloid, verbose = T, features = all.genes)
MouseLungMyeloid <- RunPCA(MouseLungMyeloid, npcs = 30, verbose = FALSE)

#Calculate 90% Variance:
st_dev <- MouseLungMyeloid@reductions$pca@stdev
var <- st_dev^2
sum(var[1:18])/ sum(var)

MouseLungMyeloid <- FindNeighbors(object = MouseLungMyeloid, dims = 1:18)
MouseLungMyeloid <- FindClusters(object = MouseLungMyeloid, resolution = 2)
MouseLungMyeloid <- RunUMAP(MouseLungMyeloid, reduction = "pca", dims = 1:18, verbose = F)
DimPlot(MouseLungMyeloid, reduction = "umap", label = T)
DimPlot(MouseLungMyeloid, reduction = "umap", label = T, group.by = "manual_clusters")

DotPlot(MouseLungMyeloid, features=c("Ptprc", "Cd79a", "Cd19", "Ms4a1","Cd3e", "Cd4","Ccr7","Il2ra", "Foxp3","Ctla4", "Nkg7","Klrb1c","Ccr5", "Cd8a", "Itgam","Cd14","Ly6g", "Ly6c2", "S100a8","Acod1","Pdgfra","Il6", "Clec3b", "Il33", "Col1a1","Col1a2","Pdpn","Cd68","Adgre1","C1qa","Apoe", 
                                     "Mrc1","Fcgr3","H2-Eb1","Siglecf","Itgax","Mertk","Spp1","Ear2","Igha", "Sdc1","Itgae","Flt3", "Batf3","Xcr1","Tagln","Acta2","Pde5a","Clu","Pecam1", "Cdh5","Flt1","Kit", "Snca","Scgb1a1","Lamp3","Cxcl15","Aqp5", "Sftpd","Lyve1","Krt8", "Krt18", "Krt19","Sec14l3",
                                     "Foxj1","Msln", "Wt1", "Myl7", "Tnni3","Tnnt2","G0s2","Mcpt4","Mcpt8","Gata2", "Ms4a2","Ccna2", "Mki67","Hbb-bt"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

#Remove non-myeloid cells
MouseLungMyeloid <- subset(MouseLungMyeloid, idents = c("13","15","16"), invert = T) #remove contamination

MouseLungMyeloid <- FindVariableFeatures(MouseLungMyeloid, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(MouseLungMyeloid)
MouseLungMyeloid <- ScaleData(MouseLungMyeloid, verbose = T, features = all.genes)
MouseLungMyeloid <- RunPCA(MouseLungMyeloid, npcs = 30, verbose = FALSE)

#Calculate 90% Variance:
st_dev <- MouseLungMyeloid@reductions$pca@stdev
var <- st_dev^2
sum(var[1:18])/ sum(var)

MouseLungMyeloid <- FindNeighbors(object = MouseLungMyeloid, dims = 1:18)
MouseLungMyeloid <- FindClusters(object = MouseLungMyeloid, resolution = 2)
MouseLungMyeloid <- RunUMAP(MouseLungMyeloid, reduction = "pca", dims = 1:18, verbose = F)
DimPlot(MouseLungMyeloid, reduction = "umap", label = T)
DimPlot(MouseLungMyeloid, reduction = "umap", label = T, group.by = "manual_clusters")


#Idenfity Clusters
DotPlot(MouseLungMyeloid, features=c("Ptprc", "Cd79a", "Cd19", "Ms4a1","Cd3e", "Cd4","Ccr7","Il2ra", "Foxp3","Ctla4", "Nkg7","Klrb1c","Ccr5", "Cd8a", "Itgam","Cd14","Ly6g", "Ly6c2", "S100a8","Acod1","Pdgfra","Il6", "Clec3b", "Il33", "Col1a1","Col1a2","Pdpn","Cd68","Adgre1","C1qa","Apoe", 
                                     "Mrc1","Fcgr3","H2-Eb1","Siglecf","Itgax","Mertk","Spp1","Ear2","Igha", "Sdc1","Itgae","Flt3", "Batf3","Xcr1","Tagln","Acta2","Pde5a","Clu","Pecam1", "Cdh5","Flt1","Kit", "Snca","Scgb1a1","Lamp3","Cxcl15","Aqp5", "Sftpd","Lyve1","Krt8", "Krt18", "Krt19","Sec14l3",
                                     "Foxj1","Msln", "Wt1", "Myl7", "Tnni3","Tnnt2","G0s2","Mcpt4","Mcpt8","Gata2", "Ms4a2","Ccna2", "Mki67","Hbb-bt"), cols = "RdYlBu", dot.scale = 6) + RotatedAxis()

#Dendritic Cell
DotPlot(MouseLungMyeloid, features = c("Siglech", "Bst2", "Irf8","Batf3", "Clec9a", "Itgae", "Itgax", "Cadm1","Id2", "Irf4", "Klf4", "Cx3cr1", "Cd163", "Cd2"), cols = c("RdYlBu"), dot.scale = 8) + RotatedAxis()
#Active Dendritic Cell
DotPlot(MouseLungMyeloid, features = c("Ccl22", "Ccl17", "Ccr7"), cols = c("RdYlBu"), dot.scale = 8) + RotatedAxis()
#Neutrophil
DotPlot(MouseLungMyeloid, features = c("Itgam","Ptprc","Ly6g","Itga2b"), cols = c("RdYlBu"), dot.scale = 8) + RotatedAxis() 
#Macrophage
DotPlot(MouseLungMyeloid, features = c("Saa3","Ptprc", "Il1rl1","Arg1", "Retnla","Cd68","Fcgr3","Adgre1", "Itgam","Cd14","Mrc1", "S100a8", "Cd33","H2-Eb1",
                                     "Itgae","Clec9a","Batf3", "Ly6c2", "Ly6g","Ccna2", "Mki67", "Kit", "Hbb-bs"), cols = "RdYlBu", dot.scale = 8) + RotatedAxis()

#Label Manual Clusters:
Idents(MouseLungMyeloid) <- "seurat_clusters"
MouseLungMyeloid <- RenameIdents(MouseLungMyeloid,
                                  "0" = "Neutrophil", 
                                  "1" = "Neutrophil", 
                                  "2" = "Neutrophil", 
                                  "3" = "Alveolar Macrophage", 
                                  "4" = "Neutrophil", 
                                  "5" = "Macrophage", 
                                  "6" = "Neutrophil", 
                                  "7" = "Interstitial Macrophage", 
                                  "8" = "Macrophage",
                                  "9" = "Dendritic Cell", 
                                  "10" = "Neutrophil",
                                  "11" = "Interstitial Macrophage",
                                  "12" = "Neutrophil", 
                                  "13" = "Active Dendritic Cell")
DimPlot(MouseLungMyeloid, label = T)

new_order <- c("Macrophage",
               "Interstitial Macrophage",
               "Alveolar Macrophage",
               "Dendritic Cell",
               "Active Dendritic Cell",
               "Neutrophil")
MouseLungMyeloid@active.ident <- factor(MouseLungMyeloid@active.ident, levels = new_order)
MouseLungMyeloid[["manual_clusters"]] <- MouseLungMyeloid@active.ident

MouseLungMacrophage <- subset(MouseLungMyeloid, idents = c("Macrophage","Interstitial Macrophage","Alveolar Macrophage"))
MouseLungNeutrophil <- subset(MouseLungMyeloid, idents = "Neutrophil")



#### Figures ----------------------------------------------------------------------------------------------------####


### Figure 2A ###
DimPlot(MouseLung, cols = c("#DF4945","#D7A8FF","#528AFF","#52E7ED",
                               "#73D279","#C1D274","#FFB300","#FFE18F",
                               "#FFB2D8","#94435F","#7FA4CA","#6A5D99",
                               "#3A9255","#009193","#FF9360","#FFBEFF",
                               "#FF7E79"), pt.size = 0.5)+ NoLegend()


### Figure 2B ###
Idents(MouseLung) <- "Group"
DimPlot(MouseLung, cols = c("#EEC25A","gray"))

### Figure 2D ###
DimPlot(MouseEpithelial, pt.size = 2, cols = c("red1","coral","orange2","firebrick4")) + NoLegend()

### Figure 2E & Suplemental Figures 2A-B ###

#DE Total Epithelial
#ON vs OFF
MouseEpithelial_DE <- FindMarkers(MouseEpithelial, ident.1 = "Kras ON", ident.2 = "Kras OFF", group.by = "Group")

Idents(MouseEpithelial) <- "Group"

VlnPlot(MouseEpithelial, features = "Cxcl2", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseEpithelial, features = "Wnt4", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseEpithelial, features = "Ccl4", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseEpithelial, features = "Timp3", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseEpithelial, features = "Vegfa", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseEpithelial, features = "Fgf1", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseEpithelial, features = "Cxcl15", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseEpithelial, features = "Cxcl1", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseEpithelial, features = "Tgfa", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseEpithelial, features = "Tgfb1", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

### Figure 2F ###

genes <- c("Cxcl2","Ccl4","Wnt4","Mmp15","Fgf1","Igf1r","Vegfa","Hdgf","Hbegf","Cxcl15","Ccl2","Il13ra1",
                "Timp3","Timp2","Col4a1","Col1a2","Col4a2","Col4a3")

data <- FetchData(MouseEpithelial, vars = c(genes, "Group"))

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

### Figure 2G & Suplemental Figure 2C ###
#Import gene sets from msigdb. 
h_gene_sets = msigdbr(species = "mouse", category = "H")
msigdbr_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)

#Run DE with 0 thresholds:
Epithelial_DE <- FindMarkers(MouseEpithelial, ident.1 = "Kras ON", ident.2 = "Kras OFF", group.by = "Group", logfc.threshold = 0, min.pct = 0)

#Rank your DE by fold change:
ranks <- Epithelial_DE$avg_log2FC
names(ranks) <- rownames(Epithelial_DE)

#Run GSEA through fgsea:
fgseaRes <- fgsea(pathways = msigdbr_list, 
                  stats = ranks)

#Retrieve GSEA result:
result <- apply(fgseaRes,2,as.character)

#Plot and export a PDF of the GSEA visualization:
x <- print(plotEnrichment(msigdbr_list[["HALLMARK_KRAS_SIGNALING_DN"]], ranks, ticksSize = 0.5) + labs(title = "Lung Epithelial −  Kras ON vs Kras OFF", subtitle="HALLMARK_KRAS_SIGNALING_DN"))

x <- print(plotEnrichment(msigdbr_list[["HALLMARK_KRAS_SIGNALING_UP"]], ranks, ticksSize = 0.5) + labs(title = "Lung Epithelial −  Kras ON vs Kras OFF", subtitle="HALLMARK_KRAS_SIGNALING_UP"))


### Figure 3D ###
DimPlot(MouseLungFB, cols = c("#D7A8FF","#531B93","#CA79C8",
                              "#942193","#FFA3FF","#9085FF"), pt.size = 1.5)+NoLegend()


### Figure 3E & Supplemental Figure 3D-E ###
#DE
#ON vs OFF
MouseLungFB_DE <- FindMarkers(MouseLungFB_noPericytes, ident.1 = "Kras ON", ident.2 = "Kras OFF", group.by = "Group")

Idents(MouseLungFB_noPericytes) <- "Group"

VlnPlot(MouseLungFB_noPericytes, features = "Ccl4", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseLungFB_noPericytes, features = "Gas6", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseLungFB_noPericytes, features = "Cxcl2", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseLungFB_noPericytes, features = "Timp3", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseLungFB_noPericytes, features = "Tgfb3", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseLungFB_noPericytes, features = "Cxcl14", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseLungFB_noPericytes, features = "Col4a1", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseLungFB_noPericytes, features = "Fbln5", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseLungFB_noPericytes, features = "Cxcl1", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseLungFB_noPericytes, features = "Il6", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseLungFB_noPericytes, features = "Il33", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()


### Figure 3F ###
genes <- c("Ccl4","Gas6","Cxcl2","Timp3","Adamts2","Col1a1","Col1a2","Mmp3","Ccl8","Fgfr1","Col14a1",
           "Tgfb3","Fgfr3","Cxcl10","Col4a1","Fbln2","Pdgfrb","Mmp2","Serpine2","Ccl11","Fgf10",
           "Cxcl14","Tgfb2","Ccl5","Pdgfra","Ccl2","Fbn1","Fbln5","Lif","Vegfd","Tcf21","Tagln")

MouseLungFB_noPericytes$Group <- factor(MouseLungFB_noPericytes$Group, levels = c("Kras ON","Kras OFF"))

LungFB_data <- FetchData(MouseLungFB_noPericytes, vars = c(genes, "Group"))

df_avg <- data.frame()

for (id in levels(factor(LungFB_data$Group))) {
  data_subset <- LungFB_data %>% filter(Group == id)
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

colnames(df_avg_ratios) <- colnames(LungFB_data)[1:(ncol(LungFB_data)-1)]
rownames(df_avg_ratios) <- levels(factor(LungFB_data$Group))
pheatmap(as.matrix(df_avg_ratios), fontsize = 14, color = brewer.pal(9, "BuPu"), scale = 'none', cluster_rows = F, cluster_cols = T)


### Supplemental Figure 3B ###
DotPlot(MouseLungFB, features = c("Tagln","Acta2","Plin2","Fgf10","Col13a1","Tcf21",
                                  "Col14a1","Clec3b","Cyp1b1","Apod","Pdgfrb","Myl9"), cols = "PRGn")+RotatedAxis()

### Supplemental Figure 3C ###
#Cluster defining heatmap - top 10
Idents(MouseLungFB) <- "manual_clusters"

MouseLungFB.markers <- FindAllMarkers(MouseLungFB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- MouseLungFB.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(MouseLungFB, features = top10$gene, group.colors = c("#D7A8FF","#531B93","#CA79C8",
                                                               "#942193","#FFA3FF","#9085FF")) + NoLegend()

### Figure 4C ###
DimPlot(MouseLungMyeloid, cols = c("#52E7ED","#00A7A8","#005353","#7FA4CA","#0050A2","#528AFF"))+NoLegend()

### Figure 4D & Supplemental Figure 4E ###
#Macrophage DE On vs OFF:
MouseLungMac_DE <- FindMarkers(MouseLungMacrophage, ident.1 = "Kras ON", ident.2 = "Kras OFF",group.by = "Group")

Idents(MouseLungMacrophage) <- "Group"

VlnPlot(MouseLungMacrophage, features = "Apoe", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseLungMacrophage, features = "Mrc1", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseLungMacrophage, features = "Csf1r", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseLungMacrophage, features = "C1qa", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseLungMacrophage, features = "C1qb", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseLungMacrophage, features = "C1qc", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseLungMacrophage, features = "Il1a", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseLungMacrophage, features = "Tlr2", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseLungMacrophage, features = "Ifngr1", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseLungMacrophage, features = "Il1b", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

### Figure 4E ###
genes <- c("C1qa","C1qb","C1qc","Fos","Ccl4","Ccl8","Cx3cr1","Ccl7","Ccl2","Ccl22","Ccl6",
           "Ccr2","Csf1r","Ifnb1","Ifngr1","Ccl3","Il1b","Il1a","Tnfsf9","Mrc1","S100a9",
           "Tgfbi","Lifr","Apoe","Csf2rb","S100a1","Il6st","Cxcr4","Mmp19","Tlr2","Hif1a",
           "Il11ra1","S100a4","Csf1","Vegfa","Mmp12","Il6","Il7r")


data <- FetchData(MouseLungMacrophage, vars = c(genes, "Group"))

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

### Figure 4F ###
#Neutrophil DE On vs OFF:
MouseLungNeutrophil_DE <- FindMarkers(MouseLungNeutrophil, ident.1 = "Kras ON", ident.2 = "Kras OFF",group.by = "Group")

Idents(MouseLungNeutrophil) <- "Group"

VlnPlot(MouseLungNeutrophil, features = "Ccl3", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseLungNeutrophil, features = "Cd274", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

VlnPlot(MouseLungNeutrophil, features = "Icam1", split.by = "Group", cols = c("#EEC25A","gray")) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95) +NoLegend()

### Figure 4G ###
genes <- c("Mmp13","Il7r","Il11ra1","Il6st","Il1rn","Il1b","Il1a","Il6","Ccl8",
                  "Ccl2","Ccl4","Ccl6","Ccl3","Ccl7","Lif","Fn1","Tgfbi","Ifnb1","Ifngr1",
                  "Vegfa","Col14a1","Col6a3","Mmp19","Mmp12","Mmp14","Tnfsf9","Cxcr4",
                  "Ccr2","Cd274","Cd86","Icam1","Fcgrt","")

data <- FetchData(MouseLungNeutrophil, vars = c(genes, "Group"))

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
pheatmap(as.matrix(df_avg_ratios), fontsize = 14, color = brewer.pal(9, "PuBu"), scale = 'none', cluster_rows = F, cluster_cols = T)

### Figure 4H ###
VlnPlot(MouseLung, features = c("Egfr","Tgfbr1","Cxcr2"), ncol = 1, cols = c("#DF4945","#D7A8FF","#528AFF","#52E7ED",
                                                                             "#73D279","#C1D274","#FFB300","#FFE18F",
                                                                             "#FFB2D8","#94435F","#7FA4CA","#6A5D99",
                                                                             "#3A9255","#009193","#FF9360","#FFBEFF",
                                                                             "#FF7E79"))
