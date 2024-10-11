# libyary R package
library(Seurat)
library(dplyr)
library(readxl)

data <- Read10X(data.dir = "~/CROP-seq/new_data/Count_output/CROP_scRNA/outs/filtered_feature_bc_matrix/")
# Create Seurat object
seurat_obj <- CreateSeuratObject(data, project = "scCROP-seq", assay = "scRNA")
# Functions for preprocessing of KO barcode data
# Can use with UMIs (umi = TRUE) or reads (UMI=FALSE)
barcode_association <- read.csv("~/CROP-seq/new_data/Count_output/ko_barcodes.txt", sep = "\t")
barcodes <- dplyr::as_tibble(barcode_association)
barcodes$count_column <- barcodes$umi_count

# Prefiltering
barcodes <- dplyr::as_tibble(barcode_association)
barcodes <- barcodes %>% dplyr::filter(umi_count > 0)
barcodes <- barcodes %>% filter(!grepl('unprocessed', barcode))

barcodes = barcodes %>%
  dplyr::group_by(cell) %>%
  dplyr::mutate(proportion=umi_count / sum(umi_count)) %>%
  dplyr::ungroup()


# Make sure case matches
barcodes = barcodes %>% mutate(barcode=toupper(barcode))

proportion_threshold=0.05
reads_threshold=2

barcodes.initial <- barcodes %>%
  dplyr::filter(proportion > proportion_threshold & umi_count > reads_threshold)

# guide_metadata <-  read_excel("sgRNA_genename.xlsx")
guide_metadata <- guide_metadata %>% mutate(gRNA=toupper(gRNA))
barcodes.initial <- dplyr::inner_join(barcodes.initial, guide_metadata, by=c("barcode"="gRNA"))
#check whether there is any cell with more than one guideRNA
unique_gRNA_per_cell <- barcodes.initial %>%
  group_by(cell) %>%
  summarize(n_unique_gRNA = n_distinct(`gene name`))

single_gRNA_cells <- unique_gRNA_per_cell %>%
  filter(n_unique_gRNA == 1)

kept <- intersect(Cells(seurat_obj), single_gRNA_cells$cell)
seurat_obj <- subset(seurat_obj, cells = kept)

# Rearrange data frame so we can add this information to the Seurat object
kept_cells <- barcodes.initial[match(Cells(seurat_obj), single_gRNA_cells$cell), ]

seurat_obj[["target"]] <- kept_cells$`gene name`
seurat_obj[["cell"]] <- kept_cells$cell
seurat_obj[["gRNA"]] <- kept_cells$barcode

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL")
seurat_obj <- subset(seurat_obj, subset = nFeature_scRNA > 200 & nFeature_scRNA < 10000 & percent.mt < 10)
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
# Normalise using SCtransform, with regression for percent.mt and percent.rb,
# # in addition to batch normalisation
# seurat_obj <- SCTransform(seurat_obj, assay = "RNA", new.assay.name = "SCT",
#                           vars.to.regress = c("percent.mt", "percent.rb"),
#                           conserve.memory = TRUE)
# 
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, reduction = "umap")
p1 <- FeaturePlot(seurat_obj, features = c("FOXO1","GABPA","HDAC3","HNF4A","PPARA","SCAP","YAP1"),min.cutoff="q1",max.cutoff="q90")
ggsave
saveRDS(seurat_obj, "CROP-seq_seurat.rds")
