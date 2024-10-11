# load table which is specific barcode-associated cells with the number of times a given UMI is observed with each observed guide assignment
barcode_association <- read.csv("ko_barcodes.txt", sep = "\t")
barcodes <-  dplyr::as_tibble(barcode_association)
barcodes$count_column <-  barcodes$umi_count

# Prefiltering
barcodes <-  dplyr::as_tibble(barcode_association)
barcodes <-  barcodes %>% dplyr::filter(umi_count > 0)
barcodes <-  barcodes %>% filter(!grepl('unprocessed', barcode))
barcodes <-  barcodes %>%
  dplyr::group_by(cell) %>%
  dplyr::mutate(proportion=umi_count / sum(umi_count)) %>%
  dplyr::ungroup()

reads_threshold=2
proportion_threshold=0.05
barcodes.initial <- barcodes %>%
  dplyr::filter(proportion > proportion_threshold & umi_count > reads_threshold)

unique_gRNA_per_cell <- barcodes.initial %>%
  group_by(cell) %>%
  summarize(n_unique_gRNA = n_distinct(barcode))

single_gRNA_cells <- unique_gRNA_per_cell %>%
  filter(n_unique_gRNA == 1)  

barcodes.initial <- inner_join(barcodes.initial, single_gRNA_cells, by = "cell")

guide_metadata <-  read_excel("sgRNA_genename.xlsx")
guide_metadata <- guide_metadata %>% mutate(gRNA=toupper(gRNA))
barcodes.initial <-  dplyr::inner_join(barcodes.initial, guide_metadata, by=c("barcode"="gRNA"))
write.csv(barcodes.initial, "barcode_enrichment.csv")

#seurat data analysis
data <- Read10X(data.dir = "./CROP_scRNA/outs/filtered_feature_bc_matrix/")

# Create Seurat object
seurat_obj <- CreateSeuratObject(data, project = "CROP-seq", assay = "RNA")
# only keep the cells in seurat_obj
kept <- intersect(Cells(seurat_obj), barcodes.initial$cell)
seurat_obj <- subset(seurat_obj, cells = kept)

seurat_obj[["target"]] <- barcodes.initial$`gene name`
seurat_obj[["cell"]] <- barcodes.initial$cell
seurat_obj[["gRNA"]] <- barcodes.initial$barcode

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10)
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, reduction = "umap")
p1 <- FeaturePlot(seurat_obj, features = c("FOXO1","GABPA","HDAC3","HNF4A","PPARA","SCAP","YAP1"),min.cutoff="q1",max.cutoff="q90")
file_path <- file.path("./", "feature_plot.svg")
ggsave(p1, filename = file_path, width = 15, height = 10)

saveRDS(seurat_obj, "./CROP-seq_seurat.rds")

#Violin plot (visualization)
#PPAR_1 guide sequence are not identified
subSG <- as.character(c("FOXO1_1", "FOXO1_2", "Ctrl", "PPARA_2", "HNF4A_1", "HNF4A_2", "GABPA_1", "GABPA_2", "YAP1_1", "YAP1_2", "HDAC3_1", "HDAC3_2"))

gene_mapping <- list(
  "FOXO1_1" = "FOXO1", "FOXO1_2" = "FOXO1",
  "PPARA_2" = "PPARA", "HNF4A_1" = "HNF4A", "HNF4A_2" = "HNF4A",
  "GABPA_1" = "GABPA", "GABPA_2" = "GABPA", "YAP1_1" = "YAP1", "YAP1_2" = "YAP1",
  "HDAC3_1" = "HDAC3", "HDAC3_2" = "HDAC3"
)

plot_dir <- "./sgRNA_gene_expression/"
dir.create(plot_dir, showWarnings = FALSE)
metadata <- seurat_obj@meta.data
for (i in 1:length(subSG)) {
  gene <- subSG[i]
  if (gene != "Ctrl") {
    target_gene <- gene_mapping[[gene]]
    if (target_gene %in% rownames(seurat_obj)) {
      cellTarget <- metadata %>% filter(target == !!gene) %>% pull(cell)
      cellNTC <- metadata %>% filter(grepl("Ctrl", target)) %>% pull(cell)
      if (length(cellTarget) > 0 && length(cellNTC) > 0) {
        seurat_obj$type <- "NS"
        seurat_obj$type[Cells(seurat_obj) %in% cellTarget] <- "Target"
        seurat_obj$type[Cells(seurat_obj) %in% cellNTC] <- "NTC"
        plot <- VlnPlot(seurat_obj, features = target_gene, group.by = "type", y.max = 3) +
          geom_boxplot(width = 0.25, fill = "white") +
          stat_compare_means(label = "p.signif", comparisons = list(c("NS", "NTC"), c("NTC", "Target"), c("NS", "Target")))+
          ggtitle(paste("Violin Plot for", gene))
        
        plot_filename <- paste0(plot_dir, "Violin_Plot_", gene, ".svg")
        ggsave(plot_filename, plot = plot, width = 8, height = 6)
        
        print(plot)
      } else {
        print(paste("No cells found for target:", gene))
      }
    } else {
      print(paste("Feature", target_gene, "not found in Seurat object"))
    }
  }
}

