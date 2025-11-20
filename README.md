#     Spatial-Single-Cell-Analysis-Using-Salmon-Alevin-fry
## 1 Generate USA quantification with alevin-fry

###  Environment Setup

### Create and activate conda environment

```
conda create -n Spatial_single_cell python=3.10 -y
conda activate Spatial_single_cell
```
### Install required tools
```
conda install -c bioconda -c conda-forge \
  simpleaf \
  alevin-fry \
  salmon \
  gffread \
  -y
```
### Check versions
```
simpleaf --version
alevin-fry --version
salmon --version
gffread --version
```
###  Download Reference Genome & Annotation

```
wget -O gencode.vM38.annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.annotation.gtf.gz

wget -O GRCm39.primary_assembly.genome.fa.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/GRCm39.primary_assembly.genome.fa.gz

gunzip -f gencode.vM38.annotation.gtf.gz
gunzip -f GRCm39.primary_assembly.genome.fa.gz
```
###  Generate Transcriptome FASTA

```
gffread -w gencode.vM38.transcriptome.fa -g GRCm39.primary_assembly.genome.fa gencode.vM38.annotation.gtf

```


### Optional: inspect
```
head gencode.vM38.transcriptome.fa
head GRCm39.primary_assembly.genome.fa
head gencode.vM38.annotation.gtf
```

### Create Transcript-to-Gene Mapping
```
awk '{
  if ($3 == "transcript") {
    match($0, /gene_id "([^"]+)"/, g);
    match($0, /transcript_id "([^"]+)"/, t);
    if (t[1] && g[1]) print t[1] "\t" g[1];
  }
}' gencode.vM38.annotation.gtf > t2g.tsv
```

### Inspect
```
head -5 t2g.tsv
```


### Build Salmon Transcriptome Index
```
salmon index -t ref/gencode.vM38.transcriptome.fa -i grch38_idx -p 4
```


### Quantify Reads with Salmon Alevin-fry
```
salmon alevin -l ISR -i grch38_idx -p 4 --chromiumV3 \
  -1 fastq/V1_Mouse_Brain_Sagittal_Anterior_Section_1_S7_L001_R1_001.fastq.gz \
     fastq/V1_Mouse_Brain_Sagittal_Anterior_Section_1_S7_L002_R1_001.fastq.gz \
  -2 fastq/V1_Mouse_Brain_Sagittal_Anterior_Section_1_S7_L001_R2_001.fastq.gz \
     fastq/V1_Mouse_Brain_Sagittal_Anterior_Section_1_S7_L002_R2_001.fastq.gz \
  -o alevin_out \
  --tgMap ref/t2g.tsv \
  --sketch
```


### Generate permit list
```
alevin-fry generate-permit-list -k -d fw -i alevin_out -o alevin_out/permitlist_knee/
```


### Collate
```
alevin-fry collate -i alevin_out/permitlist_knee/ -r alevin_out -t 4
```


### Quantify with cr-like resolution
```
alevin-fry quant -r cr-like --use-mtx \
  -m ref/t2g.tsv \
  -i alevin_out/permitlist_knee/ \
  -o alevin_out/quant_cr-like_usa \
  -t 4
```
## 2 Load Alevin-fry into R as a Seurat object

### load an R package
```
library(fishpond)
library(Seurat)
library(ggplot2)
library(patchwork)
```


### Load alevin-fry outputR
```
brain_sce <- loadFry(
  "alevin_out/quant_cr-like_usa/",
  outputFormat = "scRNA",
  quiet = FALSE
)

brain <- as.Seurat(brain_sce, data = NULL)
```




### Rename first assay to "Spatial"
```
names(brain@assays)[1] <- "Spatial"
DefaultAssay(brain) <- "Spatial"
```




### Add basic QC
```
brain[["nCount_Spatial"]] <- colSums(brain@assays$Spatial@counts)
brain[["nFeature_Spatial"]] <- colSums(brain@assays$Spatial@counts > 0)
```




### SCTransform (optional)
```
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
```




### Read spatial image and coordinates
```
spatial.img <- Read10X_Image("spatial/")

coords <- GetTissueCoordinates(spatial.img)
coords_raw_barcodes <- rownames(coords)
```




### Remove trailing "-1"
```
coords_barcodes <- gsub("-1$", "", coords_raw_barcodes)

brain_barcodes <- Cells(brain)
```




### Match alevin barcodes with Visium barcodes
```
common <- intersect(brain_barcodes, coords_barcodes)
cat("Matched barcodes:", length(common), "of", length(brain_barcodes), "\n")

brain <- subset(brain, cells = common)

coords_subset <- coords[match(common, coords_barcodes), ]
```




### Fix barcode suffix
```
new_brain_barcodes <- paste0(common, "-1")
brain <- RenameCells(brain, new.names = new_brain_barcodes)

spatial.img <- spatial.img[new_brain_barcodes]
```



### Attach image
```
brain[["slice1"]] <- spatial.img
```




### Plot
```
SpatialDimPlot(brain)
SpatialFeaturePlot(brain, features = "nCount_Spatial", pt.size.factor = 1.3)
```


### UMAP + Spatial combined plot
```
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + 
  NoLegend()

plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + 
  theme(legend.position = "right")

wrap_plots(plot1, plot2)

ggsave("qc_vln_spatial.png", wrap_plots(plot1, plot2),
       width = 10, height = 6, dpi = 300)
```


### Hpca / Ttr SpatialFeaturePlot
```
gtf <- rtracklayer::import("ref/gencode.vM38.annotation.gtf")

hpca_id <- unique(gtf$gene_id[gtf$gene_name == "Hpca"])
ttr_id  <- unique(gtf$gene_id[gtf$gene_name == "Ttr"])

p_gene <- SpatialFeaturePlot(brain, features = c(hpca_id, ttr_id))

ggsave("Hpca_Ttr_spatial.png", p_gene,
       width = 8, height = 4, dpi = 300)
```


### PCA, Neighbors, Clustering, UMAP
```
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)

brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, resolution = 0.45, verbose = FALSE)

brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

```


### UMAP and Spatial plots
```
p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)

p_umap <- p1 + p2


ggsave("UMAP_and_spatial.png", p_umap,
       width = 10, height = 5, dpi = 300)
```


### Highlight selected clusters on Spatial plot
```
p_highlight <- SpatialDimPlot(
  brain,
  cells.highlight = CellsByIdentities(
    object = brain,
    idents = c(2, 3, 5, 0, 6, 7)
  ),
  facet.highlight = TRUE,
  ncol = 3
)

ggsave("highlight_clusters.png", p_highlight,
       width = 10, height = 6, dpi = 300)
```
