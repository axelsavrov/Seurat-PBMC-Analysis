## Cargar paquetes de R
library("BiocFileCache") ## para descargar datos
library("dplyr") ## para filtrar datos
library("Seurat") ## paquete principal de este capitulo
library("patchwork") ## para combinar y crear gráficos, ya está.

# Usemos datos de pbmc3k tal y como lo hacen en
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html 
# pero con nuestro propio c??digo

# No necesitamos usar untar, ya que los archivos no est??n comprimidos.
# Solo especificamos la ruta donde est??n los archivos
fname <- "/Users/axelsalinas/filtered_gene_bc_matrices/hg19"

# Cargamos el dataset de PBMC
pbmc.data <- Read10X(data.dir = fname)
# cambiar los guiones bajos por los otros
rownames(pbmc.data) <- gsub("_", "-", rownames(pbmc.data))
# Inicializamos el objeto Seurat con los datos crudos (no normalizados)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

#vemos los primeros conteos de 3 genes
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

#Rrevisamos el porcentaje que corresponde a genes mitocondriales. O sea, no me aparecen como MT
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#visualizamos las m??tricas de control de calidad. Ah, tambi??n quiero saber bien a fondo para qu?? sirve STR
#EXPLICACI??N COMPLETA DEL VIOLIN PLOT (las 3 c??lulas HVG??s, ncount_RNA). OTRA DUDA, C??MO ACCEDO FORMALMENTE A LOS METADATOS
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#vamos a comparar pares de features
#para cada c??lula comparamos en n??mero de reads con el de mapeos mitocondriales
#A QUE SE REFIERE SOBRE ESCRIBIR OBJETOS
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Normalizacion de los datos -DUDA
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)

#Nos quedamos con los 2,000 genes m??s variables dentro del dataset
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#veamos los 10 genes m??s variables
top10 <- head(x= VariableFeatures(pbmc), 
              n= 10)

#graficamos los genes de acuerdo a su nivel de expresi??n
scater1 <- VariableFeaturePlot(pbmc)
#visualizaci??n
scater1
#destacamos los top 10 genes
scater2 <- LabelPoints( plot = scater1, 
                       points = top10, 
                       repel = TRUE ) 
#visualizar el scater2
scater2 
#obtener datos finales del objeto de seurat.                 DUDA CON ESTA PARTE
GetAssayData(pbmc)

# Normalizamos y escalamos los datos
pbmc <- ScaleData(pbmc)

# Realizamos PCA sobre los 2,000 genes más variables
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Visualizamos la varianza explicada por cada componente principal
ElbowPlot(pbmc)

# Calculamos los vecinos para el clustering (esto crea el grafo de vecinos)
pbmc <- FindNeighbors(pbmc, dims = 1:20)  # Usamos las primeras 20 componentes principales

# Aplicamos el algoritmo de clustering K-means
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Realizamos t-SNE para visualización
pbmc <- RunTSNE(pbmc, dims = 1:20)

# Visualizamos los clusters con etiquetas en el gráfico t-SNE
DimPlot(pbmc, reduction = "tsne", group.by = "seurat_clusters", label = TRUE)

# Guardamos la información de los clusters en los metadatos
head(pbmc@meta.data)

file.path(getwd(), "Seurat-PBMC-Analysis.R")


