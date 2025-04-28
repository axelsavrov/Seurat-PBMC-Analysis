# Análisis de Datos de Células PBMC con Seurat

Este proyecto realiza un análisis de datos de células PBMC utilizando el paquete **Seurat** en R. El análisis incluye la carga de datos, la normalización, el escalado, el análisis de variabilidad, la reducción de dimensiones (PCA), la identificación de clusters, y la visualización de los resultados con t-SNE.

## Requisitos

Para ejecutar este análisis, necesitas tener instalado R y los siguientes paquetes:

- **BiocFileCache**
- **dplyr**
- **Seurat**
- **patchwork**

Puedes instalar los paquetes necesarios con los siguientes comandos en R:

```R
install.packages("dplyr")
install.packages("patchwork")
BiocManager::install("BiocFileCache")
install.packages("Seurat")


