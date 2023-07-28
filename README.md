# CD4map
This package provides wrapper functions to annotate helper T cells in scRNA-Seq data based on reference dataset to be published by Chudakov lab
## Installation

```R
install.packages('devtools')
devtools::install_github("kriukovav/CD4map")
```
### Important: please download the reference
Please download the reference dataset full_reference_return_model.rds, available from [here](https://figshare.com/projects/T_helper_subsets_Kriukova_et_al_/173466), and save it inside the "data" subdirectory of the CD4map package localy on your device. To find the location of the CD4map package, simply run

```R
.libPaths()
```

## Usage

### Annotate your data
With the CD4map_annotate() wrapper you may easily apply [Seurat](https://satijalab.org/seurat/index.html) reference mapping algorithm to map your data to our integrated reference dataset. Just provide your Seurat object (with PCA already performed) as a single argument

```R
name_of_your_annotated_seurat_object <- CD4map_annotate(query = name_of_your_seurat_object)
```
This function will return you your Seurat object with two additional metadata columns: predicted.CD4map.annotation, and predicted.CD4map.annotation.score

### Plot your data
With the CD4map_plot() wrapper you may plot your scRNA-Seq data on top of our Th reference. 

```R
CD4map_plot(annotated_query = name_of_your_annotated_seurat_object)

CD4map_plot(annotated_query = name_of_your_annotated_seurat_object, coloring_based_on == "score") # will plot the annotation score

```
This function will return you a ggplot object, which you may further adjust according to your needs.

