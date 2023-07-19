#' Annotate your scRNA-Seq data based on integrated helper T cell reference dataset
#'
#' @param query query Seurat object with PCA reduction already performed
#' @return Seurat object with two new metadata columns: predicted.CD4map.annotation, and predicted.CD4map.annotation.score


CD4map_annotate <- function(query) {
  path_to_kriukova_data <- "data/full_reference_return_model.rds"
  kriukova <- readRDS(path_to_kriukova_data) # my reference data. Will change in future

  anchors <- Seurat::FindTransferAnchors(
    reference = kriukova,
    query = query,
    reference.reduction = "pca",
    k.anchor = 5)

  query_scRNA2 <- Seurat::MapQuery(anchorset = anchors, reference = kriukova, query = query,
                           refdata = list(CD4map.annotation = "Subset"), reference.reduction = "pca", reduction.model = "umap")

  query_scRNA2
}

#' Plot your query scRNAseq data on top of the integrated helper T cell reference dataset
#' @param annotated_query Seurat object annotated with CD4map::CD4map_annotate() function
#' @param coloring_based_on Allowed values: "cluster", "score". If set to "cluster" (default), color encodes predicted cell annotations; if set to "score", color encodes the prediction score
#' @return ggplot object

CD4map_plot <- function(annotated_query, coloring_based_on = "cluster"){

  umap_tx_reference <- readRDS("data/umap_tx_reference.rds")

  color_code <- data.frame("predicted.CD4map.annotation" = c("EffMem_Th2", "Treg", "Naive_SOX4", "EffMem_Th1",
                                        "EffMem_Th17", "Naive_activated", "EffMem_Th22", "Naive",
                                        "CentMem", "Tfh", "Naive_IFN_response", "SelfStopping",
                                        "EffMem_Th2a", "EffMem_IFN_response", "Temra_cytotoxic_Th1", "Cycling"),
                           "colors" = c("#AA0DFE", "#3283FE", "#85660D", "#782AB6", "#565656", "#1C8356", "#16FF32",
                                                 "#F7E1A0", "red", "#1CBE4F", "#C4451C", "#DEA0FD", "#FE00FA", "#325A9B",
                                                 "#FEAF16", "#F8A19F"))

  umap_tx_query <- cbind(as.data.frame(annotated_query@reductions$ref.umap@cell.embeddings), annotated_query[[]])
  umap_tx_query$cellbarcode <- rownames(umap_tx_query)
  umap_tx_query <- dplyr::left_join(umap_tx_query, color_code, by = "predicted.CD4map.annotation")

  if (coloring_based_on == "cluster")
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = umap_tx_reference, ggplot2::aes(x=UMAP_1, y=UMAP_2), size = 0.2, stroke = 0.1, shape = 21, color = "grey90", alpha = 1) +
    ggplot2::theme_void() +
    ggplot2::geom_point(data = umap_tx_query, ggplot2::aes(x=refUMAP_1, y=refUMAP_2, color = colors), size = 0.2) +
    ggplot2::scale_colour_identity(guide = "legend", labels = color_code$predicted.CD4map.annotation, breaks = color_code$colors) +
    ggplot2::coord_fixed() +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=4, alpha = 1), title="CD4map annotation"))

  if (coloring_based_on == "score")
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = umap_tx_reference, ggplot2::aes(x=UMAP_1, y=UMAP_2), size = 0.2, stroke = 0.1, shape = 21, color = "grey90", alpha = 1) +
    ggplot2::theme_void() +
    ggplot2::geom_point(data = umap_tx_query, ggplot2::aes(x=refUMAP_1, y=refUMAP_2, color = predicted.CD4map.annotation.score), size = 0.2) +
    ggplot2::scale_color_continuous(type = "viridis") +
    ggplot2::coord_fixed() +
    ggplot2::guides(colour = ggplot2::guide_colourbar(title = "CD4map annotation score"))

  p
}

