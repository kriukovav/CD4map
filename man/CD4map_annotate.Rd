% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{CD4map_annotate}
\alias{CD4map_annotate}
\title{Annotate your scRNA-Seq data based on integrated helper T cell reference dataset}
\usage{
CD4map_annotate(query)
}
\arguments{
\item{query}{query Seurat object with PCA reduction already performed}
}
\value{
Seurat object with two new metadata columns: predicted.CD4map.annotation, and predicted.CD4map.annotation.score
}
\description{
Annotate your scRNA-Seq data based on integrated helper T cell reference dataset
}
