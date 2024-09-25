#' @importFrom stats sd mean
#' @importFrom utils read.csv
#' @importFrom dimRed HLLE Isomap kPCA LEIM
#' @importFrom Rtsne Rtsne
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal labs
#' @importFrom DDoutlier LOF
#' @importFrom diffudist diffusionMap
#' @importFrom fitdistrplus fitdist
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom BiocParallel bplapply
#' @importFrom mclust Mclust
#' @importFrom cluster pam
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
#' @importFrom dplyr mutate
#' @importFrom ggrepel geom_text_repel
#' @importFrom aricode ARI
#' @importFrom phateR phate
#' @importFrom uwot umap

utils::globalVariables(c("Oihane", "Oihane.info.cellType", "Gutierrez", 
                         "Quake_Smartseq2_Lung", "mir", "gene", "cells", 
                         "cell.type", "anno"))

#' Get label mapping for a dataset
#'
#' This function returns the appropriate label mapping based on the provided dataset name.
#' It currently supports four datasets: "Oihane", "mir", "gene", and "Quake".
#'
#' @param dataset A string specifying the name of the dataset.
#'   Possible values are "Oihane", "mir", "gene", or "Quake".
#'
#' @return A named vector where names are the original labels and values are
#'   the corresponding abbreviations or simplified labels.
#'   Returns NULL if the provided dataset name is not recognized.
#'
#' @examples
#' get_mapping("Oihane")
#' get_mapping("Quake")
#'
#' @export
get_mapping <- function(dataset) {
  mapping_Oihane <- c(
    "Astrocytes" = "Astro",
    "Endothelial" = "Endo",
    "Ependymal" = "Epend",
    "Hybrid" = "Hyb",
    "Microglia" = "Micro",
    "Neurons" = "Neur",
    "Oligodendrocytes" = "Oligo"
  )

  mapping_Quake <- c(
    "B cell" = "BC",
    "ciliated columnar cell of tracheobronchial tree" = "CCCT",
    "classical monocyte" = "ClMono",
    "epithelial cell of lung" = "EpLung",
    "leukocyte" = "Leuk",
    "lung endothelial cell" = "LungEnd",
    "monocyte" = "Mono",
    "myeloid cell" = "Myel",
    "natural killer cell" = "NK",
    "stromal cell" = "Strom",
    "T cell" = "TC"
  )

  mapping_mir <- c(
    "AUTONOMIC_GANGLIA" = "AG",
    "BONE" = "Bo",
    "BREAST" = "Br",
    "CENTRAL_NERVOUS_SYSTEM" = "CNS",
    "ENDOMETRIUM" = "En",
    "FIBROBLAST" = "Fb",
    "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE" = "HL",
    "KIDNEY" = "Ki",
    "LARGE_INTESTINE" = "LIN",
    "LIVER" = "LIV",
    "LUNG" = "LUN",
    "OESOPHAGUS" = "Oe",
    "OVARY" = "Ov",
    "PANCREAS" = "Pa",
    "SKIN" = "Sk",
    "SOFT_TISSUE" = "STI",
    "STOMACH" = "STO",
    "THYROID" = "Th",
    "UPPER_AERODIGESTIVE_TRACT" = "UAT",
    "URINARY_TRACT" = "UT"
  )
    mapping_gene <- c(
    "AUTONOMIC_GANGLIA" = "AG",
    "BONE" = "Bo",
    "BREAST" = "Br",
    "CENTRAL_NERVOUS_SYSTEM" = "CNS",
    "ENDOMETRIUM" = "En",
    "FIBROBLAST" = "Fb",
    "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE" = "HL",
    "KIDNEY" = "Ki",
    "LARGE_INTESTINE" = "LIN",
    "LIVER" = "LIV",
    "LUNG" = "LUN",
    "OESOPHAGUS" = "Oe",
    "OVARY" = "Ov",
    "PANCREAS" = "Pa",
    "SKIN" = "Sk",
    "SOFT_TISSUE" = "STI",
    "STOMACH" = "STO",
    "THYROID" = "Th",
    "UPPER_AERODIGESTIVE_TRACT" = "UAT",
    "URINARY_TRACT" = "UT"
  )
  switch(dataset,
         "Oihane" = mapping_Oihane,
         "mir" = mapping_mir,
         "gene" = mapping_gene,
         "Quake" = mapping_Quake,
         NULL)
}

#' Load and Prepare Dataset for Analysis
#'
#' @param dataset A string specifying the dataset to load. Supported datasets are:
#'   "Oihane", "Gutierrez", "Spleen", "Quake", "pbmc", "metabolism", "kidney", "mir", "gene".
#' @param base_dir A string specifying the base directory for data files. Default is ".".
#'
#' @return A list containing:
#'   \item{dat}{The main data matrix or data frame}
#'   \item{cell.type}{A factor of cell types}
#'   \item{info}{A factor of additional information (often identical to cell.type)}
#'
#' @details
#' This function loads and prepares various datasets for analysis. Each dataset
#' is processed differently based on its structure and content. The function sets
#' the working directory to a specific path for each dataset before loading.
#'
#' @note
#' - The function assumes a specific directory structure and file naming convention.
#' - Some datasets require additional files to be present in the specified directories.
#' - For "mir" and "gene" datasets, extensive preprocessing is performed.
#'
#' @examples
#' \dontrun{
#' # Load the Oihane dataset
#' oihane_data <- dataloader("Oihane")
#'
#' # Load the PBMC dataset
#' pbmc_data <- dataloader("pbmc")
#'
#' # Load the miRNA dataset
#' mir_data <- dataloader("mir")
#' }
#'
#' @export
dataloader <- function(dataset, base_dir = ".") {
    path = paste0("../data/",dataset,"/")
    print(path)
  setwd(path)

  load_common <- function() {
    load(paste0(path,"all embeded dim 3.bin"))
    load(file.path("Data", paste(dataset, "metaspec_out.Rdata")))  # ensemble.out
    load(file.path("Data", paste(dataset, "adm_out.Rdata")))  # mev.out
  }

  result <- switch(dataset,
    "Oihane" = {
      load("Oihane.Rdata")
      list(
        dat = Oihane,
        cell.type = as.factor(Oihane.info.cellType),
        info = as.factor(Oihane.info.cellType)
      )
    },
    "Gutierrez" = {
      load("Gutierrez.Rdata")
      list(
        dat = Gutierrez$data,
        cell.type = as.factor(Gutierrez$data.cellType),
        info = as.factor(Gutierrez$data.cellType)
      )
    },
     "Spleen" = {
       dat <- utils::read.csv("./Spleen_pro_data.csv", header = TRUE, row.names = 1)
       lab <- utils::read.csv("./Spleen_pro_label.csv", header = TRUE, row.names = 1)
       lab <- lab$SpatialGlue
      list(
         dat = dat,
         cell.type = lab,
         info = lab
      )
    },
    "Quake" = {
      load("Quake_Smartseq2_Lung.Rdata")
      list(
        dat = Quake_Smartseq2_Lung$data,
        cell.type = as.factor(Quake_Smartseq2_Lung$data.cellType),
        info = as.factor(Quake_Smartseq2_Lung$data.cellType)
      )
    },
    "Brain5k" = {
      dat <- utils::read.csv("Brain5k_data.csv", header = TRUE, row.names = 1)
      lab <- utils::read.csv("Brain5k_metadata.csv", header = TRUE, row.names = 1)[,1]
      list(
        dat = dat,
        cell.type = as.factor(lab),
        info = as.factor(lab)
      )
    },
     "metabolism" = {
      dat <- utils::read.csv("./metabolism_data.csv", header = TRUE, row.names = 1)
      lab <- utils::read.csv("./metabolism_label.csv", header = TRUE, row.names = 1)
      lab<-lab$subtype
      list(
        dat = dat,
        cell.type = as.factor(lab),
        info = as.factor(lab)
      )
    },
    "kidney" = {
      dat <- utils::read.csv("kidney_data.csv", header = TRUE, row.names = 1)
      lab <- utils::read.csv("kidney_metadata.csv", header = TRUE, row.names = 1)$cell_type
      list(
        dat = dat,
        cell.type = as.factor(lab),
        info = as.factor(lab)
      )
    },
     "mir" = {
      load("CCLE_miRNA_20181103.bin")
      load("CCLE_RNAseq_genes_rpkm_20180929.bin")

      mir <- log10(1+mir)
      gene <- log10(1+gene)
      gene <- as.matrix(gene)
      mir <- as.matrix(mir)

      mir <- mir[,which(colnames(mir) %in% colnames(gene))]
      gene <- gene[,which(colnames(gene) %in% colnames(mir))]

      mir <- mir[,order(colnames(mir))]
      gene <- gene[,order(colnames(gene))]

      sum(colnames(mir)==colnames(gene))
      dim(mir)

      mir <- t(mir)
      gene <- t(gene)

      cv.mir <- apply(mir,2,stats::sd)/apply(mir,2,mean)
      mir <- mir[,cv.mir>=0.1]

      ze <- apply(gene==0,2,sum)
      gene <- gene[,which(ze<=0.25*nrow(gene))]
      cv.gene <- apply(gene,2,stats::sd)/apply(gene,2,mean)
      gene <- gene[,cv.gene>=0.5]

      cells <- rownames(gene)
      cell.type <- cells
      for(i in 1:length(cells))
      {
          this <- strsplit(cells[i], "_")[[1]][-1]
          this <- paste(this, collapse="_")
          cell.type[i] <- this
      }

      ttt <- table(cell.type)
      sel <- which(cell.type %in% names(ttt)[ttt<10])
      cell.type[sel] <- NA

      anno <- utils::read.table("Cell_lines_annotations_20181226.txt",header=T,sep="\t")
      anno <- anno[which(anno[,1] %in% cells),]
      anno <- anno[order(anno[,1]),]

      s.anno <- which(cells %in% anno[,1] & !is.na(cell.type))
      cells <- cells[s.anno]
      cell.type <- cell.type[s.anno]
      mir <- mir[s.anno,]
      gene <- gene[s.anno,]

      mir2 <- mir
      gene2 <- gene

      for(i in 1:ncol(mir2)) mir2[,i] <- (mir2[,i]-mean(mir2[,i]))/stats::sd(mir2[,i])
      for(i in 1:ncol(gene2)) gene2[,i] <- (gene2[,i]-mean(gene2[,i]))/stats::sd(gene2[,i])

      all.colors <- c("grey50","green","blue","cyan", "yellow","orange","red","black", "wheat","purple","darkblue","dodgerblue4","darkred","darkorange","darkcyan","magenta","firebrick","khaki4","cornsilk3","darkgoldenrod4")

      dat <- mir
      info <- as.factor(cell.type)

      list(
        dat = dat,
        cell.type = cell.type,
        info = info
      )
    },
     "gene" = {
      load("CCLE_miRNA_20181103.bin")
      load("CCLE_RNAseq_genes_rpkm_20180929.bin")

      mir <- log10(1+mir)
      gene <- log10(1+gene)
      gene <- as.matrix(gene)
      mir <- as.matrix(mir)

      mir <- mir[,which(colnames(mir) %in% colnames(gene))]
      gene <- gene[,which(colnames(gene) %in% colnames(mir))]

      mir <- mir[,order(colnames(mir))]
      gene <- gene[,order(colnames(gene))]

      sum(colnames(mir)==colnames(gene))
      dim(mir)

      mir <- t(mir)
      gene <- t(gene)

      cv.mir <- apply(mir,2,stats::sd)/apply(mir,2,mean)
      mir <- mir[,cv.mir>=0.1]

      ze <- apply(gene==0,2,sum)
      gene <- gene[,which(ze<=0.25*nrow(gene))]
      cv.gene <- apply(gene,2,stats::sd)/apply(gene,2,mean)
      gene <- gene[,cv.gene>=0.5]

      cells <- rownames(gene)
      cell.type <- cells
      for(i in 1:length(cells))
      {
          this <- strsplit(cells[i], "_")[[1]][-1]
          this <- paste(this, collapse="_")
          cell.type[i] <- this
      }

      ttt <- table(cell.type)
      sel <- which(cell.type %in% names(ttt)[ttt<10])
      cell.type[sel] <- NA

      anno <- utils::read.table("Cell_lines_annotations_20181226.txt",header=T,sep="\t")
      anno <- anno[which(anno[,1] %in% cells),]
      anno <- anno[order(anno[,1]),]

      s.anno <- which(cells %in% anno[,1] & !is.na(cell.type))
      cells <- cells[s.anno]
      cell.type <- cell.type[s.anno]
      mir <- mir[s.anno,]
      gene <- gene[s.anno,]

      mir2 <- mir
      gene2 <- gene

      for(i in 1:ncol(mir2)) mir2[,i] <- (mir2[,i]-mean(mir2[,i]))/stats::sd(mir2[,i])
      for(i in 1:ncol(gene2)) gene2[,i] <- (gene2[,i]-mean(gene2[,i]))/stats::sd(gene2[,i])

      all.colors <- c("grey50","green","blue","cyan", "yellow","orange","red","black", "wheat","purple","darkblue","dodgerblue4","darkred","darkorange","darkcyan","magenta","firebrick","khaki4","cornsilk3","darkgoldenrod4")

      dat <- gene
      info <- as.factor(cell.type)

      list(
        dat = dat,
        cell.type = cell.type,
        info = info
      )
    },
    stop("Unsupported dataset")
  )
  return(result)
}