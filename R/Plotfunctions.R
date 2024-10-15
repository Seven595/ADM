#' @importFrom grDevices rainbow
#' @importFrom graphics pairs
#' @importFrom methods new
#' @importFrom stats as.dist cmdscale cor density dist embed quantile rgamma sd
#' @importFrom utils read.table
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import ggrepel
#' @importFrom uwot umap
#' @importFrom mclust adjustedRandIndex
#' @importFrom cluster silhouette
#' @importFrom aricode NMI
#' @importFrom Rtsne Rtsne
#' @importFrom BiocParallel bplapply

#' Get label mapping for a dataset
#'
#' This function returns the appropriate label mapping based on the provided dataset name.
#' It currently supports four datasets: "Oihane", "mir", "gene", and "Quake".
#'
#' @param dataset A string specifying the name of the dataset.
#' Possible values are "Oihane", "mir", "gene", or "Quake".
#'
#' @return A named vector where names are the original labels and values are
#' the corresponding abbreviations or simplified labels.
#' Returns NULL if the provided dataset name is not recognized.
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

GLOBAL_PARAMS <- list(
  n_components = 3,
  k_values = c(1, 2, 5, 10, 20),
  seed = 2024,
  figure_path = "./figures/",
  font_sizes = list(
    title = 36,
    axis = 26,
    axis_cci = 26,
    legend_title = 24,
    legend_text = 24,
    caption = 12
  )
)
color_list =c("#FB6A4A","#54278F","#006635","#3182BD","#DE2D26","#72A34F","#5D7AD3", "#756BB1","#FCAE91","#fe87ac","#AFABAB","#67A9CF","#CBC9E2","#4d982e","#E6873E","#545454","#aa3474","#ee8c7d","#2e5fa1","#FDD0A3","#C22F2F","#036f73")
names_list = c("PCA","MDS","iMDS","Sammon", "HLLE","Isomap","kPCA1","kPCA2","LEIM" , "UMAP1" , "UMAP2" , "tSNE1", "tSNE2", "PHATE1", "PHATE2","KEF")

#' Visualize data using ggplot2
#'
#' This function creates a scatter plot of the input data, optionally coloring points by group.
#'
#' @param data A matrix or data frame with at least two columns for x and y coordinates.
#' @param method_name A string specifying the name of the visualization method (used for plot title).
#' @param color_list An optional named vector of colors for each group in 'info'.
#' @param info An optional vector of group labels for each data point.
#'
#' @return A ggplot object representing the scatter plot.
#'
#' @keywords internal
visualization_func <- function(data, method_name, color_list = NULL, info = NULL) {
  data <- as.data.frame(data)
  colnames(data)[1:2] <- c("x", "y")
  data$info <- info
  plot <- ggplot(data, aes(x = .data$x, y = .data$y)) +
    labs(title = method_name, x = "", y = "") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )

  if (!is.null(info) && !is.null(color_list)) {
    data$info <- info
    plot <- plot +
      geom_point(aes(color = .data$info), alpha = 0.75) +
      scale_color_manual(values = color_list) +
      guides(color = guide_legend(override.aes = list(size = 6, shape = 15), ncol = 2)) +
      theme(legend.key.height = unit(1, "cm"))
  } else {
    plot <- plot + geom_point(alpha = 0.75)
  }

  return(plot)
}

#' Create a labeled scatter plot from UMAP visualization
#'
#' This function generates a scatter plot from UMAP coordinates with labeled cluster centers.
#'
#' @param umap_viz A matrix or data frame with UMAP coordinates (at least two columns).
#' @param method A string specifying the visualization method name (used for plot title).
#' @param info A vector of cluster labels for each data point.
#' @param dataset A string specifying the dataset name (not used in the current implementation).
#' @param color_list A named vector of colors for each cluster.
#'
#' @return A ggplot object representing the labeled scatter plot.
#'
#' @details
#' The function creates a scatter plot of UMAP coordinates, colors points by cluster,
#' calculates and labels cluster centers, and applies a custom theme. It uses ggrepel
#' for non-overlapping text labels.
#'
#' @examples
#' umap_coords <- matrix(rnorm(200), ncol = 2)
#' cluster_info <- rep(c("A", "B"), each = 50)
#' color_list <- c("A" = "red", "B" = "blue")
#' p <- visualization_with_label(umap_coords, "UMAP", cluster_info, "Example", color_list)
#' if (interactive()) {
#'   print(p)
#' }
#'
#' @export
visualization_with_label <- function(umap_viz, method, info, dataset, color_list) {
  set.seed(2024)
  umap_df <- as.data.frame(umap_viz)
  colnames(umap_df)[1:2] <- c("V1", "V2")
  umap_df$cluster <- info

  cluster_centers <- umap_df %>%
    group_by(.data$cluster) %>%
    summarise(center_x = mean(.data$V1),
              center_y = mean(.data$V2)) %>%
    mutate(
      label_x = .data$center_x + (.data$center_x - mean(.data$center_x)) * 0.2,
      label_y = .data$center_y + (.data$center_y - mean(.data$center_y)) * 0.2
    )

  p <- ggplot(umap_df, aes(x = .data$V1, y = .data$V2, color = .data$cluster)) +
    geom_point(alpha = 0.75) +
    scale_color_manual(values = color_list) +
    geom_segment(data = cluster_centers,
                 aes(x = .data$center_x, y = .data$center_y, xend = .data$label_x, yend = .data$label_y),
                 color = "grey50", linetype = "dashed") +
    geom_text_repel(data = cluster_centers,
                    aes(x = .data$label_x, y = .data$label_y, label = .data$cluster),
                    color = "black",
                    fontface = "bold",
                    box.padding = 0.5,
                    point.padding = 0.5,
                    segment.color = NA, size = 8) +
    labs(title = paste(method),
         x = "",
         y = "") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.key.height = unit(0.8, "cm"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "bottom"
    ) +
    guides(color = guide_legend(override.aes = list(size = 6, shape = 15), ncol = 2))

  return(p)
}

#' Create a scatter plot with a custom legend
#'
#' This function generates a scatter plot from UMAP coordinates with a customized legend.
#'
#' @param umap_viz A matrix or data frame with UMAP coordinates (at least two columns).
#' @param method_name A string specifying the visualization method name (used for plot title).
#' @param info A vector of original cluster labels for each data point.
#' @param dataset A string specifying the dataset name (not used in the current implementation).
#' @param color_list A named vector of colors for each cluster.
#' @param label_mapping An optional named vector for mapping original labels to simplified ones.
#'
#' @return A ggplot object representing the scatter plot with a customized legend.
#'
#' @details
#' The function creates a scatter plot of UMAP coordinates, colors points by cluster,
#' and applies a custom theme with a specially formatted legend. It uses the 'simplify_labels'
#' function to potentially simplify cluster labels.
#'
#' @examples
#' umap_coords <- matrix(rnorm(200), ncol = 2)
#' cluster_info <- rep(c("Type A", "Type B"), each = 50)
#' color_list <- c("Type A" = "red", "Type B" = "blue")
#' label_map <- c("Type A" = "A", "Type B" = "B")
#' p <- plot_legend(umap_coords, "UMAP", cluster_info, "Example", color_list, label_map)
#' if (interactive()) {
#'   print(p)
#' }
#'
#' @export
plot_legend <- function(umap_viz, method_name, info, dataset, color_list, label_mapping = NULL) {
  set.seed(2024)
  umap_df <- as.data.frame(umap_viz)
  simplified_info <- simplify_labels(info, label_mapping)
  umap_df$cluster <- simplified_info

  p <- ggplot(umap_df, aes(x = .data[[1]], y = .data[[2]], color = .data$cluster)) +
    geom_point(alpha = 0.75) +
    scale_color_manual(values = color_list) +
    labs(title = method_name,
         x = "",
         y = "") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right",
      legend.key.height = unit(1.3, "cm"),
      legend.spacing.x = unit(1, "cm"),
      legend.text = element_text(margin = margin(l = -0.5, unit = "cm")),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 8, shape = 15), ncol = 2))

  return(p)
}

#' Plot Comparison of CCI Results
#'
#' This function processes and visualizes the results of Cell-Cell Interaction (CCI) analyses
#' for different methods and data transformations.
#'
#' @param cci_list A list containing two matrices: one for visualization method and one for diffusion method.
#' @param dataset A string specifying the name of the dataset (used for plot title).
#'
#' @return A ggplot object representing a boxplot of CCI results.
#'
#' @details
#' The function processes CCI matrices for visualization and diffusion methods,
#' combines the data, and creates a boxplot comparing the results across different
#' data types (raw, UMAP, PCA) and methods (metaspec, ADM).
#'
#' @examples
#' viz_matrix <- matrix(rnorm(30), nrow = 10,
#'                      dimnames = list(NULL, c("viz_raw", "viz_umap", "viz_pca")))
#' diffu_matrix <- matrix(rnorm(30), nrow = 10,
#'                        dimnames = list(NULL, c("diffu_raw", "diffu_umap", "diffu_pca")))
#' cci_list <- list(viz_matrix, diffu_matrix)
#' p <- plot_cci_results(cci_list, "Example Dataset")
#' if (interactive()) {
#'   print(p)
#' }
#'
#' @export
plot_cci_results <- function(cci_list, dataset) {
  process_matrix <- function(matrix, method) {
    df <- as.data.frame(matrix)
    df$n_neighbors <- rownames(df)
    colnames(df) <- gsub(paste0(method, "_"), "", colnames(df))
    df_long <- pivot_longer(df, cols = c("raw", "umap", "pca"),
                            names_to = "type", values_to = "value")
    df_long$method <- method
    return(df_long)
  }

  viz_data <- process_matrix(cci_list[[1]], "viz")
  diffu_data <- process_matrix(cci_list[[2]], "diffu")
  all_data <- rbind(viz_data, diffu_data)
  all_data$type <- factor(all_data$type, levels = c("raw", "umap", "pca"))

  p <- ggplot(all_data, aes(x = .data$type, y = .data$value, color = .data$method)) +
    geom_boxplot() +
    scale_fill_manual(values = c("viz" = "#8669A9", "diffu" = "#4B8537"),
                      labels = c("viz" = "metaspec", "diffu" = "ADM")) +
    labs(x = " ", y = " ", fill = "Method", title = dataset) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )

  return(p)
}

#' Simplify labels based on a provided mapping
#'
#' This function simplifies a vector of labels using an optional mapping. If no mapping
#' is provided, it returns the original labels.
#'
#' @param info A vector of original labels to be simplified.
#' @param mapping An optional named vector where names are original labels and values are
#'   simplified labels. Default is NULL.
#'
#' @return A vector of simplified labels. If no mapping is provided, returns the original
#'   labels. If mapping is provided, returns mapped labels where possible, falling
#'   back to original labels where no mapping exists.
#'
#' @details
#' If a mapping is provided, the function replaces each label in info with its
#' corresponding value in the mapping. For any labels not found in the mapping,
#' the original label is retained.
#'
#' @examples
#' # Without mapping
#' original_labels <- c("Type A", "Type B", "Type C")
#' simplify_labels(original_labels)
#' # Returns: c("Type A", "Type B", "Type C")
#'
#' # With mapping
#' mapping <- c("Type A" = "A", "Type B" = "B")
#' simplify_labels(original_labels, mapping)
#' # Returns: c("A", "B", "Type C")
#'
#' @export
simplify_labels <- function(info, mapping = NULL) {
  if (is.null(mapping)) {
    return(info)
  } else {
    simplified <- mapping[as.character(info)]
    simplified[is.na(simplified)] <- info[is.na(simplified)]
    return(simplified)
  }
}

#' Visualize Individual Methods
#'
#' This function creates scatter plots for multiple matrices using different methods,
#' and calculates ARI, NMI, and silhouette scores for each method. It prints the
#' average silhouette width and the plot for each method during execution.
#'
#' @param matrices A list of matrices, each representing a different method's output.
#' @param method_names A character vector specifying the names of the visualization methods.
#' @param info An optional vector of group labels for each data point. If not provided, all points will be the same color.
#' @param color_list An optional named vector of colors for each group in 'info'. If not provided, a default color palette will be used.
#' @param k The number of clusters for k-means clustering.
#' @param seed An optional seed for reproducibility (default is 42).
#'
#' @return A list of results for each method. Each element of the list contains:
#'   \item{plot}{A ggplot object representing the scatter plot}
#'   \item{ari}{The Adjusted Rand Index}
#'   \item{nmi}{The Normalized Mutual Information}
#'   \item{silhouette}{The average silhouette width}
#'
#' @details This function will print the average silhouette width and display
#'          the plot for each method during execution.
#'
#' @examples
#' matrices <- list(matrix(rnorm(200), 100, 2), matrix(rnorm(200), 100, 2))
#' method_names <- c("Method1", "Method2")
#' info <- sample(c("A", "B", "C"), 100, replace = TRUE)
#' color_list <- c("A" = "red", "B" = "blue", "C" = "green")
#' results <- visualize_individual_methods(matrices, method_names, info, color_list, k = 3)
#'
#' @importFrom gridExtra grid.arrange
#' @importFrom cluster silhouette
#' @importFrom stats kmeans dist
#' @export

visualize_individual_methods <- function(matrices, method_names, info = NULL, color_list = NULL, k, seed = 42) {
  if (length(matrices) != length(method_names)) {
    stop("The number of matrices must match the number of method names.")
  }
  
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Package 'cluster' is needed for this function to work. Please install it.", call. = FALSE)
  }

  results <- lapply(seq_along(matrices), function(i) {
    ari_nmi <- cal_ari_nmi(matrices[[i]], k, method_names[i], seed, info)
    
    set.seed(seed)
    cluster_viz <- stats::kmeans(matrices[[i]], centers = k)
    sil <- cluster::silhouette(cluster_viz$cluster, dist(matrices[[i]]))
    avg_sil <- mean(sil[, 3])
    
    plot <- visualization_func(matrices[[i]], method_names[i], color_list, info)
    
    cat("Average Silhouette Width:", avg_sil, "\n\n")
    
    print(plot)
    
    list(
      plot = plot,
      ari = ari_nmi$ARI,
      nmi = ari_nmi$NMI,
      silhouette = avg_sil
    )
  })
  
  return(results)
}



#' Process and Visualize Meta-Methods
#'
#' @param ensemble.out Output from an ensemble method (assumed to contain $ensemble.dist.mat). Can be NULL.
#' @param mev.out Output from MEV method (assumed to contain $diffu.dist).
#' @param info A vector of true class labels.
#' @param k The number of clusters for k-means clustering.
#' @param color_list A named vector of colors for visualization, names should match unique values in info.
#' @param seed An integer seed for reproducibility. Default is 2024.
#'
#' @return A list containing:
#'   \item{ARI_list}{A list of data frames with ARI and NMI scores for each method}
#'   \item{ASW_list}{A list of average silhouette widths for each method}
#'   \item{umap_adm}{UMAP coordinates for ADM method}
#'   \item{plot}{A list of ggplot objects for each method}
#'   \item{umap_viz}{UMAP coordinates for meta-spec method (if ensemble.out is provided)}
#'
#' @details
#' This function processes and visualizes results from two meta-methods: meta-spec and ADM.
#' If ensemble.out is NULL, only ADM method is executed.
#' It performs UMAP dimensionality reduction, k-means clustering, and calculates various metrics.
#'
#' @note
#' This function requires the following packages: uwot, mclust, aricode, cluster
#'
#' @importFrom uwot umap
#' @importFrom stats dist kmeans
#' @importFrom cluster silhouette
#' @importFrom mclust adjustedRandIndex
#' @importFrom aricode NMI
#'
#' @examples
#' # This function requires specific input structures and external functions
#' # An example cannot be easily provided without those dependencies
#'
#' @export
process_and_visualize_meta_methods <- function(mev.out, ensemble.out = NULL, info, k, color_list, seed = 2024) {
  # Input validation
  if (!is.null(ensemble.out) && !all(c("ensemble.dist.mat") %in% names(ensemble.out))) {
    stop("When provided, ensemble.out must contain 'ensemble.dist.mat'")
  }
  if (!all(c("diffu.dist") %in% names(mev.out))) {
    stop("mev.out must contain 'diffu.dist'")
  }
  if (!is.null(ensemble.out) && length(info) != nrow(ensemble.out$ensemble.dist.mat)) {
    stop("Length of info must match the number of rows in ensemble.dist.mat")
  }

  set.seed(seed)
  ARI_list <- list()
  ASW_list <- list()
  plots <- list()
  print(paste("Running R version:", R.version$major, ".", R.version$minor, sep = ""))

  # Meta-spec method (only if ensemble.out is not NULL)
  if (!is.null(ensemble.out)) {
    method <- "meta-spec"
    umap_viz0 <- uwot::umap(ensemble.out$ensemble.dist.mat)
    ARI_list[[1]] <- cal_ari_nmi(umap_viz0, k, method, seed, info)
    ASW_list[[1]] <- cluster::silhouette(as.numeric(factor(info)), dist = stats::dist(umap_viz0))[, 3]
    umap_viz <- as.data.frame(umap_viz0)
    plots[[1]] = visualization_func(umap_viz, method, color_list, info)
  }

  # ADM method
  method <- "ADM"
  umap_adm0 <- uwot::umap(mev.out$diffu.dist)
  ARI_list[[length(ARI_list) + 1]] <- cal_ari_nmi(umap_adm0, k, method, seed, info)
  ASW_list[[length(ASW_list) + 1]] <- cluster::silhouette(as.numeric(factor(info)), dist = stats::dist(umap_adm0))[, 3]
  umap_adm <- as.data.frame(umap_adm0)
  plots[[length(ASW_list) + 1]] = visualization_func(umap_adm, method, color_list, info)

  result <- list(
    ARI_list = ARI_list,
    ASW_list = ASW_list,
    umap_adm = umap_adm0,
    plot = plots
  )

  if (!is.null(ensemble.out)) {
    result$umap_viz <- umap_viz0
  }

  return(result)
}


#' Visualize Silhouette Width Comparison
#'
#' This function creates a boxplot to compare Silhouette Width values across different clusters and methods.
#'
#' @param ASW_list A list of vectors containing Silhouette Width values for each method.
#' @param info A vector of cluster labels for each data point.
#' @param data_name A string specifying the name of the dataset (used for plot title).
#' @param label_mapping An optional named vector for mapping original cluster labels to simplified ones.
#'
#' @return A ggplot object representing a boxplot of Silhouette Width values.
#'
#' @examples
#' ASW_list <- list(Method1 = rnorm(100, 0.7, 0.1), Method2 = rnorm(100, 0.6, 0.1))
#' info <- sample(c("A", "B", "C"), 100, replace = TRUE)
#' label_mapping <- c("A" = "Cluster 1", "B" = "Cluster 2", "C" = "Cluster 3")
#' p <- visualize_silhouette_width(ASW_list, info, "Example Dataset", label_mapping)
#' if (interactive()) {
#'   print(p)
#' }
#'
#' @export
visualize_silhouette_width <- function(ASW_list, info, data_name, label_mapping = NULL) {
  # Create data frame
  data_df <- data.frame(Cluster = info)
  for (i in seq_along(ASW_list)) {
    data_df[[paste0("Method", i)]] <- ASW_list[[i]]
  }

  # Pivot data to long format
  data_long <- tidyr::pivot_longer(data_df,
                                   cols = starts_with("Method"),
                                   names_to = "Method",
                                   values_to = "SilhouetteWidth")

  # Apply label mapping if provided
  if (!is.null(label_mapping)) {
    data_long$SimplifiedCluster <- label_mapping[data_long$Cluster]
  } else {
    data_long$SimplifiedCluster <- data_long$Cluster
  }

  # Create the plot
  p <- ggplot(data_long, aes(x = .data$SimplifiedCluster, y = .data$SilhouetteWidth, fill = .data$Method)) +
    geom_boxplot(position = position_dodge(0.8)) +
    scale_fill_manual(values = c("#4B8537", "#8669A9")) +
    labs(title = paste("Silhouette Width Comparison -", data_name),
         x = "Cluster",
         y = "Silhouette Width") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right",
      legend.key.height = unit(1.2, "cm"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )

  return(p)
}
#' Visualize ARI and NMI Comparison with Bar Plot
#'
#' This function generates a bar plot to compare ARI and NMI metrics across methods.
#'
#' @param ARI_list A list containing one or two elements: ARI and NMI scores for ADM method, and optionally for meta-spec method.
#' @param dataset A string naming the dataset being visualized.
#'
#' @return Invisibly returns a ggplot object visualizing ARI and NMI scores.
#'
#' @details
#' The function adapts to whether one (ADM only) or two (meta-spec and ADM) methods are present in the input.
#' It creates a bar plot using ggplot2, with metrics on the x-axis and their values on the y-axis.
#' Bars are grouped and colored by the 'group' variable. The plot uses a minimal theme and predefined colors.
#'
#' @examples
#' # Example with only ADM method
#' ARI_list_single <- list(data.frame(ARI = 0.7, NMI = 0.8))
#' visualize_ari_nmi(ARI_list_single, "Dataset A")
#'
#' # Example with both methods
#' ARI_list_double <- list(
#'   data.frame(ARI = 0.8, NMI = 0.9),
#'   data.frame(ARI = 0.7, NMI = 0.8)
#' )
#' visualize_ari_nmi(ARI_list_double, "Dataset B")
#'
#' @export
visualize_ari_nmi <- function(ARI_list, dataset) {
  # Input validation
  if (!is.list(ARI_list) || length(ARI_list) < 1 || length(ARI_list) > 2) {
    stop("ARI_list must be a list with one or two elements")
  }
  if (!all(sapply(ARI_list, function(x) all(c("ARI", "NMI") %in% names(x))))) {
    stop("Each element of ARI_list must be a data frame with 'ARI' and 'NMI' columns")
  }
  if (!is.character(dataset) || length(dataset) != 1) {
    stop("dataset must be a single string")
  }

  # Determine if we have one or two methods
  has_two_methods <- length(ARI_list) == 2

  # Prepare data
  if (has_two_methods) {
    data_combined <- dplyr::bind_rows(
      cbind(ARI_list[[1]], group = "meta-spec"),
      cbind(ARI_list[[2]], group = "ADM")
    )
  } else {
    data_combined <- cbind(ARI_list[[1]], group = "ADM")
  }

  # Convert to long format
  data_long <- tidyr::pivot_longer(data_combined,
                                   cols = c("ARI", "NMI"),
                                   names_to = "Metric",
                                   values_to = "Value")

  # Create the plot
  plot <- ggplot(data_long, aes(x = .data$Metric, y = .data$Value, fill = .data$group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
    labs(title = dataset,
         x = NULL, y = "") +
    scale_fill_manual(values = c("#4B8537", "#8669A9")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      legend.position = "bottom"
    )

  # Adjust legend if only one method
  if (!has_two_methods) {
    plot <- plot + theme(legend.position = "none")
  }

  # Print the plot and return it invisibly
  print(plot)
  invisible(plot)
}




