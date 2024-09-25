pkgname <- "ADM"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('ADM')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("cal_ari_nmi")
### * cal_ari_nmi

flush(stderr()); flush(stdout())

### Name: cal_ari_nmi
### Title: Calculate Adjusted Rand Index (ARI) and Normalized Mutual
###   Information (NMI)
### Aliases: cal_ari_nmi

### ** Examples

data <- matrix(rnorm(200), 100, 2)
info <- sample(1:3, 100, replace = TRUE)  # Simulated cluster information
result <- cal_ari_nmi(data, k = 3, method_name = "Example Method", seed = 123, info = info)




cleanEx()
nameEx("candidate.visual")
### * candidate.visual

flush(stderr()); flush(stdout())

### Name: candidate.visual
### Title: Perform Multiple Dimensionality Reduction Methods
### Aliases: candidate.visual

### ** Examples

## Not run: 
##D # Assuming 'data' is your input dataset
##D result <- candidate.visual(data, dim = 2,
##D                            methods = c("PCA", "MDS", "UMAP"),
##D                            umap.k = c(15, 30))
##D pca_result <- result$embed.list[[1]]
##D method_names <- result$method_name
## End(Not run)




cleanEx()
nameEx("dataloader")
### * dataloader

flush(stderr()); flush(stdout())

### Name: dataloader
### Title: Load and Prepare Dataset for Analysis
### Aliases: dataloader

### ** Examples

## Not run: 
##D # Load the Oihane dataset
##D oihane_data <- dataloader("Oihane")
##D 
##D # Load the PBMC dataset
##D pbmc_data <- dataloader("pbmc")
##D 
##D # Load the miRNA dataset
##D mir_data <- dataloader("mir")
## End(Not run)




cleanEx()
nameEx("get_mapping")
### * get_mapping

flush(stderr()); flush(stdout())

### Name: get_mapping
### Title: Get label mapping for a dataset
### Aliases: get_mapping

### ** Examples

get_mapping("Oihane")
get_mapping("Quake")

get_mapping("Oihane")
get_mapping("Quake")




cleanEx()
nameEx("plot_cci_results")
### * plot_cci_results

flush(stderr()); flush(stdout())

### Name: plot_cci_results
### Title: Plot Comparison of CCI Results
### Aliases: plot_cci_results

### ** Examples

viz_matrix <- matrix(rnorm(30), nrow = 10,
                     dimnames = list(NULL, c("viz_raw", "viz_umap", "viz_pca")))
diffu_matrix <- matrix(rnorm(30), nrow = 10,
                       dimnames = list(NULL, c("diffu_raw", "diffu_umap", "diffu_pca")))
cci_list <- list(viz_matrix, diffu_matrix)
p <- plot_cci_results(cci_list, "Example Dataset")
if (interactive()) {
  print(p)
}




cleanEx()
nameEx("plot_legend")
### * plot_legend

flush(stderr()); flush(stdout())

### Name: plot_legend
### Title: Create a scatter plot with a custom legend
### Aliases: plot_legend

### ** Examples

umap_coords <- matrix(rnorm(200), ncol = 2)
cluster_info <- rep(c("Type A", "Type B"), each = 50)
color_list <- c("Type A" = "red", "Type B" = "blue")
label_map <- c("Type A" = "A", "Type B" = "B")
p <- plot_legend(umap_coords, "UMAP", cluster_info, "Example", color_list, label_map)
if (interactive()) {
  print(p)
}




cleanEx()
nameEx("process_and_visualize_meta_methods")
### * process_and_visualize_meta_methods

flush(stderr()); flush(stdout())

### Name: process_and_visualize_meta_methods
### Title: Process and Visualize Meta-Methods
### Aliases: process_and_visualize_meta_methods

### ** Examples

# This function requires specific input structures and external functions
# An example cannot be easily provided without those dependencies




cleanEx()
nameEx("simplify_labels")
### * simplify_labels

flush(stderr()); flush(stdout())

### Name: simplify_labels
### Title: Simplify labels based on a provided mapping
### Aliases: simplify_labels

### ** Examples

# Without mapping
original_labels <- c("Type A", "Type B", "Type C")
simplify_labels(original_labels)
# Returns: c("Type A", "Type B", "Type C")

# With mapping
mapping <- c("Type A" = "A", "Type B" = "B")
simplify_labels(original_labels, mapping)
# Returns: c("A", "B", "Type C")




cleanEx()
nameEx("visualization_with_label")
### * visualization_with_label

flush(stderr()); flush(stdout())

### Name: visualization_with_label
### Title: Create a labeled scatter plot from UMAP visualization
### Aliases: visualization_with_label

### ** Examples

umap_coords <- matrix(rnorm(200), ncol = 2)
cluster_info <- rep(c("A", "B"), each = 50)
color_list <- c("A" = "red", "B" = "blue")
p <- visualization_with_label(umap_coords, "UMAP", cluster_info, "Example", color_list)
if (interactive()) {
  print(p)
}




cleanEx()
nameEx("visualize_ari_nmi")
### * visualize_ari_nmi

flush(stderr()); flush(stdout())

### Name: visualize_ari_nmi
### Title: Visualize ARI and NMI Comparison with Bar Plot
### Aliases: visualize_ari_nmi

### ** Examples

# Example with only ADM method
ARI_list_single <- list(data.frame(ARI = 0.7, NMI = 0.8))
visualize_ari_nmi(ARI_list_single, "Dataset A")

# Example with both methods
ARI_list_double <- list(
  data.frame(ARI = 0.8, NMI = 0.9),
  data.frame(ARI = 0.7, NMI = 0.8)
)
visualize_ari_nmi(ARI_list_double, "Dataset B")




cleanEx()
nameEx("visualize_individual_methods")
### * visualize_individual_methods

flush(stderr()); flush(stdout())

### Name: visualize_individual_methods
### Title: Visualize Individual Methods
### Aliases: visualize_individual_methods

### ** Examples

matrices <- list(matrix(rnorm(200), 100, 2), matrix(rnorm(200), 100, 2))
method_names <- c("Method1", "Method2")
# Example with default coloring (all points same color)
plots1 <- visualize_individual_methods(matrices, method_names)

# Example with group information
info <- sample(c("A", "B", "C"), 100, replace = TRUE)
color_list <- c("A" = "red", "B" = "blue", "C" = "green")
plots2 <- visualize_individual_methods(matrices, method_names, info, color_list)




cleanEx()
nameEx("visualize_silhouette_width")
### * visualize_silhouette_width

flush(stderr()); flush(stdout())

### Name: visualize_silhouette_width
### Title: Visualize Silhouette Width Comparison
### Aliases: visualize_silhouette_width

### ** Examples

ASW_list <- list(Method1 = rnorm(100, 0.7, 0.1), Method2 = rnorm(100, 0.6, 0.1))
info <- sample(c("A", "B", "C"), 100, replace = TRUE)
label_mapping <- c("A" = "Cluster 1", "B" = "Cluster 2", "C" = "Cluster 3")
p <- visualize_silhouette_width(ASW_list, info, "Example Dataset", label_mapping)
if (interactive()) {
  print(p)
}




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
