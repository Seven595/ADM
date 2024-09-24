# ADMPackage: Adaptive Graph Diffusion for Meta-Dimension Reduction

## Overview

ADMPackage is an R package that implements Adaptive Graph Diffusion for Meta-Dimension Reduction (ADM). This method performs meta-analysis on a list of dimension reduction results, providing a robust approach to visualize and analyze high-dimensional data in lower-dimensional spaces. It offers tools for comparing different dimensionality reduction methods, calculating cluster conservation indices, and evaluating cluster quality.

## Installation

You can install the development version of ADMPackage from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("Seven595/ADMPackage")
```

## Features

- Adaptive Graph Diffusion for Meta-Dimension Reduction (ADM)
- Calculation of Cluster Conservation Index (CCI)
- Support for multiple dimensionality reduction methods including PCA, MDS, UMAP, t-SNE, and more
- Visualization tools for comparing different methods
- Cluster analysis and evaluation metrics

## Usage

Here's a basic example of how to use ADMPackage:

```r
library(ADMPackage)

# Assuming you have your data in 'my_data' and cluster information in 'my_clusters'
result <- mev(my_data, my_clusters)

# Calculate CCI
cci_result <- cal_cci(result$ensemble.out, result$adm.out)

# Visualize results
plot_adm(result)
```

## Main Functions

- `adm()`: Performs Adaptive Graph Diffusion for Meta-Dimension Reduction
- `cal_cci()`: Calculates Cluster Conservation Index
- `visualization_func()`: Creates scatter plots of dimensionality reduction results
- `plot_bar_func()`: Creates bar plots for comparing metrics across methods

## Data

The package includes functions to work with various datasets, including:

- Oihane dataset
- Quake dataset
- miR dataset
- Gene expression dataset

Use the `get_mapping()` function to get the appropriate label mapping for each dataset.