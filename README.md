# BE-Resistome
# AIRWAY 'RESISTOTYPES' AND CLINICAL OUTCOMES IN BRONCHIECTASIS

Accompanying code repository for the scientific manuscript "AIRWAY 'RESISTOTYPES' AND CLINICAL OUTCOMES IN BRONCHIECTASIS." by Mac Aogáin & Ivan et al. 2023 (under review).

All associated data can be found in the folder [Data](./Data/), while the [Analysis](./Analysis/) folder contains R code required to re-run analyses. Here, analysis supporting the results documented in the paper are presented as R-markdown files with accompanying rendered HTML files, are found [Here](./Data/R_output_files/knit_RMD_HTML) and can be downloaded and viewed for ease of review.

Summary of key R packages / Pyton modules used in the analysis: 
| Method           | Packages/modules  | Version | Description                                                                 | Main Functions Used                                          |
|------------------|-------------------|---------|-----------------------------------------------------------------------------|--------------------------------------------------------------|
| Data Wrangling   | tidyverse         | 2.0.0   | An opinionated collection of R packages designed for data science.           | `read_csv`, `as_tibble`, `filter`, `select`, `mutate`, `arrange`, `gather`, `pivot_longer` |
|                  | phyloseq          | 1.42.0  | Handling and analysis of high-throughput microbiome census data.            | `transform_sample_counts`, `otu_table`, `tax_table`, `readRDS` |
|                  | plyr              | 1.8.8   | Tools for splitting, applying, and combining data.                          | `ddply`                                                       |
| Visualization    | ggplot2           | 3.4.2   | A system for declaratively creating graphics, based on The Grammar of Graphics. | `ggplot`, `geom_bar`, `geom_point`, `geom_segment`, `facet_wrap`, `theme`, `guides` |
|                  | pheatmap          | 1.0.12  | Pretty heatmaps with more sensible behavior than the default heatmap function. | `pheatmap`                                                    |
|                  | colorspace        | 2.1-0   | A toolbox for manipulating and assessing colors and palettes.               | `sequential_hcl` (in heatmap color specification)            |
|                  | RColorBrewer      | 1.1-3   | Provides color schemes for maps and other graphics.                         | `brewer.pal`, `brewer.pal.info`                               |
|                  | ggpubr            | 0.6.0   | Provides some easy-to-use functions for creating and customizing 'ggplot2'-based publication ready plots. | `ggarrange`                                                   |
| Statistical Analysis | vegan          | 2.6-4   | Community ecology package with tools for diversity analysis, ordination, and environmental fits. | `adonis2`, `vegdist`                                          |
|                  | vcd               | 1.4-11  | Visualization techniques, data sets, summary and inference procedures for categorical data. | `associates`                                                  |
| Others           | pacman            | 0.5.1   | Tools to more conveniently manage the installation and loading of packages. | `p_load`                                                      |
|                  | reticulate        | 1.32.0  | Interface to 'Python' modules, classes, and functions.                      |                                                              |
|                  | knitr             | 1.44    | Provides a general-purpose tool for dynamic report generation in R.         |                                                              |
|                  | scikit-learn (python) | 1.3.0 | A machine learning library used for various data mining and data analysis. In this context, it's used for spectral clustering. | `spectralClustering`, `silhouette_score`                     |
|                  | numpy (python)    | 1.24.3  | A fundamental package for scientific computing with Python. It's used for numerical operations and array manipulations in the silhouette score calculation. | Array and Mathematical Operations                            |

