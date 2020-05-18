# REdiscoverTEdata
Author: *Haiyin Chen*

Date: *2019-09*

# Overview
This R package enables reproducibility of main figures from the paper TRANSPOSABLE ELEMENT EXPRESSION IN TUMORS IS ASSOCIATED WITH IMMUNE INFILTRATION AND INCREASED ANTIGENICITY by `Kong Y, ..., Chen-Harris H (2019)`.


# Contents

`REdiscoverTEdata` is an `R` package that can be installed locally. It contains source R Markdown (`.Rmd`) files for five main figures from the paper. Each `Rmd` file can generate one set of figures, either in HTML format or PDF format.

# Installation Process
  1. Install `R` as described at [https://www.r-project.org/](https://www.r-project.org/). `REdiscoverTEdata` has been tested with `R 3.4.3`, `R 3.5.1`.
  2. Install `RStudio`, which we recommend using to process .Rmd files into final figures. `RStudio` can be downloaded from [https://www.rstudio.com/](https://www.rstudio.com/).
  3. Install the `REdiscoverTEdata` package (library) locally with the following:
     * `cd /your/DOWNLOADS/directory/`
     * `tar -xvf REdiscoverTEdata.tar.gz` (uncompress the downloaded archive)
     * `R CMD INSTALL REdiscoverTEdata`
  * That command should automatically install the prerequisite libraries. You will find the `.Rmd` files for each figure in the `REdiscoverTEdata/inst/` subdirectory.

# List of R modules used by REdiscoverTE
  * Biobase
  * ComplexHeatmap
  * circlize
  * dplyr
  * DT
  * GenomicRanges
  * ggplot2
  * ggrepel
  * grid
  * gridExtra
  * knitr
  * multiGSEA (not a CRAN package, see below)
  * plyr
  * RColorBrewer
  * rmarkdown
  * stats

`multiGSEA` is available on github at
 https://github.com/lianos/multiGSEA (you may need to install this manually depending on the version of R and bioconductor you are using)

# Reproducing figures from the paper
   * Once `REdiscoverTEdata` is installed to your `R` library directory, you can launch `RStudio` and, from within `RStudio`, open any of the following `.Rmd` (R Markdown) files.
      * `REdiscoverTEdata/inst/Figure_1.Rmd`
      * `REdiscoverTEdata/inst/Figure_2.Rmd`
      * `REdiscoverTEdata/inst/Figure_3.Rmd`
      * `REdiscoverTEdata/inst/Figure_4.Rmd`
      * `REdiscoverTEdata/inst/Figure_5.Rmd`
   * Once an `.Rmd` file is open in `RStudio`, you can `knit` it (convert it to a figure) as follows:
      * `File Menu` -> `Knit Document`
      * alternatively, you can click on the "Knit" submenu and select `Knit to HTML` or `Knit to PDF`. 
   * The output from the `knit` operation will be an HTML (or PDF) file with the same name as the `.Rmd`, in the same directory. 
   
# TCGA TE expression matrix
   * TCGA TE expression matrix for 7000+ samples and 1000+ TE subfamilies can be found in the subdirectory `REdiscoverTEdata/inst/Fig4_data/` under the filename `eset_TCGA_TE_intergenic_logCPM.RDS`. The matrix is stored in the form of ExpressionSet, a data structure class defined in the `Biobase` package. 
   

