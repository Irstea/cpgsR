cpgsR: uniform sampling in a convex polytope
==================================================
---
License: LGPL-3

# Introduction #
This package is based on an initial matlab code called CPRND that can be found here
(https://ch.mathworks.com/matlabcentral/fileexchange/34208-uniform-distribution-over-a-convex-polytope

# Installation #
The easiest solution is to use the 'devtools' packages, this will installed the package and all required dependencies. To use 'devtools' , you need to install first:
* on Windows: [Rtools](http://cran.r-project.org/bin/windows/Rtools/)  
* on Mac: [Xcode command line tools](https://developer.apple.com/downloads)  
* on Linux: the R development package, usually called r-devel or r-base-dev  
  
To generate the vignette, you also need to have pandoc and pandoc-citeproc installed. Instructions can be found [here](https://pandoc.org/installing.html).    
  
Then, on a R console:

    > install.packages("devtools")
    > library(devtools)
    > install_github("Irstea/cpgsR",build=TRUE,build_opts = c("--no-resave-data","--no-manual"),build_vignettes = TRUE)

# Usage #
A vignette is included in the package to explain the usage.  

# Bug Reporting #
Please, report bug on the [github site](https://github.com/Irstea/cpgsR/issues).
