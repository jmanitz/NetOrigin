
## `NetOrigin` package

[![Build Status
(Linux)](https://travis-ci.org/jmanitz/NetOrigin.svg?branch=master)](https://travis-ci.org/jmanitz/NetOrigin)
[![Build status
(Windows)](https://ci.appveyor.com/api/projects/status/github/jmanitz/NetOrigin?branch=master&svg=true)](https://ci.appveyor.com/project/jmanitz/NetOrigin/branch/master)
[![Coverage
Status](https://coveralls.io/repos/github/jmanitz/NetOrigin/badge.svg?branch=master)](https://coveralls.io/github/jmanitz/NetOrigin?branch=master)
![GitHub repo
size](https://img.shields.io/github/repo-size/jmanitz/NetOrigin)
![GitHub issues](https://img.shields.io/github/issues/jmanitz/NetOrigin)

[![CRAN Status
Badge](http://www.r-pkg.org/badges/version/NetOrigin)](https://CRAN.R-project.org/package=NetOrigin)
[![CRAN
Downloads](http://cranlogs.r-pkg.org/badges/NetOrigin)](https://CRAN.R-project.org/package=NetOrigin)
![CRAN License](https://img.shields.io/cran/l/NetOrigin) ![CRAN
dependencies
status](https://img.shields.io/librariesio/release/CRAN/NetOrigin)
![Website](https://img.shields.io/website?url=http%3A%2F%2FNetOrigin.manitz.org%2F)

Performs network-based source estimation. Different approaches are
available: effective distance median, recursive backtracking, and
centrality-based source estimation. Additionally, we provide public
transportation network data as well as methods for data preparation,
source estimation performance analysis and visualization.

### Installation

You can install the latest production version from CRAN

``` r
install.packages("NetOrigin", dependencies = TRUE)
```

or the current development version from GitHub

``` r
library("devtools")
install_github("jmanitz/NetOrigin")
```

Then, load the package

``` r
library("NetOrigin")
```

### Example

<TODO>

### Contributions

<TODO>

### References

  - Manitz, J., J. Harbering, M. Schmidt, T. Kneib, and A. Schoebel (2017): 
  Source Estimation for Propagation Processes on Complex Networks with an 
  Application to Delays in Public Transportation Systems. 
  Journal of Royal Statistical Society C (Applied Statistics), 66: 521–536.

  - Manitz, J., Kneib, T., Schlather, M., Helbing, D. and Brockmann, D.
    (2014) Origin detection during food-borne disease outbreaks - a case
    study of the 2011 EHEC/HUS outbreak in Germany. PLoS Currents
    Outbreaks, 1. \<DOI:
    10.1371/currents.outbreaks.f3fdeb08c5b9de7c09ed9cbcef5f01f2\>

  - Comin, C. H. and da Fontoura Costa, L. (2011) Identifying the
    starting point of a spreading process in complex networks. Physical
    Review E, 84. \<DOI: 10.1103/PhysRevE.84.056105\>

To cite package ‘NetOrigin’ in publications use:

Juliane Manitz (2018). NetOrigin: Origin Estimation for Propagation
Processes on Complex Networks. R package version 1.0-3.
<https://CRAN.R-project.org/package=NetOrigin>

Use `toBibtex(citation("NetOrigin"))` in R to extract BibTeX references.
