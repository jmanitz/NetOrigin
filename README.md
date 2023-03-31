
## `NetOrigin` package

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

### Example: Effective Distance

``` r
data(delayGoe)

# compute effective distance
data(ptnGoe)
goenet <- igraph::as_adjacency_matrix(ptnGoe, sparse=FALSE)
p <- goenet/rowSums(goenet)
eff <- eff_dist(p)
```

    ## Computing the effective distance between 257 nodes:
    ##  1...................................................................................................
    ##  100...................................................................................................
    ##  200.........................................................done

``` r
# apply effective distance median source estimation
om <- origin(events=delayGoe[10,-c(1:2)], type='edm', distance=eff)
summary(om)
```

    ## Effective distance median origin estimation:
    ## 
    ## estimated node of origin 91: X.Gotthelf.Leimbach.Strasse 
    ## 
    ## auxiliary variables:
    ##        id          events            wmean             wvar       
    ##  Min.   :  1   Min.   : 0.0000   Min.   : 5.482   Min.   :0.3987  
    ##  1st Qu.: 65   1st Qu.: 0.0000   1st Qu.:21.572   1st Qu.:2.2761  
    ##  Median :129   Median : 0.0000   Median :27.345   Median :2.4050  
    ##  Mean   :129   Mean   : 0.6459   Mean   :26.948   Mean   :2.4989  
    ##  3rd Qu.:193   3rd Qu.: 0.0000   3rd Qu.:33.359   3rd Qu.:2.9986  
    ##  Max.   :257   Max.   :46.0000   Max.   :47.762   Max.   :6.2052  
    ##      mdist      
    ##  Min.   :14.34  
    ##  1st Qu.:20.75  
    ##  Median :24.23  
    ##  Mean   :24.92  
    ##  3rd Qu.:28.88  
    ##  Max.   :39.16

``` r
plot(om, 'mdist', start=1)
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
plot(om, 'wvar', start=1)
```

![](README_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
performance(om, start=1, graph=ptnGoe)
```

    ##                   start                         est  hitt rank spj dist
    ## 1 X.Adolf.Hoyer.Strasse X.Gotthelf.Leimbach.Strasse FALSE    2   2 1332

### Example: Backtracking

``` r
# backtracking origin estimation (Manitz et al., 2016)
ob <- origin(events=delayGoe[10,-c(1:2)], type='backtracking', graph=ptnGoe)
summary(ob)
```

    ## Backtracking origin estimation:
    ## 
    ## estimated node of origin 87: X.Gesundbrunnen 
    ## 
    ## auxiliary variables:
    ##        id          events            bcount       
    ##  Min.   :  1   Min.   : 0.0000   Min.   :0.00000  
    ##  1st Qu.: 65   1st Qu.: 0.0000   1st Qu.:0.00000  
    ##  Median :129   Median : 0.0000   Median :0.00000  
    ##  Mean   :129   Mean   : 0.6459   Mean   :0.03891  
    ##  3rd Qu.:193   3rd Qu.: 0.0000   3rd Qu.:0.00000  
    ##  Max.   :257   Max.   :46.0000   Max.   :3.00000

``` r
plot(ob, start=1)
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
performance(ob, start=1, graph=ptnGoe)
```

    ##                   start             est  hitt rank spj dist
    ## 1 X.Adolf.Hoyer.Strasse X.Gesundbrunnen FALSE    4   8 5328

### Example: Multiple Origins

``` r
data(ptnAth)
origin_multiple(events=delayAth[10,-c(1:2)], type='backtracking', graph=ptnAth, no=2)
```

    ## [[1]]
    ## Backtracking origin estimation:
    ## 
    ## estimated node of origin 6: 6 
    ## 
    ## [[2]]
    ## Backtracking origin estimation:
    ## 
    ## estimated node of origin 1: 1

``` r
# edm
athnet <- igraph::as_adjacency_matrix(ptnAth, sparse=FALSE)
p <- athnet/rowSums(athnet)
eff <- eff_dist(p)
```

    ## Computing the effective distance between 51 nodes:
    ##  1...................................................done

``` r
origin_multiple(events=delayAth[10,-c(1:2)], type='edm', graph=ptnAth, no=2, distance=eff)
```

    ## [[1]]
    ## Effective distance median origin estimation:
    ## 
    ## estimated node of origin 3: 3 
    ## 
    ## [[2]]
    ## Effective distance median origin estimation:
    ## 
    ## estimated node of origin 2: 2

### References

  - Li, J., J. Manitz, E. Bertuzzo, and E.D. Kolaczyk (2021):
    Sensor-based localization of epidemic sources on human mobility
    networks. PLoS Comput Biol 17(1): e1008545. \<DOI:
    10.1371/journal.pcbi.1008545\>

  - Manitz, J., J. Harbering, M. Schmidt, T. Kneib, and A. Schoebel
    (2017): Source Estimation for Propagation Processes on Complex
    Networks with an Application to Delays in Public Transportation
    Systems. Journal of Royal Statistical Society C (Applied
    Statistics), 66: 521–536. \<DOI: 10.1111/rssc.12176\>

  - Manitz, J., T. Kneib, M. Schlather, J. Helbing, and D. Brockmann
    (2014): Origin detection during food-borne disease outbreaks - a
    case study of the 2011 EHEC/HUS outbreak in Germany. PLoS Currents
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
