# SpatialPosition

<!-- badges: start -->
[![Version](http://www.r-pkg.org/badges/version/SpatialPosition)](https://CRAN.R-project.org/package=SpatialPosition)
[![R-CMD-check](https://github.com/riatelab/SpatialPosition/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/riatelab/SpatialPosition/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

![Stewart Potentials ](http://rgeomatic.hypotheses.org/files/2015/12/potentials.png)

R package for computing spatial position models:  

* Stewart potentials :warning: 
* Reilly catchment areas
* Huff catchment areas

:warning: **Functions related to Stewart's potential are deprecated. Please use the [`potential`](https://riatelab.github.io/potential/) package instead.**


## Installation
### From CRAN
Stable version
```{r}
install.packages("SpatialPosition")
```

### From GitHub
Development version / unstable
```{r}
require(remotes)
remotes::install_github("riatelab/SpatialPosition")
```

### From Archive
Version 2.0.0 of the package is a major release that should not break old code.  
Anyway, in case of problem you can use these lines to get a previous version of the package. 
```{r}
packageurl <- "https://cran.r-project.org/src/contrib/Archive/SpatialPosition/SpatialPosition_1.2.0.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
```

## Demo
Vignettes contain commented scripts on how to use `SpatialPostion`.

* Introduction to the SpatialPosition package :
```{r}
vignette(topic = "SpatialPosition")
```

* Stewart Potentials: a Use Case :
```{r}
vignette(topic = "StewartExample")
```
