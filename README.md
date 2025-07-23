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

```{r}
install.packages("SpatialPosition")
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
