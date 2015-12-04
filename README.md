# SpatialPosition

![Stewart Potentials ](http://rgeomatic.hypotheses.org/files/2015/12/potentials.png)

R package for computing spatial position models:  

* Stewart potentials
* Reilly catchment areas
* Huff catchment areas



## Installation
### From CRAN
Stable version
```{r}
install.packages("SpatialPosition")
```

### From GitHub
Development version
```{r}
require(devtools)
devtools::install_github("Groupe-ElementR/SpatialPosition")
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
