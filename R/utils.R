# Internal functions

prepdata <- function(knownpts, unknownpts, matdist, bypassctrl, longlat, mask, 
                     resolution){
  if(is(knownpts, "Spatial")){knownpts <- st_as_sf(knownpts)}
  if (!missing(unknownpts)){  
    if(is(unknownpts, "Spatial")){unknownpts <- st_as_sf(unknownpts)}
    
    if (!missing(matdist)){
      matdist <- UseDistMatrix(matdist = matdist, knownpts = knownpts, 
                               unknownpts =  unknownpts) 
    }else{
      matdist <- CreateDistMatrix(knownpts = knownpts, unknownpts = unknownpts, 
                                  bypassctrl = bypassctrl, longlat = longlat)
    }
  }else{
    if(missing(mask)){
      mask <- knownpts
    } else {
      if(is(mask, "Spatial")){unknownpts <- st_as_sf(mask)}
    }
    unknownpts <- CreateGrid(w = mask, resolution = resolution, 
                             returnclass = "sf") 
    matdist <- CreateDistMatrix(knownpts = knownpts, unknownpts = unknownpts, 
                                bypassctrl = bypassctrl, longlat = longlat) 
  }
  return(list(knownpts=knownpts, unknownpts = unknownpts, matdist = matdist))
}


UseDistMatrix <- function(matdist, knownpts, unknownpts){
  i <- factor(row.names(knownpts), levels = row.names(knownpts))
  j <- factor(row.names(unknownpts), levels = row.names(unknownpts))
  matdist <- matdist[levels(i), levels(j)]
  return(round(matdist, digits = 8))
}

ComputeInteractDensity <- function(matdist, typefct, beta, span)
{
  if(typefct == "pareto") {
    alpha  <- (2 ^ (1 / beta) - 1) / span
    matDens <- (1 + alpha * matdist) ^ (-beta)
  } else if(typefct == "exponential") {
    alpha  <- log(2) / span ^ beta
    matDens <- exp(- alpha * matdist ^ beta)
  } else {
    stop("Please choose a valid interaction function argument (typefct)")
  }
  matDens <- round(matDens, digits = 8)
  return(matDens)
}

ComputeOpportunity <- function(knownpts, matdens, varname = varname)
{
  matOpport <- knownpts[[varname]] * matdens
  return(round(matOpport, digits = 8))
}

ComputePotentials <- function(unknownpts, matopport)
{
  unknownpts$OUTPUT <- apply(matopport, 2, sum, na.rm = TRUE)
  return(unknownpts)
}

ComputeReilly <- function(unknownpts, matopport)
{
  unknownpts$OUTPUT <- row.names(matopport)[apply(matopport, 2, which.max)]
  return(unknownpts)
}

ComputeHuff <- function(unknownpts, matopport)
{
  sumCol <- colSums(x = matopport, na.rm = TRUE)
  matOpportPct <- 100 * t(t(matopport) / sumCol)
  matOpportPct[is.na(matOpportPct) | is.infinite(matOpportPct)] <- 0
  unknownpts$OUTPUT <- apply(matOpportPct, 2, max, na.rm = TRUE)
  return(unknownpts)
}


ComputeSmooth<- function(unknownpts, matopport, matdens)
{
  unknownpts$OUTPUT <- apply(matopport, 2, sum, na.rm = TRUE) / 
    colSums(matdens, na.rm = TRUE)
  return(unknownpts)
}


projError <- function(x,y){
  if (is.na(x@proj4string)){
    stop("Your input does not have a valid coordinate reference system.",
         call. = F)
  }
  if(!missing(y)){
    if (is.na(y@proj4string)){
      stop("Your input does not have a valid coordinate reference system.",
           call. = F)
    }
    if(identicalCRS(x,y) == FALSE){
      stop("Inputs do not use the same coordinate reference system.",
           call. = FALSE)
    }
  }
}


