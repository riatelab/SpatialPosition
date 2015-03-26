#' @title Potentials Package
#' @name Potentials
#' @description This package contains various functions.
#' @docType package
NULL

#' @title SpatialUnits
#' @name spatMask
#' @docType data
NULL

#' @title SpatialUnits
#' @name spatPts
#' @docType data
NULL

#' @title SpatialUnits
#' @name spatUnits
#' @docType data
NULL

#' @title SpatialUnits
#' @name FRdep
#' @docType data
NULL

#' @title SpatialUnits
#' @name FRdep.spdf
#' @docType data
NULL

#' @title SpatialUnits
#' @name spuniv
#' @docType data
NULL

#' @title Create a Regular SpatialPointsDataFrame
#' @name CreateGrid
#' @description This function creates a regular grid of SpatialPointDataFrame 
#' from the extent of a given Spatial*DataFrame and a given resolution.
#' @param w object of class sp (SpatialPointsDataFrame or 
#' SpatialPolygonsDataFrame). The spatial extent of this object is used to 
#' create the regular SpatialPointsDataFrame
#' @param resolution INTEGER, resolution of the output grid (in map units). 
#' @details The output of the function is a SpatialPointsDataFrame of regularly
#' spaced points with the same extent as \code{w}. 
#' @examples 
#' # Calculate potentials
#' data(spatData)
#' @import sp
#' @export
CreateGrid <- function (w, resolution)
{
  TestSp(w)
  boundingBox <- bbox(w)
  rounder <- boundingBox %% resolution
  boundingBox[,1] <- boundingBox[,1] - rounder[,1]
  roundermax <- resolution - rounder[,2]
  boundingBox[,2] <- boundingBox[,2] + resolution - rounder[,2]
  boxCoordX <- seq(from = boundingBox[1,1], to = boundingBox[1,2], 
                   by = resolution)
  boxCoordY <- seq(from = boundingBox[2,1], to = boundingBox[2,2], 
                   by = resolution)
  spatGrid <- expand.grid(boxCoordX, boxCoordY)
  idSeq <- seq(1, nrow(spatGrid), 1)
  spatGrid <- data.frame(ID = idSeq, 
                         COORDX = spatGrid[, 1], 
                         COORDY = spatGrid[, 2])
  
  spatGrid <- SpatialPointsDataFrame(coords = spatGrid[ , c(2, 3)], 
                                     data = spatGrid, 
                                     proj4string = CRS(proj4string(w)))
  return(spatGrid)
}


#' @title Create a Distance Matrix Between Two SpatialPointsDataFrame Objects
#' @name CreateDistMatrix
#' @description This function creates a distace matrix between two 
#' SpatialPointsDataFrame Objects
#' @param knownpts object of class SpatialPointsDataFrame. 
#' @param unknownpts object of class SpatialPointsDataFrame. 
#' @param longlat if FALSE, Euclidean distance, if TRUE Great Circle distance 
#' @details The function returns a full matrix of distances in the metric of the
#'  points if longlat=FALSE, or in kilometers if longlat=TRUE. This is a wrapper
#'   around the \code{\link{spDists}} function.
#' @examples 
#' # Calculate potentials
#' data(spatData)
#' @import sp
#' @export
CreateDistMatrix  <- function(knownpts, unknownpts, longlat = FALSE)
{
  TestSp(knownpts)
  TestSp(unknownpts)
  if(identicalCRS(knownpts,unknownpts) == FALSE){
    stop(paste("Inputs (",quote(knownpts), " and ",quote(unknownpts),
               ") do not use the same projection", sep = ""),call. = FALSE)
  }
  
  nk <- nrow(knownpts)
  nu <- nrow(unknownpts)
  if(nk * nu > 100000000 | nu > 10000000 | nk > 10000000){
    if (interactive()){
      cat("Do you really want to compute potentials values (from", nk , 
          "known points to", nu,"estimated values) ? \n 
          (Your computer is likely to collapse under the charge.) [y/n]" )
      z <- readLines(con = stdin(), n = 1) 
      while (!z %in% c("n","y")){
        cat ("Enter y or n")
        z <- readLines(con = stdin(), n = 1)  
      }
      if (z == "y"){
        cat("Ok, YOLO!")
      } else {
        stop("Computation aborted.", call. = F)
      }
    } else {
      stop("Computation aborted. Matrix would be too big", call. = F)
    }
  }
  
  matDist <- spDists(x = knownpts, y = unknownpts, longlat = longlat)
  dimnames(matDist) <- list(row.names(knownpts), row.names(unknownpts))
  return(round(matDist, digits = 8))
}


#' @title Calculate Stewart's Potentials
#' @name stewart
#' @description This function compute the potentials of spatial interaction as 
#' defined by J.Q. Stewart (1950).
#' @param knownpts object of class sp (SpatialPointsDataFrame or 
#' SpatialPolygonsDataFrame). It contains the known points to calculate 
#' potentials from.
#' @param unknownpts object of class sp (SpatialPointsDataFrame or 
#' SpatialPolygonsDataFrame). It contains points to calculate the potentials on. 
#' When used \code{resolution} is not used.
#' @param matdist NUMERIC MATRIX, a distance matrix. Row names match the first 
#' column of the knownpts object dataframe. Column names match the first column 
#' of the knownpts object dataframe. if used \code{longlat} is not used
#' @param varname CHARACTER, name of the variable (in the dataframe of the 
#' Spatial*DataFrame) from which potentials are computed.
#' Quantitative variable with no negative values. 
#' @param typefct CHARACTER, spatial interaction function. Options are "pareto" 
#' (default) or "exponential".
#' If "pareto" the interaction is defined as: (1 + alpha * mDistance) ^ (-beta).
#' If "exponential" the interaction is defined as: 
#' exp(- alpha * mDistance ^ beta).
#' The alpha parameter is computed from parameters given by the user 
#' (beta and reach).
#' @param span NUMERIC, distance where the density of probability of the spatial 
#' interaction function equals 0.5.
#' @param beta NUMERIC, impedance factor for the spatial interaction function.  
#' @param resolution INTEGER, resolution of the output SpatialPointsDataFrame
#'  (in map units). 
#' @param longlat BOOLEAN, type of distance. If TRUE longitudz and latitude are 
#' expected in decimal degrees, in WGS84 reference system.
#' @param mask object of class sp (SpatialPolygonsDataFrame) to clip the map of 
#' potentials.
#' @details The ouput value returned by the function is a SpatialPointsDataFrame 
#' of points with the computed potential in the "POTENTIALS" variable
#' @examples 
#' # Calculate potentials
#' data(spatData)
#' @references Stewart J. Q., Demographic gravitation: evidence and applications
#' , Sociometry, 11(1-2), 31-58.
#' @import sp
#' @import raster
#' @export
stewart <- function(knownpts,
                    unknownpts = NULL,
                    matdist = NULL,
                    varname,
                    typefct = "exponential", 
                    span,
                    beta,
                    resolution = 2000,
                    longlat = FALSE, 
                    mask = NULL)
{
  TestSp(knownpts)
  if (!is.null(unknownpts)){  
    TestSp(unknownpts)
    if(identicalCRS(knownpts,unknownpts) == FALSE){
      stop(paste("Inputs (",quote(knownpts), " and ",quote(unknownpts),
                 ") do not use the same projection", sep = ""),call. = FALSE)
    }
    if (!is.null(matdist)){
      matdist <- UseDistMatrix(matdist =matdist, knownpts = knownpts, 
                               unknownpts =  unknownpts) 
    }else{
      matdist <- CreateDistMatrix(knownpts = knownpts, unknownpts = unknownpts, 
                                  longlat = longlat) 
    }
  } else {
    unknownpts <- CreateGrid(w = if(is.null(mask)){knownpts} else {mask}, 
                             resolution = resolution) 
    matdist <- CreateDistMatrix(knownpts = knownpts, unknownpts = unknownpts, 
                                longlat = longlat) 
  }
  
  
  matdens <- ComputeInteractDensity(matdist = matdist, typefct = typefct,
                                    beta = beta, span = span)
  
  matopport <- ComputeOpportunity(knownpts = knownpts, matdens = matdens, 
                                  varname = varname)
  
  unknownpts <- ComputePotentials(unknownpts = unknownpts, 
                                  matopport = matopport)
  
  return(unknownpts)
}


#' @title Calculate Huff
#' @name huff
#' @description This function compute the potentials of spatial interaction as 
#' defined by J.Q. Stewart (1950).
#' @param knownpts object of class sp (SpatialPointsDataFrame or 
#' SpatialPolygonsDataFrame). It contains the known points to calculate 
#' potentials from.
#' @param unknownpts object of class sp (SpatialPointsDataFrame or 
#' SpatialPolygonsDataFrame). It contains points to calculate the potentials on. 
#' When used \code{resolution} is not used.
#' @param matdist NUMERIC MATRIX, a distance matrix. Row names match the first 
#' column of the knownpts object dataframe. Column names match the first column 
#' of the knownpts object dataframe. if used \code{typedist} is not used
#' @param varname CHARACTER, name of the variable (in the dataframe of the 
#' Spatial*DataFrame) from which potentials are computed.
#' Quantitative variable with no negative values. 
#' @param typefct CHARACTER, spatial interaction function. Options are "pareto" 
#' (default) or "exponential".
#' If "pareto" the interaction is defined as: (1 + alpha * mDistance) ^ (-beta).
#' If "exponential" the interaction is defined as: 
#' exp(- alpha * mDistance ^ beta).
#' The alpha parameter is computed from parameters given by the user 
#' (beta and reach).
#' @param span NUMERIC, distance where the density of probability of the spatial 
#' interaction function equals 0.5.
#' @param beta NUMERIC, impedance factor for the spatial interaction function.  
#' @param resolution INTEGER, resolution of the output grid (in map units). 
#' @param longlat BOOLEAN, type of distance. If TRUE longitudz and latitude are 
#' expected in decimal degrees, in WGS84 reference system.
#' @param mask object of class sp (SpatialPolygonsDataFrame) to clip the map of 
#' potentials.
#' @details The ouput value returned by the function is a SpatialPointsDataFrame 
#' of points with the computed potential in the "POTENTIALS" variable
#' @examples 
#' # Calculate potentials
#' data(spatData)
#' @references Stewart J. Q., Demographic gravitation: evidence and applications
#' , Sociometry, 11(1-2), 31-58.
#' @import sp
#' @import raster
#' @export
huff <- function(knownpts,
                 unknownpts = NULL,
                 matdist = NULL,
                 varname,
                 typefct = "exponential", 
                 span,
                 beta,
                 resolution = 2000,
                 longlat = FALSE, 
                 mask = NULL)
{
  TestSp(knownpts)
  if (!is.null(unknownpts)){  
    TestSp(unknownpts)
    if(identicalCRS(knownpts,unknownpts) == FALSE){
      stop(paste("Inputs (",quote(knownpts), " and ",quote(unknownpts),
                 ") do not use the same projection", sep = ""),call. = FALSE)
    }
    if (!is.null(matdist)){
      matdist <- UseDistMatrix(matdist =matdist, knownpts = knownpts, 
                               unknownpts =  unknownpts) 
    }else{
      matdist <- CreateDistMatrix(knownpts = knownpts, unknownpts = unknownpts, 
                                  longlat = longlat) 
    }
  } else {
    unknownpts <- CreateGrid(w = if(is.null(mask)){knownpts} else {mask}, 
                             resolution = resolution) 
    matdist <- CreateDistMatrix(knownpts = knownpts, unknownpts = unknownpts, 
                                longlat = longlat) 
  }
  
  
  matdens <- ComputeInteractDensity(matdist = matdist, typefct = typefct,
                                    beta = beta, span = span)
  
  matopport <- ComputeOpportunity(knownpts = knownpts, matdens = matdens, 
                                  varname = varname)
  
  unknownpts <- ComputeHuff(unknownpts = unknownpts, 
                            matopport = matopport)
  
  return(unknownpts)
}


#' @title Calculate reilly
#' @name reilly
#' @description This function compute the potentials of spatial interaction as 
#' defined by J.Q. Stewart (1950).
#' @param knownpts object of class sp (SpatialPointsDataFrame or 
#' SpatialPolygonsDataFrame). It contains the known points to calculate 
#' potentials from.
#' @param unknownpts object of class sp (SpatialPointsDataFrame or 
#' SpatialPolygonsDataFrame). It contains points to calculate the potentials on. 
#' When used \code{resolution} is not used.
#' @param matdist NUMERIC MATRIX, a distance matrix. Row names match the first 
#' column of the knownpts object dataframe. Column names match the first column 
#' of the knownpts object dataframe. if used \code{typedist} is not used
#' @param varname CHARACTER, name of the variable (in the dataframe of the 
#' Spatial*DataFrame) from which potentials are computed.
#' Quantitative variable with no negative values. 
#' @param typefct CHARACTER, spatial interaction function. Options are "pareto" 
#' (default) or "exponential".
#' If "pareto" the interaction is defined as: (1 + alpha * mDistance) ^ (-beta).
#' If "exponential" the interaction is defined as: 
#' exp(- alpha * mDistance ^ beta).
#' The alpha parameter is computed from parameters given by the user 
#' (beta and reach).
#' @param span NUMERIC, distance where the density of probability of the spatial 
#' interaction function equals 0.5.
#' @param beta NUMERIC, impedance factor for the spatial interaction function.  
#' @param resolution INTEGER, resolution of the output grid (in map units). 
#' @param longlat BOOLEAN, type of distance. If TRUE longitudz and latitude are 
#' expected in decimal degrees, in WGS84 reference system.
#' @param mask object of class sp (SpatialPolygonsDataFrame) to clip the map of 
#' potentials.
#' @details The ouput value returned by the function is a SpatialPointsDataFrame 
#' of points with the computed potential in the "POTENTIALS" variable
#' @examples 
#' # Calculate potentials
#' data(spatData)
#' @references Stewart J. Q., Demographic gravitation: evidence and applications
#' , Sociometry, 11(1-2), 31-58.
#' @import sp
#' @import raster
#' @export
reilly <- function(knownpts,
                   unknownpts = NULL,
                   matdist = NULL,
                   varname,
                   typefct = "exponential", 
                   span,
                   beta,
                   resolution = 2000,
                   longlat = FALSE, 
                   mask = NULL)
{
  TestSp(knownpts)
  if (!is.null(unknownpts)){  
    TestSp(unknownpts)
    if(identicalCRS(knownpts,unknownpts) == FALSE){
      stop(paste("Inputs (",quote(knownpts), " and ",quote(unknownpts),
                 ") do not use the same projection", sep = ""),call. = FALSE)
    }
    if (!is.null(matdist)){
      matdist <- UseDistMatrix(matdist =matdist, knownpts = knownpts, 
                               unknownpts =  unknownpts) 
    }else{
      matdist <- CreateDistMatrix(knownpts = knownpts, unknownpts = unknownpts, 
                                  longlat = longlat) 
    }
  } else {
    unknownpts <- CreateGrid(w = if(is.null(mask)){knownpts} else {mask}, 
                             resolution = resolution) 
    matdist <- CreateDistMatrix(knownpts = knownpts, unknownpts = unknownpts, 
                                longlat = longlat) 
  }
  
  
  matdens <- ComputeInteractDensity(matdist = matdist, typefct = typefct,
                                    beta = beta, span = span)
  
  matopport <- ComputeOpportunity(knownpts = knownpts, matdens = matdens, 
                                  varname = varname)
  
  unknownpts <- ComputeReilly(unknownpts = unknownpts, 
                              matopport = matopport)
  
  return(unknownpts)
}





#' @title Create a Raster from a SpatialPointsDataFrame Potential Layer
#' @name CreateRaster
#' @description This function create a raster layer from SpatialPointsDataFrame 
#' potential layer outputed by the \code{\link{stewart}} function. 
#' @param x object of class sp (SpatialPointsDataFrame) outputed by the
#'  stewart, huff or reilly function.
#' @param mask object of class sp (SpatialPolygonsDataFrame) to clip the map of 
#' potentials.
#' @details The ouput of the function is a raster of the potentials
#' @examples 
#' # Calculate potentials
#' data(spatData)
#' @import sp
#' @import raster
#' @export
CreateRaster <- function(x, mask = NULL){
  gridded(x) <- TRUE
  r <- raster(x)
  rasterx <- rasterize(x, r, field = 'OUTPUT')
  if(!is.null(mask)){
    TestSp(mask)
    rasterx <- mask(rasterx, mask = mask)
  }
  return(rasterx)
}



#' @title Plot a Potentials Raster
#' @name PlotPotentialsRaster 
#' @description This function plot a potential raster outputed by the 
#' \code{\link{CreateRaster}} function.
#' @param potentials.raster object of class raster outputed by the 
#' \code{\link{CreateRaster}} function.
#' @param add BOOLEAN, classic add parameter 
#' @param potentials.breaks NUMERIC, a vector of break values to map. If used, 
#' this parameter override \code{typec} and \code{nclass} params 
#' @param typec CHARACTER, either "equal" or "quantile", how to discretize the 
#' map.
#' @param nclass NUMERIC, number of class to map.
#' @param legend.rnd NUMERIC, number of digits used to round the figures 
#' displayed in the legend  
#' @param col.raster FUNCTION, colors to display the map like 
#' \code{\link{colorRampPalette}}. 
#' @details The ouput of the function is a plotted raster and the (invisible) 
#' list of break values
#' @examples 
#' # Calculate potentials
#' data(spatData)
#' 
#' @import sp
#' @import raster
#' @export
PlotPotentialsRaster <- function(potentials.raster, add = TRUE, 
                                 potentials.breaks = NULL, typec = "equal", 
                                 nclass = 5, legend.rnd = 0, 
                                 col.raster =  colorRampPalette(c("#FEA3A3", 
                                                                  "#980000"))){
  if (!is.null(potentials.breaks)){
    bks <- unique(potentials.breaks[order(potentials.breaks)])
  } else if (typec == "equal"){
    bks <- seq(from = cellStats(potentials.raster, min), 
               to = cellStats(potentials.raster, max), length.out = nclass+1)
  } else if (typec == "quantile"){
    bks <- quantile (potentials.raster, probs = seq(0,1, by = 1/nclass))
  } else {
    stop('Enter a proper discretisation type: "equal" or "quantile"')
  }
  bks <- unique(bks)
  col.raster <- col.raster(length(bks)-1)
  plot(potentials.raster, breaks = bks, legend = FALSE, axes = FALSE,
       box = FALSE, col = col.raster,  add = add)
  plot(potentials.raster, legend.only=TRUE, col = col.raster, 
       breaks=round(bks,legend.rnd))
  potentials.breaks <- bks
  return(invisible(potentials.breaks))
}

#' @title Plot Potentials Contour
#' @name PlotPotentialsContour
#' @description This function plot a contour from from SpatialPointsDataFrame 
#' potential layer outputed by the \code{\link{stewart}} function. 
#' @param potentials object of class sp (SpatialPointsDataFrame or 
#' SpatialPolygonsDataFrame) outputed by the potentials function.
#' @param potentials.breaks NUMERIC, a vector of break values to map. If used, 
#' this parameter override \code{typec} and \code{nclass} params 
#' @param legend.rnd NUMERIC, number of digits used to round the figures 
#' displayed in the legend  
#' @param col.contour CHARACTER, a color
#' @details This function returns only a contour on the current plot
#' @examples 
#' # Calculate potentials
#' data(spatData)
#' @import sp
#' @import raster
#' @export
PlotPotentialsContour <- function(potentials, potentials.breaks, 
                                  legend.rnd = 0, col.contour = "black"){
  contour.mat <- xtabs(POTENTIALS ~ COORDX + COORDY, data = potentials@data)
  x <- row.names(contour.mat)
  y <- colnames(contour.mat)
  contour(as.numeric(x), as.numeric(y), contour.mat, add=T, col = col.contour,
          levels = potentials.breaks[-1], 
          labels = round(potentials.breaks[-1], legend.rnd))
}

#' @title Create a SpatialLinesDataFrame of the Potential Raster Contour Lines
#' @name CreatePotentialsContourLines 
#' @description This function creates a SpatialLinesDataFrame from the potential 
#' raster outputed by the \code{\link{CreateRaster}} function.
#' @param potentials.raster object of class raster outputed by the 
#' \code{\link{CreateRaster}} function.
#' @param potentials.breaks NUMERIC, a vector of break values to map. 
#' @details The ouput of the function is a SpatialLinesDataFrame
#' @examples 
#' # Calculate potentials
#' data(spatData)
#' @import raster
#' @export
CreatePotentialsContourLines <- function(potentials.raster, potentials.breaks){
  rasterToContour(potentials.raster, levels = potentials.breaks)
}

#' @title Create a SpatialPolygonsDataFrame of the Potential Raster Contour Lines
#' @name CreatePotentialsContourPoly 
#' @description This function creates a SpatialPolygonsDataFrame from the 
#' potential raster outputed by the \code{\link{CreateRaster}} 
#' function.
#' @param potentials.raster object of class raster outputed by the 
#' \code{\link{CreateRaster}} function.
#' @param potentials.breaks NUMERIC, a vector of break values to map. 
#' @details The ouput of the function is a SpatialPolygonsDataFrame
#' @examples 
#' # Calculate potentials
#' data(spatData)
#' @import sp
#' @import raster
#' @export
CreatePotentialsContourPoly <- function(potentials.raster, potentials.breaks){
  if (!requireNamespace("rgeos", quietly = T)) {
    stop("'rgeos' package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  require(rgeos)
  cl <- rasterToContour(potentials.raster, levels = potentials.breaks)
  cl$level <- as.numeric (as.character(cl$level))
  SPlist <- NULL
  SPlist <- list()
  SPlevels <- character()
  for (i in cl$level){ 
    linex <- cl[cl@data$level == i,]
    linex <- linex@lines
    linex <- linex[[1]]
    linex <- linex@Lines
    Plist <- NULL
    Plist <- list()
    for (j in 1:length(linex)){
      x <- linex[[j]]@coords
      x <- Polygon(coords =  x, hole = F)
      x <- Polygons(srl = list(x), ID = j)
      Plist[j] <- x
    }  
    x <- SpatialPolygons(Srl = Plist)
    x <- rgeos::union(x = x)
    if (class(x) != "SpatialPolygonsDataFrame"){
      x <- SpatialPolygonsDataFrame(Sr = x, 
                                    data = data.frame(
                                      level = rep(i, length(x))))
    } else {
      x <- x[x@data$count < 2,]
      x@data <- data.frame(level = rep(i, dim(x)[1]))
    }
    SPlist <- c(SPlist , x@polygons  )
    SPlevels <- c(SPlevels,x@data$level)
  }
  
  for (i in 1:length(SPlist)){
    SPlist[[i]]@ID <- as.character(i)
  }
  
  x <- SpatialPolygons(Srl = SPlist)
  x <- SpatialPolygonsDataFrame(Sr = x, data = data.frame(levels = SPlevels))
  return (x)
}


# Internal functions
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
  matOpport <- knownpts@data[, varname] * matdens
  return(round(matOpport, digits = 8))
}

ComputePotentials <- function(unknownpts, matopport)
{
  unknownpts@data$OUTPUT <- apply(matopport, 2, sum, na.rm = TRUE)
  return(unknownpts)
}

ComputeReilly <- function(unknownpts, matopport)
{
  unknownpts@data$OUTPUT <- row.names(matopport)[apply(matopport, 2, which.max)]
  return(unknownpts)
}

ComputeHuff <- function(unknownpts, matopport)
{
  sumCol <- colSums(x = matopport, na.rm = TRUE)
  matOpportPct <- 100 * t(t(matopport) / sumCol)
  matOpportPct[is.na(matOpportPct) | is.infinite(matOpportPct)] <- 0
  unknownpts@data$OUTPUT <- apply(matOpportPct, 2, max, na.rm = TRUE)
  return(unknownpts)
}

TestSp <- function(x){
  if (substr(class(x),1,7) != "Spatial"){
    stop(paste("Your input (",quote(x),") is not a spatial object.", sep=""),
         call. = F)
  }
  if (is.na(x@proj4string)){
    stop(
      paste(
        "Your input (", quote(x),
        ") does not have a valid coordinate reference system.", sep=""),
      call. = F)
  }
}

