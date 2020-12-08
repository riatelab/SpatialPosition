#' @title Huff Catchment Areas
#' @name huff
#' @description This function computes the catchment areas as defined by D. Huff (1964).
#' @param knownpts sp or sf object; 
#' this is the set of known observations to estimate the catchment areas from.
#' @param unknownpts sp or sf object; 
#' this is the set of unknown units for which the function computes the estimates. 
#' Not used when \code{resolution} is set up. (optional)
#' @param matdist matrix; distance matrix between known observations and unknown 
#' units for which the function computes the estimates. Row names match the row 
#' names of \code{knownpts} and column names match the row names of 
#' \code{unknownpts}. \code{matdist} can contain any distance metric (time 
#' distance or euclidean distance for example). If \code{matdist} is not set, the distance 
#' matrix is automaticly built with \code{\link{CreateDistMatrix}}. (optional)
#' @param varname character; name of the variable in the \code{knownpts} dataframe 
#' from which values are computed. Quantitative variable with no negative values. 
#' @param typefct character; spatial interaction function. Options are "pareto" 
#' (means power law) or "exponential".
#' If "pareto" the interaction is defined as: (1 + alpha * mDistance) ^ (-beta).
#' If "exponential" the interaction is defined as: 
#' exp(- alpha * mDistance ^ beta).
#' The alpha parameter is computed from parameters given by the user 
#' (\code{beta} and \code{span}).
#' @param span numeric; distance where the density of probability of the spatial 
#' interaction function equals 0.5.
#' @param beta numeric; impedance factor for the spatial interaction function.  
#' @param resolution numeric; resolution of the output grid (in map units). 
#' If resolution is not set, the grid will contain around 7000 points. (optional)
#' @param mask sp or sf object; the spatial extent of this object is used to 
#' create the regularly spaced points output. (optional)
#' @param bypassctrl logical; bypass the distance matrix size control (see 
#' \code{\link{CreateDistMatrix}} Details).
#' @param longlat	logical; if FALSE, Euclidean distance, if TRUE Great Circle 
#' (WGS84 ellipsoid) distance.
#' @param returnclass "sp" or "sf"; class of the returned object.
#' @return Point object with the computed catchment areas in a new 
#' field named \code{OUTPUT}.
#' @seealso \link{huff}, \link{rasterHuff}, \link{plotHuff}, \link{CreateGrid}, 
#' \link{CreateDistMatrix}.
#' @examples
#' # Create a grid of paris extent and 200 meters
#' # resolution
#' data(hospital)
#' mygrid <- CreateGrid(w = paris, resolution = 200, returnclass = "sf")
#' # Create a distance matrix between known points (hospital) and mygrid
#' mymat <- CreateDistMatrix(knownpts = hospital, unknownpts = mygrid, 
#'                           longlat = FALSE)
#' # Compute Huff catchment areas from known points (hospital) on a given
#' # grid (mygrid) using a given distance matrix (mymat)
#' myhuff <- huff(knownpts = hospital, unknownpts = mygrid,
#'                matdist = mymat, varname = "capacity",
#'                typefct = "exponential", span = 1250,
#'                beta = 3, mask = paris, returnclass = "sf")
#' # Compute Huff catchment areas from known points (hospital) on a
#' # grid defined by its resolution
#' myhuff2 <- huff(knownpts = hospital, varname = "capacity",
#'                 typefct = "exponential", span = 1250, beta = 3,
#'                 resolution = 200, mask = paris, returnclass= "sf")
#' # The two methods have the same result
#' identical(myhuff, myhuff2)
#' # the function output an sf object
#' class(myhuff)
#' @references HUFF D. (1964) Defining and Estimating a Trading Area. Journal of Marketing, 28: 34-38.
#' @importFrom methods is as
#' @importFrom sf st_as_sf
#' @export
huff <- function(knownpts, unknownpts, matdist, varname,
                 typefct = "exponential", span, beta, resolution, mask, 
                 bypassctrl = FALSE, longlat = TRUE, returnclass="sp"){
  res <- prepdata(knownpts = knownpts, unknownpts = unknownpts, 
                  matdist = matdist, bypassctrl = bypassctrl, longlat = longlat,
                  mask = mask, resolution = resolution) 
  matdens <- ComputeInteractDensity(matdist = res$matdist, typefct = typefct,
                                    beta = beta, span = span)
  matopport <- ComputeOpportunity(knownpts = res$knownpts, matdens = matdens, 
                                  varname = varname)
  unknownpts <- ComputeHuff(unknownpts = res$unknownpts, matopport = matopport)
  if(returnclass=="sp"){unknownpts <- as(unknownpts, "Spatial")}
  return(unknownpts)
}

#' @title Create a Raster from a Huff SpatialPointsDataFrame
#' @name rasterHuff
#' @description This function creates a raster from a regularly spaced 
#' Huff grid (output of the \code{\link{huff}} function). 
#' @param x sp or sf object; output of the \code{huff} function.
#' @param mask sp or sf object; this object is used to clip 
#' the raster. (optional)
#' @return Raster of catchment areas values.
#' @seealso \link{huff}, \link{plotHuff}.
#' @examples
#' library(raster)
#' data(hospital)
#' # Compute Huff catchment areas from known points (hospital) on a
#' # grid defined by its resolution
#' myhuff <- huff(knownpts = hospital, varname = "capacity",
#'                typefct = "exponential", span = 750, beta = 2,
#'                resolution = 100, mask = paris, returnclass = "sf")
#' # Create a raster of huff values
#' myhuffraster <- rasterHuff(x = myhuff, mask = paris)
#' plot(myhuffraster)
#' @import sp
#' @import raster
#' @export
rasterHuff <- function(x, mask = NULL){
  if(is(x, "sf")){x <-   suppressWarnings(as(x, "Spatial"))}
  gridded(x) <- TRUE
  r <- raster(x)
  rasterx <- rasterize(x, r, field = 'OUTPUT')
  if(!is.null(mask)){
    if(is(mask, "sf")){mask <-   suppressWarnings(as(mask, "Spatial"))}
    projError(x, mask)
    rasterx <- mask(rasterx, mask = mask)
  }
  return(rasterx)
}

#' @title Plot a Huff Raster
#' @name plotHuff
#' @description This function plots the raster produced by the 
#' \code{\link{rasterHuff}} function.
#' @param x raster; output of the \code{\link{rasterHuff}} function.
#' @param add logical; if TRUE the raster is added to the current plot, if FALSE 
#' the raster is displayed in a new plot.
#' @return Display the raster nicely.
#' @seealso \link{huff}, \link{rasterHuff}.
#' @examples
#' data(hospital)
#' # Compute Huff catchment areas from known points (hospital) on a
#' # grid defined by its resolution
#' myhuff <- huff(knownpts = hospital, varname = "capacity",
#'                typefct = "exponential", span = 750, beta = 2,
#'                resolution = 100, mask = paris, returnclass = "sf")
#' # Create a raster of huff values
#' myhuffraster <- rasterHuff(x = myhuff, mask = paris)
#' plotHuff(myhuffraster)
#' @import raster
#' @import sp
#' @export
plotHuff <- function(x, add = FALSE){
  bks <- seq(from = cellStats(x, min), 
             to = cellStats(x, max), length.out = 11)
  col <- c("#543005", "#8C510A", "#BF812D", "#DFC27D", 
           "#F6E8C3","#C7EAE5", "#80CDC1", "#35978F", 
           "#01665E", "#003C30")
  plot(x, breaks = bks, legend = FALSE, axes = FALSE,
       box = FALSE, col = col,  add = add)
  plot(x, legend.only=TRUE, col = col, 
       breaks=round(bks, 0) )
}




