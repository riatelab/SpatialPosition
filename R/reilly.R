#' @title Reilly Catchment Areas
#' @name reilly
#' @description This function computes the catchment areas as defined by W.J. Reilly (1931).
#' @param knownpts sp or sf object; 
#' this is the set of known observations to estimate the catchment areas from.
#' @param unknownpts sp or sf object; 
#' this is the set of unknown units for which the function computes the estimates. 
#' Not used when \code{resolution} is set up. (optional)
#' @param matdist matrix; distance matrix between known observations and unknown 
#' units for which the function computes the estimates. Row names match the row 
#' names of \code{knownpts} and column names match the row names of 
#' \code{unknownpts}. \code{matdist} can contain any distance metric (time 
#' distance or euclidean distance for example). If \code{matdist} is not set, 
#' the distance matrix is built with \code{\link{CreateDistMatrix}}. (optional)
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
#' If resolution is not set, the grid will contain around 7250 points. (optional)
#' @param mask sp or sf object; the spatial extent of this object is used to 
#' create the regularly spaced points output. (optional)
#' @param bypassctrl logical; bypass the distance matrix size control (see 
#' \code{\link{CreateDistMatrix}} Details).
#' @param longlat	logical; if FALSE, Euclidean distance, if TRUE Great Circle 
#' (WGS84 ellipsoid) distance.
#' @param returnclass "sp" or "sf"; class of the returned object.
#' @return Point object with the computed catchment areas in a new 
#' field named \code{OUTPUT}. Values match the row names of \code{knownpts}.
#' @seealso \link{reilly}, \link{rasterReilly}, \link{plotReilly}, \link{CreateGrid}, 
#' \link{CreateDistMatrix}.
#' @examples 
#' # Create a grid of paris extent and 200 meters
#' # resolution
#' data(hospital)
#' mygrid <- CreateGrid(w = hospital, resolution = 200, returnclass = "sf")
#' # Create a distance matrix between known points (hospital) and mygrid
#' mymat <- CreateDistMatrix(knownpts = hospital, unknownpts = mygrid)
#' # Compute Reilly catchment areas from known points (hospital) on a given
#' # grid (mygrid) using a given distance matrix (mymat)
#' myreilly2 <- reilly(knownpts = hospital, unknownpts = mygrid,
#'                     matdist = mymat, varname = "capacity",
#'                     typefct = "exponential", span = 1250,
#'                     beta = 3, mask = paris, returnclass = "sf")
#' # Compute Reilly catchment areas from known points (hospital) on a
#' # grid defined by its resolution
#' myreilly <- reilly(knownpts = hospital, varname = "capacity",
#'                    typefct = "exponential", span = 1250, beta = 3,
#'                    resolution = 200, mask = paris, returnclass = "sf")
#' # The function output an sf object
#' class(myreilly)
#' # The OUTPUT field values match knownpts row names
#' head(unique(myreilly$OUTPUT))
#' @references REILLY, W. J. (1931) The law of retail gravitation, W. J. Reilly, New York.
#' @import sp
#' @import raster
#' @export
reilly <- function(knownpts, unknownpts, matdist, varname,
                   typefct = "exponential", span, beta, resolution, mask,
                   bypassctrl = FALSE, longlat = TRUE, returnclass="sp"){
  res <- prepdata(knownpts = knownpts, unknownpts = unknownpts, 
                  matdist = matdist, bypassctrl = bypassctrl, longlat = longlat,
                  mask = mask, resolution = resolution) 
  matdens <- ComputeInteractDensity(matdist = res$matdist, typefct = typefct,
                                    beta = beta, span = span)
  matopport <- ComputeOpportunity(knownpts = res$knownpts, matdens = matdens, 
                                  varname = varname)
  unknownpts <- ComputeReilly(unknownpts = res$unknownpts, 
                              matopport = matopport)
  if(returnclass=="sp"){unknownpts <- as(unknownpts, "Spatial")}
  return(unknownpts)
}

#' @title Create a Raster from a Reilly Regular Grid
#' @name rasterReilly
#' @description This function creates a raster from a regularly spaced 
#' Reilly grid (output of the \code{\link{reilly}} function). 
#' @param x sp or sf object; output of the \code{reilly} function.
#' @param mask sp or sf object; this object is used to clip 
#' the raster. (optional)
#' @return Raster of catchment areas values.
#' The raster uses a RAT (\code{\link{ratify}}) that contains the 
#' correspondance between raster values and catchement areas values. Use \code{
#' unique(levels(rasterName)[[1]])} to see the correpondance table.
#' @seealso \link{reilly}, \link{plotReilly}.
#' @examples
#' library(raster)
#' data(hospital)
#' # Compute Reilly catchment areas from known points (hospital) on a
#' # grid defined by its resolution
#' myreilly <- reilly(knownpts = hospital, varname = "capacity",
#'                    typefct = "exponential", span = 1250, beta = 3,
#'                    resolution = 200, mask = paris, returnclass = "sf")
#' # Create a raster of reilly values
#' myreillyraster <- rasterReilly(x = myreilly, mask = paris)
#' plot(myreillyraster, col = rainbow(18))
#' # Correspondance between raster values and reilly areas
#' head(unique(levels(myreillyraster)[[1]]))
#' @import sp
#' @import raster
#' @export
rasterReilly <- function(x ,mask = NULL){
  if(is(x, "sf")){x <- suppressWarnings(as(x, "Spatial"))}
  gridded(x) <- TRUE
  r <- raster(x)
  x$OUTPUT2 <- as.factor(x$OUTPUT)
  levels(x$OUTPUT2) <- 1:length(levels(x$OUTPUT2) )
  x$OUTPUT2 <- as.numeric(x$OUTPUT2)
  rasterx <- rasterize(x, r, field = 'OUTPUT2')
  if(!is.null(mask)){
    if(is(mask, "sf")){mask <- suppressWarnings(as(mask, "Spatial"))}
    projError(x, mask)
    rasterx <- mask(rasterx, mask = mask)
  }
  levels(rasterx) <- data.frame(ID = x$OUTPUT2, idarea = x$OUTPUT)
  return(rasterx)
}

#' @title Plot a Reilly Raster
#' @name plotReilly
#' @description This function plots the raster produced by the 
#' \code{\link{rasterReilly}} function.
#' @param x raster; output of the \code{\link{rasterReilly}} function.
#' @param add logical; if TRUE the raster is added to the current plot, if FALSE 
#' the raster is displayed in a new plot.
#' @param col function; color ramp function, such as \code{\link{colorRampPalette}}.
#' @details Display the raster nicely.
#' @seealso \link{reilly}, \link{rasterReilly}.
#' @examples
#' data(hospital)
#' # Compute Reilly catchment areas from known points (hospital) on a
#' # grid defined by its resolution
#' myreilly <- reilly(knownpts = hospital, varname = "capacity",
#'                    typefct = "exponential", span = 1250, beta = 3,
#'                    resolution = 200, mask = paris, returnclass = 'sf')
#' # Create a raster of reilly values
#' myreillyraster <- rasterReilly(x = myreilly, mask = paris)
#' # Plot the raster nicely
#' plotReilly(x = myreillyraster)
#' @import sp
#' @import raster
#' @importFrom grDevices rainbow
#' @export
plotReilly <- function(x, add = FALSE, 
                       col = rainbow){
  nclass <- nrow(unique(levels(x)[[1]]))
  colorReilly <- col(n = nclass)
  plot(x, legend = FALSE, axes = FALSE,
       box = FALSE, col = colorReilly,  add = add)
}