#' @title Stewart Smooth
#' @name smoothy
#' @description This function computes a distance weighted mean. It offers the 
#' same parameters as \code{\link{stewart}}: user defined distance matrix, user 
#' defined impedance function (power or exponential), user defined exponent.
#' @param knownpts sp object (SpatialPointsDataFrame or SpatialPolygonsDataFrame);
#' this is the set of known observations to estimate the potentials from.
#' @param unknownpts sp object (SpatialPointsDataFrame or SpatialPolygonsDataFrame); 
#' this is the set of unknown units for which the function computes the estimates. 
#' Not used when \code{resolution} is set up. (optional)
#' @param matdist matrix; distance matrix between known observations and unknown 
#' units for which the function computes the estimates. Row names match the row 
#' names of \code{knownpts} and column names match the row names of 
#' \code{unknownpts}. \code{matdist} can contain any distance metric (time 
#' distance or euclidean distance for example). If \code{matdist} is NULL, the distance 
#' matrix is built with \code{\link{CreateDistMatrix}}. (optional)
#' @param varname character; name of the variable in the \code{knownpts} dataframe 
#' from which potentials are computed. Quantitative variable with no negative values. 
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
#' @param resolution numeric; resolution of the output SpatialPointsDataFrame
#'  (in map units). If resolution is not set, the grid will contain around 7250 
#'  points. (optional)
#' @param mask sp object; the spatial extent of this object is used to 
#' create the regularly spaced SpatialPointsDataFrame output. (optional)
#' @param bypassctrl logical; bypass the distance matrix size control (see 
#' \code{\link{CreateDistMatrix}} Details).
#' @param longlat	logical; if FALSE, Euclidean distance, if TRUE Great Circle 
#' (WGS84 ellipsoid) distance.
#' @return SpatialPointsDataFrame with the computed potentials in a new field 
#' named \code{OUTPUT}
#' @examples 
#' # Create a SpatialPointsDataFrame grid of spatMask extent and 200 meters 
#' # resolution
#' data(spatData)
#' mygrid <- CreateGrid(w = spatMask, resolution = 200)
#' # Create a distance matrix between known points (spatPts) and mygrid
#' mymat <- CreateDistMatrix(knownpts = spatPts, unknownpts = mygrid)
#' # Compute Stewart potentials from known points (spatPts) on a given 
#' # grid (mygrid) using a given distance matrix (mymat)
#' mystewart <- smoothy(knownpts = spatPts, unknownpts = mygrid, 
#'                      matdist = mymat, varname = "Capacite", 
#'                      typefct = "exponential", span = 1250, 
#'                      beta = 3, mask = spatMask)
#' # Compute Stewart potentials from known points (spatPts) on a 
#' # grid defined by its resolution
#' mystewart2 <- smoothy(knownpts = spatPts, varname = "Capacite", 
#'                       typefct = "exponential", span = 1250, beta = 3, 
#'                       resolution = 200, mask = spatMask)
#' # The two methods have the same result
#' identical(mystewart, mystewart2)
#' # the function output a SpatialPointsDataFrame
#' class(mystewart)
#' # Computed values
#' summary(mystewart$OUTPUT)
#' @seealso \link{stewart}.
#' @import sp
#' @import raster
#' @export
smoothy <- function(knownpts,
                    unknownpts = NULL,
                    matdist = NULL,
                    varname,
                    typefct = "exponential", 
                    span,
                    beta,
                    resolution = NULL,
                    mask = NULL,
                    bypassctrl = FALSE, longlat = TRUE)
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
                                  bypassctrl = bypassctrl, longlat = longlat)
    }
  } else {
    unknownpts <- CreateGrid(w = if(is.null(mask)){knownpts} else {mask}, 
                             resolution = resolution) 
    matdist <- CreateDistMatrix(knownpts = knownpts, unknownpts = unknownpts, 
                                bypassctrl = bypassctrl, longlat = longlat) 
  }
  
  matdens <- ComputeInteractDensity(matdist = matdist, typefct = typefct,
                                    beta = beta, span = span)
  
  matopport <- ComputeOpportunity(knownpts = knownpts, matdens = matdens, 
                                  varname = varname)
  
  unknownpts <- ComputeSmooth(unknownpts = unknownpts, matdens = matdens,
                              matopport = matopport)
  unknownpts@data
  return(unknownpts)
}


ComputeSmooth<- function(unknownpts, matopport, matdens)
{
  unknownpts@data$OUTPUT <- apply(matopport, 2, sum, na.rm = TRUE) / colSums(matdens)
  return(unknownpts)
}
