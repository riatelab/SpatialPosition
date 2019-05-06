#' @title Stewart Smooth
#' @name smoothy
#' @description This function computes a distance weighted mean. It offers the 
#' same parameters as \code{\link{stewart}}: user defined distance matrix, user 
#' defined impedance function (power or exponential), user defined exponent.
#' @param knownpts sp or sf object; this is the set of known observations to 
#' estimate the potentials from.
#' @param unknownpts sp or sf object;
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
#' @param resolution numeric; resolution of the output grid
#'  (in map units). If resolution is not set, the grid will contain around 7250 
#'  points. (optional)
#' @param mask sp or sf object; the spatial extent of this object is used to 
#' create the regularly spaced points output. (optional)
#' @param bypassctrl logical; bypass the distance matrix size control (see 
#' \code{\link{CreateDistMatrix}} Details).
#' @param longlat	logical; if FALSE, Euclidean distance, if TRUE Great Circle 
#' (WGS84 ellipsoid) distance.
#' @param returnclass "sp" or "sf"; class of the returned object.
#' @return Point object with the computed distance weighted mean in a new field 
#' named \code{OUTPUT}. 
#' @examples
#' # Create a grid of paris extent and 200 meters
#' # resolution
#' data(hospital)
#' mygrid <- CreateGrid(w = paris, resolution = 200)
#' # Create a distance matrix between known points (hospital) and mygrid
#' mymat <- CreateDistMatrix(knownpts = hospital, unknownpts = mygrid)
#' # Compute  distance weighted mean from known points (hospital) on a given
#' # grid (mygrid) using a given distance matrix (mymat)
#' mysmoothy <- smoothy(knownpts = hospital, unknownpts = mygrid,
#'                      matdist = mymat, varname = "capacity",
#'                      typefct = "exponential", span = 1250,
#'                      beta = 3, mask = paris, returnclass = "sf")
#' # Compute  distance weighted mean from known points (hospital) on a
#' # grid defined by its resolution
#' mysmoothy2 <- smoothy(knownpts = hospital, varname = "capacity",
#'                       typefct = "exponential", span = 1250, beta = 3,
#'                       resolution = 200, mask = paris, returnclass = "sf")
#' # The two methods have the same result
#' identical(mysmoothy, mysmoothy2)
#' # Computed values
#' summary(mysmoothy$OUTPUT)
#' @seealso \link{stewart}.
#' @import sp
#' @import raster
#' @export
smoothy <- function(knownpts, unknownpts, matdist, varname,
                    typefct = "exponential", span, beta, resolution , mask,
                    bypassctrl = FALSE, longlat = TRUE, returnclass="sp"){
  res <- prepdata(knownpts = knownpts, unknownpts = unknownpts, 
                  matdist = matdist, bypassctrl = bypassctrl, longlat = longlat,
                  mask = mask, resolution = resolution) 
  matdens <- ComputeInteractDensity(matdist = res$matdist, typefct = typefct,
                                    beta = beta, span = span)
  matopport <- ComputeOpportunity(knownpts = res$knownpts, matdens = matdens, 
                                  varname = varname)
  unknownpts <- ComputeSmooth(unknownpts = res$unknownpts, matdens = matdens,
                              matopport = matopport)
  if(returnclass=="sp"){unknownpts <- as(unknownpts, "Spatial")}
  return(unknownpts)
}

