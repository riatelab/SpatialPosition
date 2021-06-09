#' @title Stewart Potentials
#' @name stewart
#' @description This function computes the potentials as defined by J.Q. Stewart (1942).
#' @param knownpts sp or sf object; this is the set of known observations to 
#' estimate the potentials from.
#' @param unknownpts sp or sf object; this is the set of unknown units for which 
#' the function computes the estimates. Not used when \code{resolution} is set 
#' up. (optional)
#' @param matdist matrix; distance matrix between known observations and unknown 
#' units for which the function computes the estimates. Row names match the row 
#' names of \code{knownpts} and column names match the row names of 
#' \code{unknownpts}. \code{matdist} can contain any distance metric (time 
#' distance or euclidean distance for example). If \code{matdist} is missing, the distance 
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
#' @param resolution numeric; resolution of the output grid (in map units). 
#' If resolution is not set, the grid will contain around 7250 points. (optional)
#' @param mask sp or sf object; the spatial extent of this object is used to 
#' create the regularly spaced points output. (optional)
#' @param longlat	logical; if FALSE, Euclidean distance, if TRUE Great Circle 
#' (WGS84 ellipsoid) distance.
#' @param bypassctrl logical; bypass the distance matrix size control (see 
#' \code{\link{CreateDistMatrix}} Details).
#' @param returnclass "sp" or "sf"; class of the returned object.
#' @return Point object with the computed potentials in a new field 
#' named \code{OUTPUT}. 
#' @seealso \link{rasterStewart}, \link{plotStewart}, \link{quickStewart},
#' \link{isopoly}, \link{CreateGrid}, \link{CreateDistMatrix}.
#' @examples
#' # Create a grid of paris extent and 200 meters
#' # resolution
#' data(hospital)
#' mygrid <- CreateGrid(w = paris, resolution = 200, returnclass = "sf")
#' # Create a distance matrix between known points (spatPts) and mygrid
#' mymat <- CreateDistMatrix(knownpts = hospital, unknownpts = mygrid)
#' # Compute Stewart potentials from known points (spatPts) on a given
#' # grid (mygrid) using a given distance matrix (mymat)
#' mystewart <- stewart(knownpts = hospital, unknownpts = mygrid,
#'                      matdist = mymat, varname = "capacity",
#'                      typefct = "exponential", span = 1250,
#'                      beta = 3, mask = paris, returnclass = "sf")
#' # Compute Stewart potentials from known points (spatPts) on a
#' # grid defined by its resolution
#' mystewart2 <- stewart(knownpts = hospital, varname = "capacity",
#'                       typefct = "exponential", span = 1250, beta = 3,
#'                       resolution = 200, mask = paris, returnclass = "sf")
#' # The two methods have the same result
#' identical(mystewart, mystewart2)
#' # the function output a sf data.frame
#' class(mystewart)
#' # Computed values
#' summary(mystewart$OUTPUT)
#' @references 
#' STEWART J.Q. (1942) "Measure of the influence of a population at a distance", Sociometry, 5(1): 63-71.  
#' @importFrom methods is as
#' @importFrom sf st_as_sf
#' @export
stewart <- function(knownpts,unknownpts, matdist, varname, 
                    typefct = "exponential", span, beta, resolution, mask, 
                    bypassctrl = FALSE, longlat = TRUE,  returnclass = "sp"){
  .Deprecated(new = "potential", package = "potential", 
              msg = paste0("stewart(), rasterStewart(), plotStewart(), ",
              "quickStewart(), isopoly(), mcStewart() and smoothy() ", 
              "are deprecated. Use 'potential' package for all ", 
              "potential-related functions."))
  res <- prepdata(knownpts = knownpts, unknownpts = unknownpts, 
                  matdist = matdist, bypassctrl = bypassctrl, longlat = longlat,
                  mask = mask, resolution = resolution) 
  matdens <- ComputeInteractDensity(matdist = res$matdist, typefct = typefct,
                                    beta = beta, span = span)
  matopport <- ComputeOpportunity(knownpts = res$knownpts, matdens = matdens, 
                                  varname = varname)
  unknownpts <- ComputePotentials(unknownpts = res$unknownpts, 
                                  matopport = matopport)
  if(returnclass=="sp"){unknownpts <- suppressWarnings(as(unknownpts, "Spatial"))}
  return(unknownpts)
}


#' @title Create a Raster from a Stewart Regular Grid
#' @name rasterStewart
#' @description This function creates a raster from a regularly spaced 
#' Stewart points grid (output of the \code{\link{stewart}} function). 
#' @param x sp or sf object; output of the \code{stewart} 
#' function.
#' @param mask sp or sf object; this object is used to clip 
#' the raster. (optional)
#' @return Raster of potential values.
#' @seealso \link{stewart}, \link{quickStewart}, \link{plotStewart}, 
#' \link{CreateGrid}, \link{CreateDistMatrix}.
#' @examples
#' library(raster)
#' data(hospital)
#' # Compute Stewart potentials from known points (hospital) on a
#' # grid defined by its resolution
#' mystewart <- stewart(knownpts = hospital, varname = "capacity",
#'                      typefct = "exponential", span = 1000, beta = 3,
#'                      resolution = 100, mask = paris)
#' # Create a raster of potentials values
#' mystewartraster <- rasterStewart(x = mystewart, mask = paris)
#' plot(mystewartraster)
#' @import sp
#' @import raster
#' @export
rasterStewart <- function(x, mask = NULL){
  if(is(x, "sf")){x <- suppressWarnings(as(x, "Spatial"))}
  gridded(x) <- TRUE
  r <- raster(x)
  rasterx <- rasterize(x[!is.na(x$OUTPUT),], r, field = 'OUTPUT')
  if(!is.null(mask)){
    if(is(mask, "sf")){mask <- suppressWarnings(as(mask, "Spatial"))}
    projError(x, mask)
    rasterx <- mask(rasterx, mask = mask)
  }
  return(rasterx)
}




#' @title Plot a Stewart Raster
#' @name plotStewart
#' @description This function plots the raster produced by the 
#' \code{\link{rasterStewart}} function.
#' @param x raster; output of the \code{\link{rasterStewart}} function.
#' @param add logical; if TRUE the raster is added to the current plot, if FALSE 
#' the raster is displayed in a new plot.
#' @param breaks numeric; vector of break values to map. If used, 
#' this parameter overrides \code{typec} and \code{nclass} parameters 
#' @param typec character; either "equal" or "quantile", how to discretize the values.
#' @param nclass numeric (integer), number of classes.
#' @param legend.rnd numeric (integer); number of digits used to round the values 
#' displayed in the legend.
#' @param col function; color ramp function, such as \code{\link{colorRampPalette}}.
#' @return Display the raster nicely and return the list of break values (invisible).
#' @seealso \link{stewart}, \link{rasterStewart}, \link{quickStewart}, 
#' \link{CreateGrid}, \link{CreateDistMatrix}.
#' @examples 
#' data(hospital)
#' # Compute Stewart potentials from known points (hospital) on a
#' # grid defined by its resolution
#' mystewart <- stewart(knownpts = hospital, varname = "capacity",
#'                      typefct = "exponential", span = 1000, beta = 3,
#'                      resolution = 100, mask = paris)
#' # Create a raster of potentials values
#' mystewartraster <- rasterStewart(x = mystewart, mask = paris)
#' # Plot stewart potentials nicely
#' plotStewart(x = mystewartraster, add = FALSE, nclass = 5)
#' # Can be used to obtain break values
#' break.values <- plotStewart(x = mystewartraster, add = FALSE, nclass = 5)
#' break.values
#' @import sp
#' @import raster
#' @importFrom grDevices colorRampPalette
#' @export
plotStewart <- function(x, add = FALSE, 
                        breaks = NULL, typec = "equal", 
                        nclass = 5, legend.rnd = 0, 
                        col =  colorRampPalette(c("#FEA3A3","#980000"))){
  if (!is.null(breaks)){
    bks <- unique(breaks[order(breaks)])
  } else if (typec == "equal"){
    bks <- seq(from = cellStats(x, min), 
               to = cellStats(x, max), length.out = nclass+1)
  } else if (typec == "quantile"){
    bks <- quantile (x, probs = seq(0,1, by = 1/nclass))
  } else {
    stop('Enter a proper discretisation type: "equal" or "quantile"')
  }
  bks <- unique(bks)
  col <- col(length(bks)-1)
  plot(x, breaks = bks, legend = FALSE, axes = FALSE,
       box = FALSE, col = col,  add = add)
  
  nbks <- round(bks,legend.rnd)
  leglab <- rep(NA, (length(nbks)-1))
  for(i in 1:(length(nbks)-1)){
    leglab[i] <- paste("[", nbks[i], " - ", nbks[i+1],"[" ,sep="")
  }
  leglab[i] <- paste( substr(leglab[i],1, nchar(leglab[i])-1), "]", sep="")
  
  graphics::legend(x='topright', legend = rev(leglab), 
                   xpd=T,inset=c(-0.2,0), 
                   fill = rev(col), cex = 0.7, 
                   plot = TRUE, bty = "n", 
                   title = "Potentials")
  
  return(invisible(bks))
}

