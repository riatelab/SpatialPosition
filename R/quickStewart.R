#' @title Create Polygons of Potentials Contours
#' @name quickStewart
#' @description 
#' This function is a wrapper around \link{stewart}, and \link{isopoly} functions. 
#' Providing only the main parameters of these functions, it simplifies a lot the computation of potentials. 
#' This function creates polygons of potential values. 
#' It also allows to compute directly the ratio between the potentials of two variables. 
#' @param x sp or sf object; this is the set of known observations to 
#' estimate the potentials from.
#' @param spdf a SpatialPolygonsDataFrame.
#' @param df a data frame that contains the values to compute
#' @param spdfid name of the identifier field in spdf, default to the first column 
#' of the spdf data frame. (optional)
#' @param dfid name of the identifier field in df, default to the first column 
#' of df. (optional)
#' @param var name of the numeric field in df used to compute potentials.
#' @param var2 name of the numeric field in df used to compute potentials. 
#' This field is used for ratio computation (see Details).
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
#' @param nclass	numeric; a targeted number of classes (default to 8). Not used 
#' if breaks is set.
#' @param breaks numeric; a vector of values used to discretize the potentials. 
#' @param mask sp or sf object; the spatial extent of this object is used to 
#' create the regularly spaced points output. (optional)
#' @param bypassctrl logical; bypass the distance matrix size control (see 
#' \code{\link{CreateDistMatrix}} Details).
#' @param returnclass "sp" or "sf"; class of the returned object.
#' @return A polyfon object is returned ("sp" or "sf", see \link{isopoly} Value). 
#' @details 
#' If var2 is provided, the ratio between the potentials of var (numerator) 
#' and var2 (denominator) is computed.
#' @seealso \link{stewart}, \link{isopoly}
#' @import sp
#' @import raster
#' @export
#' @examples
#' # load data
#' data("hospital")
#' # Compute potentials
#' pot <- quickStewart(x = hospital,
#'                     var = "capacity",
#'                     span = 1000,
#'                     beta = 2, mask = paris, 
#'                     returnclass = "sf")
#' # cartography
#' if(require("cartography")){
#'   breaks <- sort(c(unique(pot$min), max(pot$max)), decreasing = FALSE)
#'   choroLayer(x = pot,
#'              var = "center", breaks = breaks,
#'              legend.pos = "topleft",
#'              legend.title.txt = "Nb. of Beds")
#' }
#' 
#' # Compute a ratio of potentials
#' hospital$dummy <- hospital$capacity + c(rep(50, 18))
#' pot2 <- quickStewart(x = hospital,
#'                      var = "capacity",
#'                      var2 = "dummy",
#'                      span = 1000,
#'                      beta = 2, 
#'                      mask = paris, 
#'                      returnclass = "sf")
#' # cartography
#' if(require("cartography")){
#'   breaks <- sort(c(unique(pot2$min), max(pot2$max)), decreasing = FALSE)
#'   choroLayer(x = pot2,
#'              var = "center", breaks = breaks,
#'              legend.pos = "topleft",legend.values.rnd = 3,
#'              legend.title.txt = "Nb. of DummyBeds")
#' }
quickStewart <- function(x, spdf, df, spdfid = NULL, dfid = NULL, var, var2, 
                         typefct = "exponential", span, beta, resolution, 
                         mask, nclass = 8, breaks, bypassctrl = FALSE, 
                         returnclass="sp"){
  # IDs  
  if(missing(x)){
    if (is.null(spdfid)){spdfid <- names(spdf@data)[1]}
    if (is.null(dfid)){dfid <- names(df)[1]}
    # Join
    spdf@data <- data.frame(spdf@data[,spdfid], 
                            df[match(spdf@data[,spdfid], df[,dfid]),])
    x <- spdf[!is.na(spdf@data[,dfid]),]
  }
  
  # sp mgmt
  if(is(x, "Spatial")){x <- st_as_sf(x)}
  # pot computation
  pot <- suppressWarnings(stewart(knownpts = x, 
                 varname = var, 
                 typefct = typefct, 
                 span = span, 
                 beta = beta, 
                 resolution = resolution, 
                 mask = mask, 
                 bypassctrl = bypassctrl, 
                 returnclass = "sf"))
  
  if(!missing(var2)){
    if(!is.null(var2)){
      pot2 <- suppressWarnings(stewart(knownpts = x, 
                      varname = var2, 
                      typefct = typefct, 
                      span = span, 
                      beta = beta, 
                      resolution = resolution, 
                      mask = mask, 
                      bypassctrl = bypassctrl, 
                      returnclass = "sf"))
      pot$OUTPUT <- pot$OUTPUT / pot2$OUTPUT
    }
  }
  # dirty stuff to accomodate cartography :-(
  if(!missing(breaks)){
    if(!is.null(breaks)){
      pot <- isopoly(x = pot, breaks = breaks, mask = mask, returnclass = "sf")
    }else{
      pot <- isopoly(x = pot, nclass = nclass, mask = mask, returnclass = "sf")
    }
  }else{
    pot <- isopoly(x = pot, nclass = nclass, mask = mask, returnclass = "sf")
  }
  
  if(returnclass=="sp"){pot <- as(pot, "Spatial")}
  
  return(pot)
  
}
