#' @title Create a SpatialPolygonsDataFrame of Potentials Contours
#' @name quickStewart
#' @description 
#' This function is a wrapper around \link{stewart}, \link{rasterStewart} and \link{contourStewart} functions. 
#' Providing only the main parameters of these functions, it simplifies a lot the computation of potentials. 
#' This function creates a SpatialPolygonsDataFrame of potential values. 
#' It also allows to compute directly the ratio between the potentials of two variables. 
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
#'  (in map units). 
#' @param nclass	numeric; a targeted number of classes (default to 8). Not used if breaks is set.
#' @param breaks numeric; a vector of values used to discretize the potentials. 
#' @param mask SpatialPolygonsDataFrame; mask used to clip contours of potentials.
#' @return A SpatialPolygonsDataFrame is returned (see \link{contourStewart} Value). 
#' @details 
#' If var2 is provided the ratio between the potentials of var (numerator) 
#' and var2 (denominator) is computed.
#' @seealso \link{stewart}, \link{rasterStewart}, \link{plotStewart}, 
#' \link{contourStewart}, \link{CreateGrid}, \link{CreateDistMatrix}.
#' @import sp
#' @import raster
#' @export
#' @examples 
#' # load data
#' data("spatData")
#' # Compute a SpatialPolygonsDataFrame of potentials
#' pot.spdf <- quickStewart(spdf = spatPts, 
#'                          df = spatPts@data, 
#'                          var = "Capacite", 
#'                          span = 1000, 
#'                          beta = 2, mask = spatMask)
#' plot(pot.spdf)
#' # cartography
#' if(require("cartography")){
#'   breaks <- c(unique(pot.spdf$min), max(pot.spdf$max))
#'   cartography::choroLayer(spdf = pot.spdf, df = pot.spdf@data,
#'                           var = "mean", breaks = breaks, 
#'                           legend.pos = "topleft",
#'                           legend.title.txt = "Nb. of Beds")
#' }
#' pot.spdf@data
#' 
#' 
#' # Compute a SpatialPolygonsDataFrame of a ratio of potentials
#' spatPts$dummy <- spatPts$Capacite + round(runif(n = nrow(spatPts), 
#'                                                 min = -6, max = 6), 0)
#' pot2.spdf <- quickStewart(spdf = spatPts, 
#'                           df = spatPts@data, 
#'                           var = "Capacite", 
#'                           var2 = "dummy",
#'                           span = 1000, 
#'                           beta = 2, mask = spatMask)
#' # cartography
#' if(require("cartography")){
#'   breaks <- c(unique(pot2.spdf$min), max(pot2.spdf$max))
#'   cartography::choroLayer(spdf = pot2.spdf, df = pot2.spdf@data,
#'                           var = "mean", breaks = breaks, 
#'                           legend.pos = "topleft",legend.values.rnd = 3,
#'                           legend.title.txt = "Nb. of Beds")
#' }
quickStewart <- function(spdf, df, spdfid = NULL, dfid = NULL, var, 
                         var2 = NULL, 
                         typefct = "exponential", span, 
                         beta, resolution = NULL, 
                         mask = NULL, 
                         nclass = 8, breaks = NULL){

  # IDs  
  if (is.null(spdfid)){spdfid <- names(spdf@data)[1]}
  if (is.null(dfid)){dfid <- names(df)[1]}
  
  # Join
  spdf@data <- data.frame(spdf@data[,spdfid], 
                          df[match(spdf@data[,spdfid], df[,dfid]),])
  spdf <- spdf[!is.na(spdf@data[,dfid]),]
  
  # Missing mask 
  if(is.null(mask)){
    bb <- bbox(spdf)
    mm <- matrix(data = c(bb[1,1], bb[2,1], 
                          bb[1,2], bb[2,1], 
                          bb[1,2], bb[2,2], 
                          bb[1,1], bb[2,2], 
                          bb[1,1], bb[2,1]), 
                 nrow = 5, byrow = T)
    mask <- SpatialPolygons(
      Srl = list(Polygons(
        srl = list(Polygon(coords = mm, hole = FALSE)), ID = "id")), 
      proj4string = spdf@proj4string)
  }
  
  # pot computation
  pot <- stewart(knownpts = spdf, 
                 varname = var, 
                 typefct = typefct, 
                 span = span, 
                 beta = beta, 
                 resolution = resolution, 
                 mask = mask)
  
  if(!is.null(var2)){
    pot2 <- stewart(knownpts = spdf, 
                    varname = var2, 
                    typefct = typefct, 
                    span = span, 
                    beta = beta, 
                    resolution = resolution, 
                    mask = mask)
    ras <- rasterStewart(pot) / rasterStewart(pot2)
  }else{
    ras <- rasterStewart(pot)
  }
  # missing break
  if(is.null(breaks)){
    breaks <- seq(from = cellStats(ras, min, na.rm = TRUE), 
                  to = cellStats(ras, max, na.rm = TRUE), 
                  length.out = (nclass+1))
  }
breaks
  # Spdf creation
  pot.spdf <- contourStewart(x = ras, 
                             breaks = unique(breaks), 
                             mask = mask, 
                             type = "poly")
  plot(pot.spdf)
  return(pot.spdf)
}
