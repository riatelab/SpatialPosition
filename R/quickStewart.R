#' @title Create a SpatialPolygonsDataFrame of Potentials
#' @name quickStewart
#' @description This function creates a SpatialPolygonsDataFrame of potential values. 
#' It allows also to compute directly the ratio between two variables. 
#' @param spdf a SpatialPolygonsDataFrame.
#' @param df a data frame that contains the values to plot.
#' @param spdfid name of the identifier field in spdf, default to the first column 
#' of the spdf data frame. (optional)
#' @param dfid name of the identifier field in df, default to the first column 
#' of df. (optional)
#' @param var name of the numeric field in df used to compute potentials.
#' @param var2 name of the numeric field in df used to compute potentials, used for ratio (numerator).
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
#' @param breaks numeric; a vector of break values. 
#' @param mask SpatialPolygonsDataFrame; mask used to clip contour shapes.
#' @return The ouput of the function is a SpatialPolygonsDataFrame.
#' @details This function is mainly a wrapper around stewart, rasterStewart and contourStewart.
#' @seealso \link{stewart}, \link{rasterStewart}, \link{plotStewart}, \link{contourStewart}, \link{CreateGrid}, \link{CreateDistMatrix}.
#' @import sp
#' @import raster
#' @export
quickStewart <- function(spdf, df, spdfid = NULL, dfid = NULL, var, 
                         var2 = NULL, 
                         typefct = "exponential", span, 
                         beta, resolution = NULL, 
                         mask = NULL, breaks = NULL){
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
    if(is.null(breaks)){
      breaks <- seq(from = cellStats(ras, min), 
                    to = cellStats(ras, max), length.out = 9)
    }
  }else{
    # missing break
    if(is.null(breaks)){
      breaks <- seq(from = min(pot$OUTPUT, na.rm = TRUE), 
                    to = max(pot$OUTPUT, na.rm = TRUE), 
                    length.out = 9)
      ras <- rasterStewart(pot)
    }
  }

  # Spdf creation
  pot.spdf <- contourStewart(x = ras, 
                             breaks = unique(breaks), 
                             mask = mask, 
                             type = "poly")

  return(pot.spdf)
}
