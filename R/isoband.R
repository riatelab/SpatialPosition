#' @title Create Spatial Polygons Contours from a Raster
#' @name isobander
#' @description 
#' This function creates spatial polygons of contours from a raster.
#' @param x sf POINT data.frame; must contain X, Y and OUTPUT fields.
#' @param nclass numeric; a number of class.
#' @param breaks numeric; a vector of break values. 
#' @param mask sf POLYGON data.frame; mask used to 
#' clip contour shapes. 
#' @return The output is an sf POLYGON data.frame.
#' The data frame contains four fields: 
#' id (id of each polygon), min and max (minimum and maximum breaks of the polygon), 
#' center (central values of classes).
#' @seealso \link{stewart}, \link{rasterStewart}, \link{plotStewart}, 
#' \link{quickStewart}, \link{CreateGrid}, \link{CreateDistMatrix}.
#' @export
isobander <- function(x, nclass = 8, breaks = NULL, mask = NULL){
  # get initial min and max values
  vmin <- min(x$OUTPUT, na.rm = TRUE)
  vmax <- max(x$OUTPUT, na.rm = TRUE)
  
  if(is.null(breaks)){
    breaks <- seq(from = vmin, to = vmax, length.out = (nclass+1))
  }else{
    breaks <- sort(unique(c(vmin, breaks[breaks > vmin & breaks < vmax], vmax)))
  }
  
  m <- matrix(data = x$OUTPUT, nrow = length(unique(x$COORDX)), 
              dimnames = list((unique(x$COORDX)), unique(x$COORDY)))
  
  lev_low = breaks[1:(length(breaks)-1)]
  lev_high = breaks[2:length(breaks)]
  raw <- isoband::isobands(x = as.numeric(rownames(m)), 
                       y = as.numeric(colnames(m)), z = t(m), 
                       levels_low = lev_low,
                       levels_high = c(lev_high[-length(lev_high)], vmax + 1e-10))
  
  bands <- isoband::iso_to_sfg(raw)
  res <- sf::st_sf(id = 1:length(bands), 
                      min = lev_low, 
                      max = lev_high,
                      geometry = sf::st_sfc(bands), 
                      crs = sf::st_crs(x))
  res$center = res$min +(res$max - res$min) / 2
  
  
  if(!missing(mask)){
      res <- sf::st_intersection(x = res, y = mask)
  }
  
  return(res)
}

