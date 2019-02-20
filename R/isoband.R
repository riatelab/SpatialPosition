#' @title Create Spatial Polygons Contours from a Raster
#' @name isoStewart
#' @description 
#' This function creates spatial polygons of contours from a raster.
#' @param x sf POINT data.frame; must contain X, Y and OUTPUT fields.
#' @param nclass numeric; a number of class.
#' @param breaks numeric; a vector of break values. 
#' @param mask sf POLYGON data.frame; mask used to 
#' clip contour shapes. 
#' @param xcoords character; name of the X coordinates field in x.
#' @param ycoords character; name of the Y coordinates field in x.
#' @param var character; name of the OUTPUT field in x.
#' @return The output is an sf POLYGON data.frame.
#' The data frame contains four fields: 
#' id (id of each polygon), min and max (minimum and maximum breaks of the polygon), 
#' center (central values of classes).
#' @seealso \link{stewart}, \link{rasterStewart}, \link{plotStewart}, 
#' \link{quickStewart}, \link{CreateGrid}, \link{CreateDistMatrix}.
#' @importFrom sf st_as_sf st_crs st_bbox st_cast st_sf st_sfc st_intersection st_union
#' @importFrom isoband isobands iso_to_sfg
#' @importFrom methods is
#' @export
isoStewart <- function(x, nclass = 8, breaks, mask, 
                       xcoords = "COORDX", ycoords = "COORDY", 
                       var = "OUTPUT"){
  # get initial min and max values
  vmin <- min(x[[var]], na.rm = TRUE)
  vmax <- max(x[[var]], na.rm = TRUE)
  
  if(missing(breaks)){
    breaks <- seq(from = vmin, to = vmax, length.out = (nclass+1))
  }else{
    breaks <- sort(unique(c(vmin, breaks[breaks > vmin & breaks < vmax], vmax)))
  }
  
  m <- matrix(data = x[[var]], nrow = length(unique(x[[xcoords]])), 
              dimnames = list(unique(x[[xcoords]]), unique(x[[ycoords]])))
  
  lev_low = breaks[1:(length(breaks)-1)]
  lev_high = breaks[2:length(breaks)]
  raw <- isobands(x = as.numeric(rownames(m)), 
                  y = as.numeric(colnames(m)), z = t(m), 
                  levels_low = lev_low,
                  levels_high = c(lev_high[-length(lev_high)], vmax + 1e-10))
  
  bands <- iso_to_sfg(raw)
  res <- st_sf(id = 1:length(bands), 
               min = lev_low, 
               max = lev_high,
               geometry = st_sfc(bands), 
               crs = st_crs(x))
  res$center = res$min +(res$max - res$min) / 2
  
  
  if(!missing(mask)){
    if(is(mask, "Spatial")){mask <- st_as_sf(mask)}
    options(warn=-1)
    res <- st_cast(st_intersection(x = res, y = st_union(mask)))
    options(warn=0)
  }
  
  return(res)
}

