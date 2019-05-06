#' @title Create Spatial Polygons Contours from a Raster
#' @name isopoly
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
#' @param returnclass "sp" or "sf"; class of the returned object.
#' @return The output is an sf POLYGON data.frame.
#' The data frame contains four fields: 
#' id (id of each polygon), min and max (minimum and maximum breaks of the polygon), 
#' center (central values of classes).
#' @seealso \link{stewart}.
#' @importFrom sf st_as_sf st_crs st_bbox st_cast st_sf st_sfc st_intersection 
#' st_union st_agr<- st_collection_extract
#' @importFrom isoband isobands iso_to_sfg
#' @importFrom methods is
#' @importFrom lwgeom st_make_valid
#' @examples
#' data(hospital)
#' # Compute Stewart potentials
#' mystewart <- stewart(knownpts = hospital, varname = "capacity",
#'                      typefct = "exponential", span = 1000, beta = 3,
#'                      mask = paris, returnclass = "sf")
#' # Create contour
#' contourpoly <- isopoly(x = mystewart,
#'                        nclass = 6,
#'                        mask = paris, returnclass = "sf")
#' library(sf)
#' plot(st_geometry(contourpoly))
#' if(require(cartography)){
#'   # Created breaks
#'   bks <- sort(unique(c(contourpoly$min, contourpoly$max)))
#'   opar <- par(mar = c(0,0,1.2,0))
#'   # Display the map
#'   choroLayer(x = contourpoly,
#'              var = "center", legend.pos = "topleft",
#'              breaks = bks, border = "grey90",
#'              lwd = 0.2,
#'              legend.title.txt = "Potential number\nof beds in the\nneighbourhood",
#'              legend.values.rnd = 0)
#'   plot(st_geometry(paris), add = TRUE)
#'   propSymbolsLayer(x = hospital, var = "capacity",
#'                    legend.pos = "right",
#'                    legend.title.txt = "Number of beds",
#'                    col = "#ff000020")
#'   layoutLayer(title = "Global Accessibility to Public Hospitals",
#'               sources = "", author = "")
#'   par(opar)
#' }
#' @export
isopoly <- function(x, nclass = 8, breaks, mask, 
                    xcoords = "COORDX", ycoords = "COORDY", 
                    var = "OUTPUT", returnclass="sp"){
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
  iso <- st_sf(id = 1:length(bands), 
               min = lev_low, 
               max = lev_high,
               geometry = st_sfc(bands), 
               crs = st_crs(x))
  iso$center = iso$min + (iso$max - iso$min) / 2
  
  
  st_geometry(iso) <- st_make_valid(st_geometry(iso))
  if(methods::is(st_geometry(iso),"sfc_GEOMETRY")){
    st_geometry(iso) <-   st_collection_extract(st_geometry(iso), "POLYGON")
  }
  
  if(!missing(mask)){
    if(is(mask, "Spatial")){mask <- st_as_sf(mask)}
    st_agr(iso) <- "constant"
    iso <- st_cast(st_intersection(x = iso, y = st_union(mask)))
  }
  if(returnclass=="sp"){iso <- as(iso, "Spatial")}
  return(iso)
}

