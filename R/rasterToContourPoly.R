#' @title Create a SpatialPolygonsDataFrame from a Raster
#' @name rasterToContourPoly
#' @description 
#' This function creates a contour SpatialPolygonsDataFrame from a raster.
#' @param r raster; the raster must contain only positive values.
#' @param nclass numeric; a number of class.
#' @param breaks numeric; a vector of break values. 
#' @param mask SpatialPolygonsDataFrame; mask used to clip contour shapes. 
#' The mask should have a smaller extent than r.
#' @return The ouput of the function is a SpatialPolygonsDataFrame. 
#' The data frame of the outputed SpatialPolygonsDataFrame contains four fields: 
#' id (id of each polygon), min and max (minimum and maximum breaks of the polygon), 
#' center (central values of classes)
#' @details This function uses the rgeos package.
#' @seealso \link{stewart}, \link{rasterStewart}, \link{plotStewart}, 
#' \link{quickStewart}, \link{CreateGrid}, \link{CreateDistMatrix}.
#' @import sp
#' @import raster
#' @examples
#' data("spatData")
#' \dontrun{
#' mystewart <- stewart(knownpts = spatPts, varname = "Capacite",
#'                      typefct = "exponential", span = 1000, beta = 3,
#'                      resolution = 50, mask = spatMask)
#' # Create a raster of potentials values
#' mystewartraster <- rasterStewart(x = mystewart)
#' # Create contour SpatialLinesDataFrame
#' contourpoly <- rasterToContourPoly(r = mystewartraster,
#'                                    nclass = 6,
#'                                    mask = spatMask)
#' # Created breaks
#' bks <- sort(unique(c(contourpoly$min, contourpoly$max)))
#' # Display the map
#' library(cartography)
#' opar <- par(mar = c(0,0,1.2,0))
#' choroLayer(spdf = contourpoly,
#'            df = contourpoly@data,
#'            var = "center", legend.pos = "topleft",
#'            breaks = bks, border = "grey90",
#'            lwd = 0.2,
#'            legend.title.txt = "Potential number\nof beds in the\nneighbourhood",
#'            legend.values.rnd = 0)
#' plot(spatMask, add = TRUE)
#' propSymbolsLayer(spdf = spatPts, df = spatPts@data, var = "Capacite",
#'                  legend.title.txt = "Number of beds",
#'                  col = "#ff000020")
#' layoutLayer(title = "Global Accessibility to Public Hospitals",
#'             south = TRUE, sources = "", author = "")
#' par(opar)
#' }
#' @export
rasterToContourPoly <- function(r, nclass = 8, breaks = NULL, mask = NULL){
  breaks <- get_bks(r = r, nclass =  nclass, breaks = breaks)
  res <- get_mask(r, mask)
  res2 <- adjust_bks(r = res$r, breaks = breaks)
  res3 <- get_poly(r = res2$r, breaks = res2$breaks, finalBreaks = res2$finalBreaks)
  res4 <- add_hole(res3)
  result <- mask_clip(res4, res$mask)
  return(result)
}




# build a mask around the raster
masker <- function(r){
  xy <- sp::coordinates(r)[which(!is.na(raster::values(r))),]
  bb <- sp::bbox(xy)
  xseq <- seq(bb[1,1], bb[1,2], length.out = 50)
  yseq <- seq(bb[2,1], bb[2,2], length.out = 50)
  b <- data.frame(x = c(xseq, rep(bb[1,2], 50), rev(xseq), rep(bb[1,1], 50)), 
                  y = c(rep(bb[2,1], 50), yseq, rep(bb[2,2], 50), rev(yseq)))
  b <- as.matrix(b)
  mask <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(b, hole = FALSE)), 
                                                ID = "1")), 
                              proj4string = sp::CRS(sp::proj4string(r)))
  return(mask)
}


# adjust breaks to fit the raster
get_bks <- function(r, nclass = 8, breaks = NULL){
  # get initial min and max values
  rmin <- raster::cellStats(r, min, na.rm = TRUE)
  rmax <- raster::cellStats(r, max, na.rm = TRUE)
  
  if(is.null(breaks)){
    breaks <- seq(from = rmin, 
                  to = rmax, 
                  length.out = (nclass+1))
  }
  breaks <- sort(unique(c(rmin, breaks[breaks > rmin & breaks < rmax], rmax)))
  return(breaks)
} 

# build a mask around the raster
get_mask <- function(r, mask){
  myres <- raster::res(r)[1]
  if(!is.null(mask)){
    bmask <- TRUE
    # Dissolve the mask
    mask <- rgeos::gUnaryUnion(mask)
    # Is the mask valid
    options(warn=-1)
    if(!rgeos::gIsValid(mask)){
      mask <- rgeos::gBuffer(mask)
    }
    options(warn=0)
    # mask too big or not valid
    if(!rgeos::gIsValid(mask)){
      textx <- "'mask' does not have a valid geometry.\nThe contour SpatialPolygonsDataFrame is built without mask."
      warning(textx, call. = FALSE)
      bmask <- FALSE
    }
  }else{
    bmask <- FALSE
  }
  maskr <- masker(r)
  r <- raster::extend(r, 
                      rgeos::gBuffer(maskr, byid = FALSE, width =  2 * myres ), 
                      value = -1)
  if(!bmask){mask <- maskr}
  return(list(r = r, mask = mask))
}

# adjust breaks to take intoccount 0 and NA values
adjust_bks <- function(r, breaks){
  # adjust breaks if necessary (adapted to mask)
  rmin <- min(r[r!=-1])
  rmax <- max(r[r!=-1])
  breaks <- sort(unique(c(rmin, breaks[breaks > rmin & breaks < rmax], rmax)))
  finalBreaks <- breaks
  
  # deal with 0
  if(breaks[1] <= 0){
    breaks <- breaks + 1
    r <- r + 1
  }
  
  # remove top break
  breaks <- breaks[-(length(breaks))]
  
  # deal with NAs
  r[is.na(r)] <- 0
  
  return(list(breaks = breaks, finalBreaks = finalBreaks, r = r))
}

# get polygons out of the raster 
get_poly <- function(r, breaks, finalBreaks){
  # test breaks
  if(length(breaks)<2){stop("breaks values do not fit the raster values", 
                            call. = FALSE)}
  # build the contour lines polygones
  cl <- rasterToContour(r, levels = breaks)
  cl$level <- as.numeric(as.character(cl$level))
  SPlist <- list()
  SPlevels <- character()
  for (i in cl$level){ 
    linex <- cl[cl@data$level == i,]@lines[[1]]@Lines
    Plist <- list()
    for (j in 1:length(linex)){
      Plist[[j]] <- sp::Polygons(srl = list(sp::Polygon(coords = linex[[j]]@coords, 
                                                        hole = F)), ID = j)
    }  
    x <- rgeos::union(x = sp::SpatialPolygons(Srl = Plist))
    
    if (class(x) != "SpatialPolygonsDataFrame"){
      x <- sp::SpatialPolygonsDataFrame(Sr = x, 
                                        data = data.frame(
                                          level = rep(i, length(x))))
    } else {
      x <- x[x@data$count < 2,]
      x@data <- data.frame(level = rep(i, dim(x)[1]))
    }
    SPlist <- c(SPlist , x@polygons)
    SPlevels <- c(SPlevels, x@data$level)
  }
  
  for (i in 1:length(SPlist)){
    SPlist[[i]]@ID <- as.character(i)
  }
  x <- sp::SpatialPolygonsDataFrame(Sr = sp::SpatialPolygons(Srl = SPlist, 
                                                             proj4string = r@crs), 
                                    data = data.frame(levels = SPlevels))
  bks <- data.frame(b =c(breaks, max(r[r!=-1])), t = finalBreaks)
  # manage attributes data of the contour spdf
  x@data <- data.frame(id = paste("id_",row.names(x),sep=""),
                       min = bks[match(x$levels, bks[,1]),2], 
                       max = bks[match(x$levels, bks[,1])+1,2],
                       center = NA, 
                       stringsAsFactors = FALSE)
  x$center <- (x$min+x$max) / 2 
  row.names(x) <- x$id
  return(x)
}

# add holes into polygons overlaid 
add_hole <- function(final){
  # ring correction
  df <- unique(final@data[,2:4])
  df$id <- 1:nrow(df)
  df <- df[order(df$center, decreasing = T),]
  
  z <- rgeos::gIntersection(final[final$center==df[1,3],],
                            final[final$center==df[1,3],], byid = F,
                            id = as.character(df[1,4]))
  
  for(i in 2:nrow(df)){
    y <- rgeos::gDifference(final[final$center==df[i,3],],
                            final[final$center==df[i-1,3],], byid = F, 
                            id = as.character(df[i,4]))
    if(!is.null(y)){ z <- rbind(z, y)}
  }
  dfx <- data.frame(id = sapply(methods::slot(z, "polygons"), 
                                methods::slot, "ID"))
  row.names(dfx) <- dfx$id
  z <- sp::SpatialPolygonsDataFrame(z, dfx)
  z@data <- df[match(x=z@data$id, table = df$id),c(4,1:3)]
  return(z)
}

# clip polygons with the mask
mask_clip <- function(x, mask){
  # clip the contour spdf with the mask
  final <- rgeos::gIntersection(spgeom1 = x, spgeom2 = mask, byid = TRUE,
                                id = row.names(x))
  df <- data.frame(id = sapply(methods::slot(final, "polygons"),
                               methods::slot, "ID"))
  row.names(df) <- df$id
  final <- sp::SpatialPolygonsDataFrame(Sr = final, data = df)
  final@data <- data.frame(id = final$id, x[match(final$id, x$id),2:4])
  final@plotOrder <- 1:nrow(final)
  return(final)
}
