#' @title Create a SpatialPolygonsDataFrame or a SpatialLinesDataFrame from a 
#' Stewart Raster
#' @name contourStewart
#' @description This function create a SpatialPolygonsDataFrame or SpatialLinesDataFrame contour from the Stewart raster.
#' @param x raster; output of the \code{\link{rasterStewart}} function. The raster must contain only positive values.
#' @param breaks numeric; a vector of break values. 
#' @param mask SpatialPolygonsDataFrame; mask used to clip contour shapes.
#' @param type character; "poly" or "line". WARNING: the poly option is experimental (see details). It needs the rgeos package.
#' @return The ouput of the function is a SpatialPolygonsDataFrame (\code{type = "poly"}) or a SpatialLinesDataFrame (\code{type = "line"}).
#' @details To obtain a correct SpatialPolygonsDataFrame of potentials follow theses steps: \itemize{
#' \item{Step 1: Create a SpatialPointsDataFrame of potentials with the 
#' stewart function. Do not enter an unknownpts layer, set a resolution, 
#' and set a SpatialPolygonsDataFrame (spmask) as mask.}
#' \item{Step 2: Create a raster from the SpatialPointsDataFrame of potentials 
#' with the rasterStewart function without using a mask.}
#' \item{Step 3: Create the SpatialPolygonsDataFrame of potential with the 
#' contourStewart function and use the spamask SpatialPolygonsDataFrame (Step1) as mask.}
#' }
#' See also the second example in the examples section.
#' @seealso \link{stewart}, \link{rasterStewart}, \link{plotStewart}, \link{contourStewart}, \link{CreateGrid}, \link{CreateDistMatrix}.
#' @import sp
#' @import raster
#' @examples
#' #### Example with type = "line"
#' mystewart <- stewart(knownpts = spatPts, varname = "Capacite",
#'                      typefct = "exponential", span = 1000, beta = 3,
#'                      resolution = 50, longlat = FALSE,
#'                      mask = spatMask)
#' # Create a raster of potentials values
#' mystewartraster <- rasterStewart(x = mystewart, mask = spatMask)
#' # Display the raster and get break values
#' break.values <- plotStewart(x = mystewartraster)
#' # Create contour SpatialLinesDataFrame
#' mystewartcontourpoly <- contourStewart(x = mystewartraster,
#'                                        breaks = break.values,
#'                                        type = "line")
#' # Display the Map
#' plot(spatMask, add=TRUE)
#' plot(mystewartcontourpoly, border = "grey40",add = TRUE)
#' plot(spatPts, cex = 0.8, pch = 20, col  = "black", add = TRUE)
#' 
#' 
#' 
#' #### Example with type = "poly"
#' \dontrun{
#' mystewart <- stewart(knownpts = spatPts, varname = "Capacite",
#'                      typefct = "exponential", span = 1000, beta = 3,
#'                      resolution = 50, longlat = FALSE,
#'                      mask = spatMask)
#' # Create a raster of potentials valuesn, no mask
#' mystewartraster <- rasterStewart(x = mystewart)
#' # Display the raster and get break values
#' break.values <- plotStewart(x = mystewartraster)
#' # Create contour SpatialLinesDataFrame
#' mystewartcontourpoly <- contourStewart(x = mystewartraster,
#'                                        breaks = break.values,
#'                                        mask = spatMask,
#'                                        type = "poly")
#' # Display the map
#' library(cartography)
#' opar <- par(mar = c(0,0,1.1,0))
#' choroLayer(spdf = mystewartcontourpoly, 
#'            df = mystewartcontourpoly@data, 
#'            var = "mean", legend.pos = "topleft",
#'            breaks = break.values, border = "grey90", 
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
contourStewart <- function(x, breaks, mask, type = "line"){
  if (type=="line"){
    return(rasterToContour(x = x, levels = breaks[1:(length(breaks)-1)]))
  } 
  if (type=="poly"){
    if (!requireNamespace("rgeos", quietly = TRUE)) {
      stop("'rgeos' package needed for this function to work. Please install it.",
           call. = FALSE)
    }
    if(!'package:rgeos' %in% search()){
      attachNamespace('rgeos')
    }
    
    
    myproj <- mask@proj4string
    # breaks
    breakss <- breaks
    # raster resolution
    myres <- res(x)[1]
    
    # union the mask
    mask <- rgeos::gUnaryUnion(spgeom = mask, id = NULL)
    
    # buffer around the mask
    maskbuff <- rgeos::gBuffer(mask, byid = FALSE, width = 5 * myres )

    # use a mask around the raster with the maskbuff
    x <- mask(x,maskbuff, updatevalue = -1)
    x[is.na(x)] <- -1
    
    # adapt the breaks to the masked raster
    minx <- min(x[x!=-1])
    maxx <- max(x[x!=-1])
    breaks <- breakss[(breakss>=minx & breakss<maxx)]
    breaks <- c(breakss[which(min(breaks)==breakss)-1],breaks)
    
    
    
    # test breaks
    if(length(breaks)<2){stop("break values do not fit the raster values", 
                              call. = FALSE)}
    # build the contour lines
    cl <- rasterToContour(x, levels = breaks)
    
    cl$level <- as.numeric (as.character(cl$level))
    SPlist <- list()
    SPlevels <- character()
    for (i in cl$level){ 
      linex <- cl[cl@data$level == i,]
      linex <- linex@lines
      linex <- linex[[1]]
      linex <- linex@Lines
      Plist <- NULL
      Plist <- list()
      for (j in 1:length(linex)){
        x <- linex[[j]]@coords
        x <- sp::Polygon(coords =  x, hole = F)
        x <- sp::Polygons(srl = list(x), ID = j)
        Plist[j] <- x
      }  
      x <- sp::SpatialPolygons(Srl = Plist)
      x <- rgeos::union(x = x)
      if (class(x) != "SpatialPolygonsDataFrame"){
        x <- sp::SpatialPolygonsDataFrame(Sr = x, 
                                          data = data.frame(
                                            level = rep(i, length(x))))
      } else {
        x <- x[x@data$count < 2,]
        x@data <- data.frame(level = rep(i, dim(x)[1]))
      }
      SPlist <- c(SPlist , x@polygons  )
      SPlevels <- c(SPlevels,x@data$level)
    }
    for (i in 1:length(SPlist)){
      SPlist[[i]]@ID <- as.character(i)
    }
    x <- sp::SpatialPolygons(Srl = SPlist, proj4string = myproj)
    x <- sp::SpatialPolygonsDataFrame(Sr = x, 
                                      data = data.frame(levels = SPlevels))
    # manage attributes data of the contour spdf
    breaks <- c(breaks, maxx)
    
    x@data <- data.frame(id = paste("id_",row.names(x),sep=""),
                         min = breaks[match(x$levels, breaks)], 
                         max = breaks[match(x$levels, breaks)+1],
                         mean = NA, 
                         stringsAsFactors = FALSE)
    x@data$mean <- (x$min+x$max) / 2 
    row.names(x) <- x$id
    
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
}