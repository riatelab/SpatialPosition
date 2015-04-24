#' @title Stewart Potentials
#' @name stewart
#' @description This function computes the potentials as defined by J.Q. Stewart (1942).
#' @param knownpts sp object (SpatialPointsDataFrame or SpatialPolygonsDataFrame);
#' this is the set of known observations to estimate the potentials from.
#' @param unknownpts sp object (SpatialPointsDataFrame or SpatialPolygonsDataFrame); 
#' this is the set of unknown units for which the function computes the estimates. 
#' Not used when \code{resolution} is set up. (optional)
#' @param matdist matrix; a distance matrix. Row names match the first 
#' column of the \code{knownpts} object dataframe. Column names match the first column 
#' of the \code{unknownpts} object dataframe. (optional)
#' @param varname character; name of the variable in the \code{knownpts} dataframe from which potentials are computed.
#' Quantitative variable with no negative values. 
#' @param typefct character; spatial interaction function. Options are "pareto" 
#' (default, means power law) or "exponential".
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
#' @param longlat logical; euclidean distance (FALSE, default) or Great Circle distance (TRUE).
#' If TRUE inputs are expected in the WGS84 reference system.
#' @param mask sp object; the spatial extent of this object is used to 
#' create the regularly spaced SpatialPointsDataFrame output. (optional)
#' @details If \code{unknownpts} is NULL then \code{resolution} must be used. 
#' @return SpatialPointsDataFrame with the computed potentials in a new field nammed \code{OUTPUT}
#' @seealso \link{stewart}, \link{rasterStewart}, \link{plotStewart}, \link{contourStewart}, \link{CreateGrid}, \link{CreateDistMatrix}.
#' @examples 
#' # Create a SpatialPointsDataFrame grid of spatMask extent and 200 meters 
#' # resolution
#' data(spatData)
#' mygrid <- CreateGrid(w = spatMask, resolution = 200)
#' # Create a distance matrix between known points (spatPts) and mygrid
#' mymat <- CreateDistMatrix(knownpts = spatPts, unknownpts = mygrid, 
#'                           longlat = FALSE)
#' # Compute Stewart potentials from known points (spatPts) on a given 
#' # grid (mygrid) using a given distance matrix (mymat)
#' mystewart <- stewart(knownpts = spatPts, unknownpts = mygrid, 
#'                      matdist = mymat, varname = "Capacite", 
#'                      typefct = "exponential", span = 1250, 
#'                      beta = 3, longlat = FALSE, mask = spatMask)
#' # Compute Stewart potentials from known points (spatPts) on a 
#' # grid defined by its resolution
#' mystewart2 <- stewart(knownpts = spatPts, varname = "Capacite", 
#'                       typefct = "exponential", span = 1250, beta = 3, 
#'                       resolution = 200, longlat = FALSE, mask = spatMask)
#' # The two methods have the same result
#' identical(mystewart, mystewart2)
#' # the function output a SpatialPointsDataFrame
#' class(mystewart)
#' # Computed values
#' summary(mystewart$OUTPUT)
#' @references 
#' STEWART J.Q. (1942) "Measure of the influence of a population at a distance", Sociometry, 5(1): 63-71.  
#' @import sp
#' @import raster
#' @export
stewart <- function(knownpts,
                    unknownpts = NULL,
                    matdist = NULL,
                    varname,
                    typefct = "exponential", 
                    span,
                    beta,
                    resolution = 2000,
                    longlat = FALSE, 
                    mask = NULL)
{
  TestSp(knownpts)
  if (!is.null(unknownpts)){  
    TestSp(unknownpts)
    if(identicalCRS(knownpts,unknownpts) == FALSE){
      stop(paste("Inputs (",quote(knownpts), " and ",quote(unknownpts),
                 ") do not use the same projection", sep = ""),call. = FALSE)
    }
    if (!is.null(matdist)){
      matdist <- UseDistMatrix(matdist =matdist, knownpts = knownpts, 
                               unknownpts =  unknownpts) 
    }else{
      matdist <- CreateDistMatrix(knownpts = knownpts, unknownpts = unknownpts, 
                                  longlat = longlat) 
    }
  } else {
    unknownpts <- CreateGrid(w = if(is.null(mask)){knownpts} else {mask}, 
                             resolution = resolution) 
    matdist <- CreateDistMatrix(knownpts = knownpts, unknownpts = unknownpts, 
                                longlat = longlat) 
  }
  
  
  matdens <- ComputeInteractDensity(matdist = matdist, typefct = typefct,
                                    beta = beta, span = span)
  
  matopport <- ComputeOpportunity(knownpts = knownpts, matdens = matdens, 
                                  varname = varname)
  
  unknownpts <- ComputePotentials(unknownpts = unknownpts, 
                                  matopport = matopport)
  
  return(unknownpts)
}

#' @title Create a Raster from a Stewart SpatialPointsDataFrame
#' @name rasterStewart
#' @description This function creates a raster from a regularly spaced 
#' Stewart SpatialPointsDataFrame (output of the \code{\link{stewart}} function). 
#' @param x sp object (SpatialPointsDataFrame); output of the \code{stewart} function.
#' @param mask sp object (SpatialPolygonsDataFrame); this object is used to clip the raster. (optional)
#' @return Raster of potential values.
#' @seealso \link{stewart}, \link{rasterStewart}, \link{plotStewart}, \link{contourStewart}, \link{CreateGrid}, \link{CreateDistMatrix}.
#' @examples
#' data(spatData)
#' # Compute Stewart potentials from known points (spatPts) on a
#' # grid defined by its resolution
#' mystewart <- stewart(knownpts = spatPts, varname = "Capacite",
#'                      typefct = "exponential", span = 1000, beta = 3,
#'                      resolution = 50, longlat = FALSE, mask = spatMask)
#' # Create a raster of potentials values
#' mystewartraster <- rasterStewart(x = mystewart, mask = spatMask)
#' plot(mystewartraster)
#' @import sp
#' @import raster
#' @export
rasterStewart <- function(x, mask = NULL){
  gridded(x) <- TRUE
  r <- raster(x)
  rasterx <- rasterize(x, r, field = 'OUTPUT')
  if(!is.null(mask)){
    TestSp(mask)
    rasterx <- mask(rasterx, mask = mask)
  }
  return(rasterx)
}




#' @title Plot a Stewart Raster
#' @name plotStewart
#' @description This function plots the raster produced by the \code{\link{rasterStewart}} function.
#' @param x raster; output of the \code{\link{rasterStewart}} function.
#' @param add logical; if TRUE the raster is added to the current plot, if FALSE the raster is displayed in a new plot.
#' @param breaks numeric; vector of break values to map. If used, 
#' this parameter overrides \code{typec} and \code{nclass} parameters 
#' @param typec character; either "equal" or "quantile", how to discretize the values.
#' @param nclass numeric (integer), number of classes.
#' @param legend.rnd numeric (integer); number of digits used to round the values displayed in the legend.
#' @param col function; color ramp function, such as \code{\link{colorRampPalette}}.
#' @return Display the raster nicely and return the list of break values (invisible).
#' @seealso \link{stewart}, \link{rasterStewart}, \link{plotStewart}, \link{contourStewart}, \link{CreateGrid}, \link{CreateDistMatrix}.
#' @examples 
#' data(spatData)
#' # Compute Stewart potentials from known points (spatPts) on a
#' # grid defined by its resolution
#' mystewart <- stewart(knownpts = spatPts, varname = "Capacite",
#'                      typefct = "exponential", span = 1000, beta = 3,
#'                      resolution = 50, longlat = FALSE, mask = spatMask)
#' # Create a raster of potentials values
#' mystewartraster <- rasterStewart(x = mystewart, mask = spatMask)
#' # Plot stewart potentials nicely
#' plotStewart(x = mystewartraster, add = FALSE, nclass = 5)
#' # Can be used to obtain break values
#' break.values <- plotStewart(x = mystewartraster, add = FALSE, nclass = 5)
#' break.values
#' @import sp
#' @import raster
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
  plot(x, legend.only=TRUE, col = col, 
       breaks=round(bks,legend.rnd))
  breaks <- bks
  return(invisible(breaks))
}


#' @title Create a SpatialPolygonsDataFrame or a SpatialLinesDataFrame from a 
#' Stewart Raster
#' @name contourStewart 
#' @description This function create a SpatialPolygonsDataFrame or a SpatialLinesDataFrame from the Stewart raster.
#' @param x raster; output of the \code{\link{rasterStewart}} function.
#' @param breaks numeric; a vector of break values. 
#' @param type character; "poly" or "line". WARNING: the poly option is experimental and needs the rgeos package.
#' @return The ouput of the function is a SpatialPolygonsDataFrame (\code{type = "poly"}) or a SpatialLinesDataFrame (\code{type = "line"}).
#' @seealso \link{stewart}, \link{rasterStewart}, \link{plotStewart}, \link{contourStewart}, \link{CreateGrid}, \link{CreateDistMatrix}.
#' @examples 
#' data(spatData)
#' # Compute Stewart potentials from known points (spatPts) on a
#' # grid defined by its resolution
#' mystewart <- stewart(knownpts = spatPts, varname = "Capacite",
#'                      typefct = "exponential", span = 1000, beta = 3,
#'                      resolution = 50, longlat = FALSE, 
#'                      mask = spatMask)
#' # Create a raster of potentials values
#' mystewartraster <- rasterStewart(x = mystewart, mask = spatMask)
#' # Display the raster and get break values
#' break.values <- plotStewart(x = mystewartraster)
#' # Create contour SpatialLinesDataFrame
#' mystewartcontourlines <- contourStewart(x = mystewartraster,
#'                                         breaks = break.values,
#'                                         type = "line")
#' mystewartcontourlines@@data
#' plot(mystewartcontourlines, col = "grey20", add=TRUE)
#' plot(spatMask, lwd = 1.25, add=TRUE)
#' # Create contour SpatialPolygonsDataFrame
#' mystewartcontourpoly<- contourStewart(x = mystewartraster,
#'                                       breaks = break.values,
#'                                       type = "poly")
#' unique(mystewartcontourpoly@@data)
#' plot(spatMask, col = "grey80", border = "grey50")
#' plot(mystewartcontourpoly, col = "#0000ff50", border = "grey40", 
#'      add = TRUE)
#' plot(spatPts, cex = 0.5, pch = 20, col  = "#ff000050", add = TRUE)
#' @import sp
#' @import raster
#' @export
contourStewart <- function(x, breaks, type = "line"){
  if (type=="line"){
    return(rasterToContour(x = x, levels = breaks))
  } 
  if (type=="poly"){
    if (!requireNamespace("rgeos", quietly = TRUE)) {
      stop("'rgeos' package needed for this function to work. Please install it.",
           call. = FALSE)
    }
    if(!'package:rgeos' %in% search()){
      attachNamespace('rgeos')
    }
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
    x <- sp::SpatialPolygons(Srl = SPlist)
    x <- sp::SpatialPolygonsDataFrame(Sr = x, 
                                      data = data.frame(levels = SPlevels))
    return (x)
  } else {
    stop(paste("type must be either 'SpatialLinesDataFrame' or", 
               "'SpatialPolygonsDataFrame'", sep =""), call. = FALSE)
  }
}
