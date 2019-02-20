#' @title Create a Regularly Spaced SpatialPointsDataFrame
#' @name CreateGrid
#' @description This function creates a regular grid of SpatialPointsDataFrame 
#' from the extent of a given sp object and a given resolution.
#' @param w sp object; the spatial extent of this object is used to 
#' create the regular SpatialPointsDataFrame.
#' @param resolution numeric; resolution of the grid (in map units). If 
#' resolution is not set, the grid will contain around 7500 points. (optional)
#' @return The output of the function is a SpatialPointsDataFrame of regularly
#' spaced points with the same extent as \code{w}. 
#' @seealso \link{CreateDistMatrix}
#' @examples 
#' # Create a SpatialPointsDataFrame grid of spatMask extent and 200 meters 
#' # resolution
#' data(spatData)
#' mygrid <- CreateGrid(w = spatMask, resolution = 200)
#' plot(mygrid, cex = 0.1, pch = ".")
#' plot(spatMask, border="red", lwd = 2, add = TRUE)
#' @import sp
#' @export
CreateGrid <- function (w, resolution)
{
  # w <- wo
  # resolution <- 5000000
  TestSp(w)
  boundingBox <- bbox(w)
  if(is.null(resolution)){
    resolution <- sqrt(((boundingBox[1,2] - boundingBox[1,1]) * 
                          (boundingBox[2,2] - boundingBox[2,1]))/4000)
  }
  rounder <- boundingBox %% resolution
  boundingBox[,1] <- boundingBox[,1] - rounder[,1]
  boundingBox[,2] <- boundingBox[,2] + resolution - rounder[,2]
  boxCoordX <- seq(from = boundingBox[1,1] - resolution * 2, 
                   to = boundingBox[1,2]+resolution * 2, 
                   by = resolution)
  boxCoordY <- seq(from = boundingBox[2,1] - resolution * 2, 
                   to = boundingBox[2,2] + resolution * 2, 
                   by = resolution)
  
  spatGrid <- expand.grid(boxCoordX, boxCoordY)
  
  idSeq <- seq(1, nrow(spatGrid), 1)
  
  spatGrid <- data.frame(ID = idSeq, 
                         COORDX = spatGrid[, 1], 
                         COORDY = spatGrid[, 2])
  
  spatGrid <- SpatialPointsDataFrame(coords = spatGrid[ , c(2, 3)], 
                                     data = spatGrid, 
                                     proj4string = CRS(proj4string(w)))
  
  result <- tryCatch({
    x <- spTransform(spatGrid, "+init=epsg:4326")
    TRUE
  }, warning = function(war) {
    return(FALSE)
  }, error = function(err) {
    return(FALSE)
  }, finally = {
  })
  
  if (result==FALSE){
    pref <- SpatialPoints(coords = data.frame(x = c(-179.9999,179.9999),
                                              y =  c(89.9999,-89.9999) ),
                          proj4string = CRS("+init=epsg:4326"))
    a <- raster::extent(sp::spTransform(pref, proj4string(w)))
    b <- raster::extent(spatGrid)
    e <- c(b[1], b[2], a[3], a[4])
    spatGrid <- raster::crop(spatGrid, e)
  }
  
  return(spatGrid)
}


#' @title Create a Distance Matrix Between Two Sp Objects
#' @name CreateDistMatrix
#' @description This function creates a distance matrix between two 
#' sp objects (SpatialPointsDataFrame or SpatialPolygonsDataFrame).
#' @param knownpts sp object; rows of the distance matrix.
#' @param unknownpts sp object; columns of the distance matrix.
#' @param bypassctrl logical; bypass the distance matrix size control (see Details).
#' @param longlat	logical; if FALSE, Euclidean distance, if TRUE Great Circle 
#' (WGS84 ellipsoid) distance.
#' @details The function returns a full matrix of distances in meters. 
#' This is a wrapper
#' for the \code{\link{spDists}} function. \cr
#' 
#' If the matrix to compute is too large (more than 100,000,000 cells, more than 
#' 10,000,000 origins or more than 10,000,000 destinations) 
#' the function sends a confirmation message to warn users about the amount of 
#' RAM mobilized. 
#' Use \code{bypassctrl} = TRUE to skip this control.
#' @return A distance matrix, row names are \code{knownpts} row names, column 
#' names are \code{unknownpts} row names.
#' @seealso \link{CreateGrid}
#' @examples 
#' # Create a SpatialPointsDataFrame grid of spatMask extent and 200 meters 
#' # resolution
#' data(spatData)
#' mygrid <- CreateGrid(w = spatMask, resolution = 200)
#' # Create a distance matrix between known spatPts and mygrid
#' mymat <- CreateDistMatrix(knownpts = spatPts, unknownpts = mygrid)
#' mymat[1:5,1:5]
#' nrow(spatPts)
#' nrow(mygrid)
#' dim(mymat)
#' @import sp
#' @import rgdal
#' @export
CreateDistMatrix  <- function(knownpts, 
                              unknownpts, 
                              bypassctrl = FALSE, 
                              longlat = TRUE)
{
  TestSp(knownpts)
  TestSp(unknownpts)
  if(identicalCRS(knownpts,unknownpts) == FALSE){
    stop(paste("Inputs (",quote(knownpts), " and ",quote(unknownpts),
               ") do not use the same projection", sep = ""),call. = FALSE)
  }
  if (bypassctrl == FALSE){
    nk <- nrow(knownpts)
    nu <- nrow(unknownpts)
    if(nk * nu > 100000000 | nu > 10000000 | nk > 10000000){
      if (interactive()){
        cat("Do you really want to this distance matrix (from", nk , 
            "known points to", nu,"estimated values) ? \n 
            (It seems to be a heavy computation.) [y/n]" )
        z <- readLines(con = stdin(), n = 1) 
        while (!z %in% c("n","y")){
          cat ("Enter y or n")
          z <- readLines(con = stdin(), n = 1)  
        }
        if (z == "y"){
          cat("Ok, YOLO!")
        } else {
          stop("Computation aborted. Matrix would probably be too big.",
               call. = F)
        }
      } else {
        stop("Computation aborted. Matrix would probably be too big.", 
             call. = F)
      }
    }
  }
  
  if(methods::is(knownpts, "SpatialPolygons")){
    knownpts <- SpatialPointsDataFrame(coordinates(knownpts), 
                                       data = knownpts@data, 
                                       proj4string = knownpts@proj4string)
    
    
  }
  if(methods::is(unknownpts, "SpatialPolygons")){
    unknownpts <- SpatialPointsDataFrame(coordinates(unknownpts), 
                                         data = unknownpts@data, 
                                         proj4string = unknownpts@proj4string)
  }
  
  if(sp::is.projected(knownpts)){
    if(longlat){
      knownpts <- sp::spTransform(knownpts,"+init=epsg:4326")
      unknownpts <- sp::spTransform(unknownpts,"+init=epsg:4326")
      matDist <- sp::spDists(x = knownpts, y = unknownpts, longlat = TRUE) * 1000
    }else{
      matDist <- sp::spDists(x = knownpts, y = unknownpts, longlat = FALSE)
    }
  }else{
    matDist <- sp::spDists(x = knownpts, y = unknownpts, longlat = TRUE) * 1000
  }
  
  dimnames(matDist) <- list(row.names(knownpts), row.names(unknownpts))
  
  return(round(matDist, digits = 8))
}
