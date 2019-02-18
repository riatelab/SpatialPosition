#' @title Create a Regularly Spaced Points Grid
#' @name CreateGrid
#' @description This function creates a regular grid of points 
#' from the extent of a given spatial object and a given resolution.
#' @param w sp or sf object; the spatial extent of this object is used to 
#' create the regular grid.
#' @param resolution numeric; resolution of the grid (in map units). If 
#' resolution is not set, the grid will contain around 7500 points. (optional)
#' @return The output of the function is a regularly spaced points grid with the 
#' extent of \code{w}. If \code{w} is an sp object, the output is a 
#' SpatialPointsDataFrame; if \code{w} is an sf object, the output is an sf 
#' POINT data.frame.
#' @seealso \link{CreateDistMatrix}
#' @examples 
#' # Create a SpatialPointsDataFrame grid of spatMask extent and 200 meters 
#' # resolution
#' data(spatData)
#' mygrid <- CreateGrid(w = spatMask, resolution = 200)
#' plot(mygrid, cex = 0.1, pch = ".")
#' plot(spatMask, border="red", lwd = 2, add = TRUE)
#' @import sp
#' @importFrom methods is
#' @export
CreateGrid <- function (w, resolution)
{
  # test sf
  if(methods::is(w, "Spatial")){w <- sf::st_as_sf(w)}
  
  w <- st_as_sf(spatPts)
  
  boundingBox <- sf::st_bbox(w)
  if(is.null(resolution)){
    resolution <- sqrt(((boundingBox[3] - boundingBox[1]) * 
                          (boundingBox[4] - boundingBox[2]))/4000)
  }
  rounder <- boundingBox %% resolution
  boundingBox[c(1,3)] <- boundingBox[c(1,3)] - rounder[c(1,3)]
  boundingBox[c(2,4)] <- boundingBox[c(2,4)] + resolution - rounder[c(2,4)]
  boxCoordX <- seq(from = boundingBox[1] - resolution * 2, 
                   to = boundingBox[3]+resolution * 2, 
                   by = resolution)
  boxCoordY <- seq(from = boundingBox[2] - resolution * 2, 
                   to = boundingBox[4] + resolution * 2, 
                   by = resolution)
  
  spatGrid <- expand.grid(boxCoordX, boxCoordY)
  spatGrid <- data.frame(ID = 1:nrow(spatGrid),
                         COORDX = spatGrid[, 1], 
                         COORDY = spatGrid[, 2])
  spatGrid <- st_as_sf(spatGrid, coords = c("COORDX", "COORDY"), crs = sf::st_crs(w), remove = FALSE)
  
  # result <- tryCatch({
  #   x <- spTransform(spatGrid, "+init=epsg:4326")
  #   TRUE
  # }, warning = function(war) {
  #   return(FALSE)
  # }, error = function(err) {
  #   return(FALSE)
  # }, finally = {
  # })
  # 
  # if (result==FALSE){
  #   pref <- SpatialPoints(coords = data.frame(x = c(-179.9999,179.9999),
  #                                             y =  c(89.9999,-89.9999) ),
  #                         proj4string = CRS("+init=epsg:4326"))
  #   a <- raster::extent(sp::spTransform(pref, proj4string(w)))
  #   b <- raster::extent(spatGrid)
  #   e <- c(b[1], b[2], a[3], a[4])
  #   spatGrid <- raster::crop(spatGrid, e)
  # }
  # 
  # if(sfsp){
  #   spatGrid <- sf::st_as_sf(spatGrid)
  # }
  
  return(spatGrid)
}


#' @title Create a Distance Matrix Between Two Spatial Objects
#' @name CreateDistMatrix
#' @description This function creates a distance matrix between two 
#' spatial objects (sp or sf objects).
#' @param knownpts sp or sf object; rows of the distance matrix.
#' @param unknownpts sp or sf object; columns of the distance matrix.
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
  # test sf
  if(methods::is(knownpts, "Spatial")){knownpts <- sf::st_as_sf(knownpts)}
  if(methods::is(unknownpts, "Spatial")){unknownpts <- sf::st_as_sf(unknownpts)}
  
  
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

  
  # polygon mngmnt
  if(!methods::is(sf::st_geometry(knownpts), "sfc_POINT")){
    st_geometry(knownpts) <- st_centroid(st_geometry(knownpts))
  }
  if(!methods::is(sf::st_geometry(unknownpts), "sfc_POINT")){
    st_geometry(unknownpts) <- st_centroid(st_geometry(unknownpts))
  }
  
  
  
  if(!st_is_longlat(knownpts)){
    if(longlat){
      knownpts <- sf::st_transform(knownpts, 4326)
      unknownpts <- sf::st_transform(unknownpts, 4326)
    }
  }
  x <- st_distance(knownpts, unknownpts)
  
  y = as.vector(x)
  dim(y) = dim(x)
  
  dimnames(y) <- list(row.names(knownpts), row.names(unknownpts))
  
  return(round(y, digits = 2))
}
