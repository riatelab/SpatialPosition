#' @title Create a Regularly Spaced Points Grid
#' @name CreateGrid
#' @description This function creates a regular grid of points 
#' from the extent of a given spatial object and a given resolution.
#' @param w sp or sf object; the spatial extent of this object is used to 
#' create the regular grid.
#' @param resolution numeric; resolution of the grid (in map units). If 
#' resolution is not set, the grid will contain around 7500 points. (optional)
#' @param returnclass "sp" or "sf"; class of the returned object.
#' @return The output of the function is a regularly spaced points grid
#'  with the extent of \code{w}.
#' @seealso \link{CreateDistMatrix}
#' @examples
#' # Create a grid of paris extent and 200 meters
#' # resolution
#' library(SpatialPosition)
#' library(sf)
#' data(hospital)
#' mygrid <- CreateGrid(w = paris, resolution = 200, returnclass = "sf")
#' plot(st_geometry(mygrid), cex = 0.1, pch = ".")
#' plot(st_geometry(paris), border="red", lwd = 2, add = TRUE)
#' @importFrom sf st_as_sf st_crs st_bbox
#' @importFrom methods is
#' @export
CreateGrid <- function (w, resolution, returnclass="sp")
{
  # test sf
  if(is(w, "Spatial")){w <- st_as_sf(w)}
  
  boundingBox <- st_bbox(w)
  if(missing(resolution)){
    resolution <- sqrt(((boundingBox[3] - boundingBox[1]) * 
                          (boundingBox[4] - boundingBox[2]))/4000)
  }
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
  spatGrid <- st_as_sf(spatGrid, coords = c("COORDX", "COORDY"),
                       crs = st_crs(w), remove = FALSE)
  if(returnclass=="sp"){spatGrid <- suppressWarnings(as(spatGrid, "Spatial"))}
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
#' If the matrix to compute is too large (more than 100,000,000 cells, more than 
#' 10,000,000 origins or more than 10,000,000 destinations) 
#' the function sends a confirmation message to warn users about the amount of 
#' RAM mobilized. 
#' Use \code{bypassctrl} = TRUE to skip this control.
#' @return A distance matrix, row names are \code{knownpts} row names, column 
#' names are \code{unknownpts} row names.
#' @seealso \link{CreateGrid}
#' @examples
#' # Create a grid of paris extent and 200 meters
#' # resolution
#' data(hospital)
#' mygrid <- CreateGrid(w = paris, resolution = 200, returnclass = "sf")
#' # Create a distance matrix between known hospital and mygrid
#' mymat <- CreateDistMatrix(knownpts = hospital, unknownpts = mygrid, 
#'                           longlat = FALSE)
#' mymat[1:5,1:5]
#' nrow(paris)
#' nrow(mygrid)
#' dim(mymat)
#' @importFrom sf st_centroid st_geometry st_geometry<- st_as_sf st_is_longlat 
#' st_distance st_transform st_is
#' @importFrom methods is
#' @export
CreateDistMatrix  <- function(knownpts, 
                              unknownpts, 
                              bypassctrl = FALSE, 
                              longlat = TRUE)
{
  # test sf
  if(is(knownpts, "Spatial")){knownpts <- st_as_sf(knownpts)}
  if(is(unknownpts, "Spatial")){unknownpts <- st_as_sf(unknownpts)}
  
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
          stop("Computation aborted. Matrix would probably be too large.",
               call. = F)
        }
      } else {
        stop("Computation aborted. Matrix would probably be too large.", 
             call. = F)
      }
    }
}
  # polygon mngmnt
  if(!is(st_geometry(knownpts), "sfc_POINT")){
    st_geometry(knownpts) <- st_centroid(st_geometry(knownpts), 
                                         of_largest_polygon = all(
                                           st_is(knownpts, "MULTIPOLYGON")))
  }
  if(!is(st_geometry(unknownpts), "sfc_POINT")){
    st_geometry(unknownpts) <- st_centroid(st_geometry(unknownpts), 
                                           of_largest_polygon = all(
                                             st_is(unknownpts, "MULTIPOLYGON")))
  }
  
  if(!st_is_longlat(knownpts)){
    if(longlat){
      knownpts <- st_transform(knownpts, 4326)
      unknownpts <- st_transform(unknownpts, 4326)
    }
  }
  x <- st_distance(knownpts, unknownpts)
  mat = as.vector(x)
  dim(mat) = dim(x)
  dimnames(mat) <- list(row.names(knownpts), row.names(unknownpts))
  return(round(mat, digits = 2))
}
