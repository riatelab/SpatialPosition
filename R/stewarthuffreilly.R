#' @title Spatial Position Package
#' @name SpatialPosition
#' @description Computes spatial position models: Stewart 
#' potentials, Reilly catchment areas, Huff catchment areas.
#' @docType package
NULL

#' @title Spatial Units of Paris 
#' @description 20 spatial units representing the Paris' districts (*arrondissements*), (SpatialPolygonsDataFrame).
#' @name spatUnits
#' @docType data
NULL

#' @title SpatialUnits
#' @description 18 points representing the public hospitals with their capacity (number of beds), (SpatialPointsDataFrame).
#' @name spatPts
#' @docType data
NULL

#' @title Paris Perimeter
#' @description 1 spatial unit representing the Paris' perimeter, (SpatialPolygonsDataFrame). 
#' @name spatMask
#' @docType data
NULL

#' @title Create a Regular SpatialPointsDataFrame
#' @name CreateGrid
#' @description This function creates a regular grid of SpatialPointDataFrame 
#' from the extent of a given Spatial*DataFrame and a given resolution.
#' @param w SP OBJECT, (SpatialPointsDataFrame or 
#' SpatialPolygonsDataFrame). The spatial extent of this object is used to 
#' create the regular SpatialPointsDataFrame
#' @param resolution INTEGER, resolution of the output grid (in map units). 
#' @details The output of the function is a SpatialPointsDataFrame of regularly
#' spaced points with the same extent as \code{w}. 
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
  TestSp(w)
  boundingBox <- bbox(w)
  rounder <- boundingBox %% resolution
  boundingBox[,1] <- boundingBox[,1] - rounder[,1]
  roundermax <- resolution - rounder[,2]
  boundingBox[,2] <- boundingBox[,2] + resolution - rounder[,2]
  boxCoordX <- seq(from = boundingBox[1,1], to = boundingBox[1,2], 
                   by = resolution)
  boxCoordY <- seq(from = boundingBox[2,1], to = boundingBox[2,2], 
                   by = resolution)
  spatGrid <- expand.grid(boxCoordX, boxCoordY)
  idSeq <- seq(1, nrow(spatGrid), 1)
  spatGrid <- data.frame(ID = idSeq, 
                         COORDX = spatGrid[, 1], 
                         COORDY = spatGrid[, 2])
  
  spatGrid <- SpatialPointsDataFrame(coords = spatGrid[ , c(2, 3)], 
                                     data = spatGrid, 
                                     proj4string = CRS(proj4string(w)))
  return(spatGrid)
}


#' @title Create a Distance Matrix Between Two SpatialPointsDataFrame Objects
#' @name CreateDistMatrix
#' @description This function creates a distance matrix between two 
#' SpatialPointsDataFrame Objects
#' @param knownpts SP OBJECT (SpatialPointsDataFrame). 
#' @param unknownpts SP OBJECT (SpatialPointsDataFrame). 
#' @param longlat LOGICAL, Euclidean distance (FALSE, default) or Great Circle distance (TRUE)
#' @details The function returns a full matrix of distances in the metric of the
#'  points if longlat=FALSE, or in kilometers if longlat=TRUE. This is a wrapper
#'   for the \code{\link{spDists}} function.
#' @examples 
#' # Create a SpatialPointsDataFrame grid of spatMask extent and 200 meters 
#' # resolution
#' data(spatData)
#' mygrid <- CreateGrid(w = spatMask, resolution = 200)
#' # Create a distance matrix between known points (spatPts) and mygrid
#' mymat <- CreateDistMatrix(knownpts = spatPts, unknownpts = mygrid, 
#'                           longlat = FALSE)
#' mymat[1:5,1:5]
#' dim(mymat)
#' @import sp
#' @export
CreateDistMatrix  <- function(knownpts, unknownpts, longlat = FALSE)
{
  TestSp(knownpts)
  TestSp(unknownpts)
  if(identicalCRS(knownpts,unknownpts) == FALSE){
    stop(paste("Inputs (",quote(knownpts), " and ",quote(unknownpts),
               ") do not use the same projection", sep = ""),call. = FALSE)
  }
  
  nk <- nrow(knownpts)
  nu <- nrow(unknownpts)
  if(nk * nu > 100000000 | nu > 10000000 | nk > 10000000){
    if (interactive()){
      cat("Do you really want to compute potentials values (from", nk , 
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
        stop("Computation aborted.", call. = F)
      }
    } else {
      stop("Computation aborted. Matrix would probably be too big", call. = F)
    }
  }
  
  matDist <- spDists(x = knownpts, y = unknownpts, longlat = longlat)
  dimnames(matDist) <- list(row.names(knownpts), row.names(unknownpts))
  return(round(matDist, digits = 8))
}


#' @title Stewart Potentials
#' @name stewart
#' @description This function compute the potentials as defined by J.Q. Stewart.
#' @param knownpts SP OBJECT, (SpatialPoints- or SpatialPolygonsDataFrame). 
#' This is the set of known points (observed variable) to estimate the potentials from.
#' @param unknownpts SP OBJECT, (SpatialPoints- or SpatialPolygonsDataFrame). 
#' This is the set of unknown points for which the function computes the estimates. 
#' Not used when \code{resolution} is set up.
#' @param matdist NUMERIC MATRIX, a distance matrix. Row names match the first 
#' column of the knownpts object dataframe. Column names match the first column 
#' of the knownpts object dataframe. 
#' @param varname CHARACTER, name of the variable in the dataframe of the 
#' Spatial*DataFrame from which potentials are computed.
#' Quantitative variable with no negative values. 
#' @param typefct CHARACTER, spatial interaction function. Options are "pareto" 
#' (default, means power law) or "exponential".
#' If "pareto" the interaction is defined as: (1 + alpha * mDistance) ^ (-beta).
#' If "exponential" the interaction is defined as: 
#' exp(- alpha * mDistance ^ beta).
#' The alpha parameter is computed from parameters given by the user 
#' (beta and span).
#' @param span NUMERIC, distance where the density of probability of the spatial 
#' interaction function equals 0.5.
#' @param beta NUMERIC, impedance factor for the spatial interaction function.  
#' @param resolution NUMERIC, resolution of the output SpatialPointsDataFrame
#'  (in map units). 
#' @param longlat LOGICAL, Euclidean distance (FALSE, default) or Great Circle distance (TRUE). 
#' If TRUE longitude and latitude are expected in decimal degrees, in WGS84 reference system.
#' @param mask SP OBJECT (SpatialPolygonsDataFrame), border of the studied region
#' to clip the map of potentials with.
#' @details SpatialPointsDataFrame with the computed potentials in a new field called \code{OUTPUT}
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
#' @references Stewart J. Q. (1948) Demographic gravitation: evidence and applications
#' , Sociometry, 11(1-2), 31-58.
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


#' @title Huff Catchment Areas
#' @name huff
#' @description This function computes the catchment areas as defined by D. Huff (1964).
#' @param knownpts SP OBJECT (SpatialPoints- or SpatialPolygonsDataFrame). 
#' This is the set of known points (observed variable) to draw the catchment areas from.
#' @param unknownpts SP OBJECT, (SpatialPoints- or SpatialPolygonsDataFrame). 
#' This is the set of unknown points for which the function computes the estimates. 
#' Not used when \code{resolution} is set up.
#' @param matdist NUMERIC MATRIX, a distance matrix. Row names match the first 
#' column of the knownpts object dataframe. Column names match the first column 
#' of the knownpts object dataframe. 
#' @param varname CHARACTER, name of the variable in the dataframe of the 
#' Spatial*DataFrame from which values are computed.
#' Quantitative variable with no negative values. 
#' @param typefct CHARACTER, spatial interaction function. Options are "pareto" 
#' (default, means power law) or "exponential".
#' If "pareto" the interaction is defined as: (1 + alpha * mDistance) ^ (-beta).
#' If "exponential" the interaction is defined as: 
#' exp(- alpha * mDistance ^ beta).
#' The alpha parameter is computed from parameters given by the user 
#' (beta and span).
#' @param span NUMERIC, distance where the density of probability of the spatial 
#' interaction function equals 0.5.
#' @param beta NUMERIC, impedance factor for the spatial interaction function.  
#' @param resolution NUMERIC, resolution of the output SpatialPointsDataFrame
#'  (in map units). 
#' @param longlat LOGICAL, Euclidean distance (FALSE, default) or Great Circle distance (TRUE). 
#' If TRUE longitude and latitude are expected in decimal degrees, in WGS84 reference system.
#' @param mask SP OBJECT (SpatialPolygonsDataFrame), border of the studied region
#' to clip the map with.
#' @details SpatialPointsDataFrame with the computed catchment areas in a new field called \code{OUTPUT}
#' @examples 
#' # Create a SpatialPointsDataFrame grid of spatMask extent and 200 meters 
#' # resolution
#' data(spatData)
#' mygrid <- CreateGrid(w = spatMask, resolution = 200)
#' # Create a distance matrix between known points (spatPts) and mygrid
#' mymat <- CreateDistMatrix(knownpts = spatPts, unknownpts = mygrid, 
#'                           longlat = FALSE)
#' # Compute Huff catchment areas from known points (spatPts) on a given 
#' # grid (mygrid) using a given distance matrix (mymat)
#' myhuff <- huff(knownpts = spatPts, unknownpts = mygrid, 
#'                matdist = mymat, varname = "Capacite", 
#'                typefct = "exponential", span = 1250, 
#'                beta = 3, longlat = FALSE, mask = spatMask)
#' # Compute Huff catchment areas from known points (spatPts) on a 
#' # grid defined by its resolution
#' myhuff2 <- huff(knownpts = spatPts, varname = "Capacite", 
#'                       typefct = "exponential", span = 1250, beta = 3, 
#'                       resolution = 200, longlat = FALSE, mask = spatMask)
#' # The two methods have the same result
#' identical(myhuff, myhuff2)
#' # the function output a SpatialPointsDataFrame
#' class(myhuff)
#' @references Huff D. (1964) Defining and Estimating a Trading Area. Journal of Marketing, 28: 34-38.
#' @import sp
#' @import raster
#' @export
huff <- function(knownpts,
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
  
  unknownpts <- ComputeHuff(unknownpts = unknownpts, 
                            matopport = matopport)
  
  return(unknownpts)
}


#' @title Reilly Catchment Areas
#' @name reilly
#' @description This function compute the catchment areas as defined by W.J. Reilly (1931).
#' @param knownpts SP OBJECT (SpatialPoints- or SpatialPolygonsDataFrame). 
#' This is the set of known points (observed variable) to draw the catchment areas from.
#' @param unknownpts SP OBJECT, (SpatialPoints- or SpatialPolygonsDataFrame). 
#' This is the set of unknown points for which the function computes the estimates. 
#' Not used when \code{resolution} is set up.
#' @param matdist NUMERIC MATRIX, a distance matrix. Row names match the first 
#' column of the knownpts object dataframe. Column names match the first column 
#' of the knownpts object dataframe. 
#' @param varname CHARACTER, name of the variable in the dataframe of the 
#' Spatial*DataFrame from which potentials are computed.
#' Quantitative variable with no negative values. 
#' @param typefct CHARACTER, spatial interaction function. Options are "pareto" 
#' (default, means power law) or "exponential".
#' If "pareto" the interaction is defined as: (1 + alpha * mDistance) ^ (-beta).
#' If "exponential" the interaction is defined as: 
#' exp(- alpha * mDistance ^ beta).
#' The alpha parameter is computed from parameters given by the user 
#' (beta and span).
#' @param span NUMERIC, distance where the density of probability of the spatial 
#' interaction function equals 0.5.
#' @param beta NUMERIC, impedance factor for the spatial interaction function.  
#' @param resolution INTEGER, resolution of the output SpatialPointsDataFrame
#'  (in map units). 
#' @param longlat LOGICAL, Euclidean distance (FALSE, default) or Great Circle distance (TRUE). 
#' If TRUE longitude and latitude are expected in decimal degrees, in WGS84 reference system.
#' @param mask SP OBJECT (SpatialPolygonsDataFrame), border of the studied region
#' to clip the map of potentials with.
#' @details SpatialPointsDataFrame with the computed catchment areas in a new field called \code{OUTPUT}. 
#' Values match the row.names of knownpts
#' @examples 
#' # Create a SpatialPointsDataFrame grid of spatMask extent and 200 meters 
#' # resolution
#' data(spatData)
#' mygrid <- CreateGrid(w = spatMask, resolution = 200)
#' # Create a distance matrix between known points (spatPts) and mygrid
#' mymat <- CreateDistMatrix(knownpts = spatPts, unknownpts = mygrid, 
#'                           longlat = FALSE)
#' # Compute Reilly catchment areas from known points (spatPts) on a given 
#' # grid (mygrid) using a given distance matrix (mymat)
#' myreilly2 <- reilly(knownpts = spatPts, unknownpts = mygrid, 
#'                matdist = mymat, varname = "Capacite", 
#'                typefct = "exponential", span = 1250, 
#'                beta = 3, longlat = FALSE, mask = spatMask)
#' row.names(spatPts) <- spatPts$CodHop
#' # Compute Reilly catchment areas from known points (spatPts) on a 
#' # grid defined by its resolution
#' myreilly <- reilly(knownpts = spatPts, varname = "Capacite", 
#'                 typefct = "exponential", span = 1250, beta = 3, 
#'                 resolution = 200, longlat = FALSE, mask = spatMask)
#' # The function output a SpatialPointsDataFrame
#' class(myreilly)
#' # The OUTPUT field values match knownpts row names
#' head(unique(myreilly$OUTPUT))
#' @references Reilly, W. J. (1931) The law of retail gravitation, W. J. Reilly, New York.
#' @import sp
#' @import raster
#' @export
reilly <- function(knownpts,
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
  
  unknownpts <- ComputeReilly(unknownpts = unknownpts, 
                              matopport = matopport)
  
  return(unknownpts)
}





#' @title Create a Raster from a Stewart SpatialPointsDataFrame Layer
#' @name rasterStewart
#' @description This function creates a raster layer from SpatialPointsDataFrame 
#' potential layer, output of the \code{\link{stewart}} function. 
#' @param x SP OBJECT (SpatialPointsDataFrame), output of the \code{stewart} function.
#' @param mask SP OBJECT (SpatialPolygonsDataFrame) to clip the raster with.
#' @details Raster layer with output values (potentials).
#' @examples
#' data(spatData) 
#' # Compute Stewart potentials from known points (spatPts) on a 
#' # grid defined by its resolution
#' mystewart <- stewart(knownpts = spatPts, varname = "Capacite", 
#'                      typefct = "exponential", span = 1250, beta = 3, 
#'                      resolution = 200, longlat = FALSE, mask = spatMask)
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


#' @title Create a Raster from a Huff SpatialPointsDataFrame Layer
#' @name rasterHuff
#' @description This function creates a raster layer from SpatialPointsDataFrame 
#' potential layer, output of the \code{\link{huff}} function. 
#' @param x SP OBJECT (SpatialPointsDataFrame), output of the \code{huff} function.
#' @param mask SP OBJECT (SpatialPolygonsDataFrame) to clip the raster with.
#' @details Raster layer with output values (catchment areas).
#' @examples 
#' data(spatData)
#' # Compute Huff catchment areas from known points (spatPts) on a 
#' # grid defined by its resolution
#' myhuff <- huff(knownpts = spatPts, varname = "Capacite", 
#'                      typefct = "exponential", span = 1250, beta = 3, 
#'                      resolution = 200, longlat = FALSE, mask = spatMask)
#' # Create a raster of huff values
#' myhuffraster <- rasterHuff(x = myhuff, mask = spatMask)
#' plot(myhuffraster)
#'                      
#' @import sp
#' @import raster
#' @export
rasterHuff <- function(x, mask = NULL){
  gridded(x) <- TRUE
  r <- raster(x)
  rasterx <- rasterize(x, r, field = 'OUTPUT')
  if(!is.null(mask)){
    TestSp(mask)
    rasterx <- mask(rasterx, mask = mask)
  }
  return(rasterx)
}

#' @title Create a Raster from a Reilly SpatialPointsDataFrame Layer
#' @name rasterReilly
#' @description This function creates a raster layer from SpatialPointsDataFrame 
#' potential layer, output of the \code{\link{reilly}} function. 
#' @param x, SP OBJECT (SpatialPointsDataFrame), output of the \code{reilly} function.
#' @param mask, SP OBJECT (SpatialPolygonsDataFrame) to clip the raster with.
#' @details Raster layer with output values (catchment areas).
#' The raster uses a RAT (\code{\link{ratify}}) that contains the 
#' correspondance between raster values and catchement areas values. Use \code{
#' unique(levels(rasterName)[[1]])} to see the correpondance table.
#' @examples 
#' data(spatData)
#' row.names(spatPts) <- spatPts$CodHop
#' # Compute Reilly catchment areas from known points (spatPts) on a 
#' # grid defined by its resolution
#' myreilly <- reilly(knownpts = spatPts, varname = "Capacite", 
#'                typefct = "exponential", span = 1250, beta = 3, 
#'                resolution = 200, longlat = FALSE, mask = spatMask)
#' # Create a raster of reilly values
#' myreillyraster <- rasterReilly(x = myreilly, mask = spatMask)
#' plot(myreillyraster)
#' # Correspondance between raster values and reilly areas
#' head(unique(levels(myreillyraster)[[1]]))
#' 
#' @import sp
#' @import raster
#' @export
rasterReilly <- function(x ,mask = NULL){
  gridded(x) <- TRUE
  r <- raster(x)
  x$OUTPUT2 <- as.factor(x$OUTPUT)
  levels(x$OUTPUT2) <- 1:length(levels(x$OUTPUT2) )
  x$OUTPUT2 <- as.numeric(x$OUTPUT2)
  rasterx <- rasterize(x, r, field = 'OUTPUT2')
  if(!is.null(mask)){
    TestSp(mask)
    rasterx <- mask(rasterx, mask = mask)
  }
  ratify(rasterx)
  levels(rasterx) <- data.frame(ID = x$OUTPUT2, idarea = x$OUTPUT)
  return(rasterx)
}



#' @title Plot a Stewart Raster
#' @name plotStewart
#' @description This function plots the raster produced by the \code{\link{rasterStewart}} function.
#' @param x RASTER OBJECT, output of the \code{\link{rasterStewart}} function.
#' @param add LOGICAL, add parameter for the \code{plot} function.
#' @param breaks NUMERIC, vector of break values to map. If used, 
#' this parameter overrides \code{typec} and \code{nclass} parameters 
#' @param typec CHARACTER, either "equal" or "quantile", how to discretize the values.
#' @param nclass INTEGER, number of classes to map.
#' @param legend.rnd INTEGER, number of digits used to round the figures displayed in the legend.
#' @param col FUNCTION, color ramp produced by functions such as \code{\link{colorRampPalette}}.
#' @details Display the raster nicely and return invisible list of break values.
#' @examples 
#' data(spatData)
#' # Compute Stewart potentials from known points (spatPts) on a
#' # grid defined by its resolution
#' mystewart <- stewart(knownpts = spatPts, varname = "Capacite",
#'                      typefct = "exponential", span = 1250, beta = 3,
#'                      resolution = 200, longlat = FALSE, mask = spatMask)
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


#' @title Plot a Huff Raster
#' @name plotHuff
#' @description This function plots the raster produced by the \code{\link{rasterHuff}} function.
#' @param x RASTER OBJECT, output of the \code{\link{rasterHuff}} function.
#' @param add LOGICAL, add parameter for the \code{plot} function.
#' @param breaks NUMERIC, vector of break values to map. If used, 
#' this parameter override \code{typec} and \code{nclass} parameters 
#' @param typec CHARACTER, either "equal" or "quantile", how to discretize the values.
#' @param nclass INTEGER, number of classes to map.
#' @param legend.rnd INTEGER, number of digits used to round the figures displayed in the legend.
#' @param col FUNCTION, color ramp function such as \code{\link{colorRampPalette}}.
#' @details Display the raster nicely and return invisible list of break values.
#' @examples 
#' data(spatData)
#' # Compute Huff catchment areas from known points (spatPts) on a 
#' # grid defined by its resolution
#' myhuff <- huff(knownpts = spatPts, varname = "Capacite", 
#'                      typefct = "exponential", span = 1250, beta = 3, 
#'                      resolution = 200, longlat = FALSE, mask = spatMask)
#' # Create a raster of huff values
#' myhuffraster <- rasterHuff(x = myhuff, mask = spatMask)
#' # Plot Huff values nicely
#' plotHuff(x = myhuffraster)
#' @import sp
#' @import raster
#' @export
plotHuff <- function(x, add = FALSE, 
                     breaks = NULL, typec = "equal", 
                     nclass = 5, legend.rnd = 0, 
                     col =  colorRampPalette(c("#E8F6A4","#445200"))){
  
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



#' @title Plot a Reilly Raster
#' @name plotReilly
#' @description This function plots the raster produced by the \code{\link{rasterReilly}} function.
#' @param x RASTER OBJECT, output of the \code{\link{rasterStewart}} function.
#' @param add LOGICAL, add parameter for the \code{plot} function.
#' @param col FUNCTION, color ramp produced by functions such as \code{\link{colorRampPalette}}.
#' @details Display the raster nicely.
#' @examples 
#' data(spatData)
#' row.names(spatPts) <- spatPts$CodHop
#' # Compute Reilly catchment areas from known points (spatPts) on a 
#' # grid defined by its resolution
#' myreilly <- reilly(knownpts = spatPts, varname = "Capacite", 
#'                typefct = "exponential", span = 1500, beta = 3, 
#'                resolution = 200, longlat = FALSE, mask = spatMask)
#' # Create a raster of reilly values
#' myreillyraster <- rasterReilly(x = myreilly, mask = spatMask)
#' # Plot the raster nicely
#' plotReilly(x = myreillyraster)
#' @import sp
#' @import raster
#' @export
plotReilly <- function(x, add = FALSE, 
                       col =  rainbow){
  nclass <- nrow(unique(levels(x)[[1]]))
  colorReilly <- col(n = nclass)
  plot(x, legend = FALSE, axes = FALSE,
       box = FALSE, col = colorReilly,  add = add)
}


#' @title Create a SpatialPolygonsDataFrame or a SpatialLinesDataFrame from the 
#' Stewart Raster
#' @name contourStewart 
#' @description This function create a SpatialPolygonsDataFrame or a SpatialLinesDataFrame from the Stewart raster
#' @param x RASTER OBJECT, output of the \code{rasterStewart} function.
#' @param breaks NUMERIC, a vector of break values. 
#' @param type CHARACTER, "poly" or "line". WARNING: the poly option is experimental and needs the rgeos package.
#' @details The ouput of the function is a SpatialPolygonsDataFrame or a SpatialLinesDataFrame.
#' @examples 
#' data(spatData)
#' # Compute Stewart potentials from known points (spatPts) on a
#' # grid defined by its resolution
#' mystewart <- stewart(knownpts = spatPts, varname = "Capacite",
#'                      typefct = "exponential", span = 1250, beta = 3,
#'                      resolution = 200, longlat = FALSE, mask = spatMask)
#' # Create a raster of potentials values
#' mystewartraster <- rasterStewart(x = mystewart, mask = spatMask)
#' break.values <- c(0,50,100,150,200,250,300)
#' # Create contour SpatialLinesDataFrame
#' mystewartcontourlines <- contourStewart(x = mystewartraster, 
#'                                        breaks = break.values,
#'                                        type = "line")
#' mystewartcontourlines@@data
#' plotStewart(x = mystewartraster, breaks = break.values)
#' plot(mystewartcontourlines, col = "grey20", add=TRUE)
#' plot(spatMask, lwd = 1.25, add=TRUE)
#' # Create contour SpatialPolygonsDataFrame
#' mystewartcontourpoly<- contourStewart(x = mystewartraster, 
#'                                         breaks = break.values,
#'                                         type = "poly")
#' unique(mystewartcontourpoly@@data)
#' plot(spatMask, col = "grey80", border = "grey50")
#' plot(mystewartcontourpoly, col = "#0000ff50", border = "grey40", add = TRUE)
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


# Internal functions
UseDistMatrix <- function(matdist, knownpts, unknownpts){
  i <- factor(row.names(knownpts), levels = row.names(knownpts))
  j <- factor(row.names(unknownpts), levels = row.names(unknownpts))
  matdist <- matdist[levels(i), levels(j)]
  return(round(matdist, digits = 8))
}

ComputeInteractDensity <- function(matdist, typefct, beta, span)
{
  if(typefct == "pareto") {
    alpha  <- (2 ^ (1 / beta) - 1) / span
    matDens <- (1 + alpha * matdist) ^ (-beta)
  } else if(typefct == "exponential") {
    alpha  <- log(2) / span ^ beta
    matDens <- exp(- alpha * matdist ^ beta)
  } else {
    stop("Please choose a valid interaction function argument (typefct)")
  }
  matDens <- round(matDens, digits = 8)
  return(matDens)
}

ComputeOpportunity <- function(knownpts, matdens, varname = varname)
{
  matOpport <- knownpts@data[, varname] * matdens
  return(round(matOpport, digits = 8))
}

ComputePotentials <- function(unknownpts, matopport)
{
  unknownpts@data$OUTPUT <- apply(matopport, 2, sum, na.rm = TRUE)
  return(unknownpts)
}

ComputeReilly <- function(unknownpts, matopport)
{
  unknownpts@data$OUTPUT <- row.names(matopport)[apply(matopport, 2, which.max)]
  return(unknownpts)
}

ComputeHuff <- function(unknownpts, matopport)
{
  sumCol <- colSums(x = matopport, na.rm = TRUE)
  matOpportPct <- 100 * t(t(matopport) / sumCol)
  matOpportPct[is.na(matOpportPct) | is.infinite(matOpportPct)] <- 0
  unknownpts@data$OUTPUT <- apply(matOpportPct, 2, max, na.rm = TRUE)
  return(unknownpts)
}

TestSp <- function(x){
  if (substr(class(x),1,7) != "Spatial"){
    stop(paste("Your input (",quote(x),") is not a spatial object.", sep=""),
         call. = F)
  }
  if (is.na(x@proj4string)){
    stop(
      paste(
        "Your input (", quote(x),
        ") does not have a valid coordinate reference system.", sep=""),
      call. = F)
  }
}

