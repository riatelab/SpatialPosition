#' @title Spatial Position Package
#' @name SpatialPosition
#' @description Computes spatial position models: Stewart 
#' potentials, Reilly catchment areas, Huff catchment areas.
#' @docType package
NULL

#' @title SpatialUnits
#' @name spatMask
#' @docType data
NULL

#' @title SpatialUnits
#' @name spatPts
#' @docType data
NULL

#' @title SpatialUnits
#' @name spatUnits
#' @docType data
NULL

#' @title SpatialUnits
#' @name FRdep
#' @docType data
NULL

#' @title SpatialUnits
#' @name FRdep.spdf
#' @docType data
NULL

#' @title SpatialUnits
#' @name spuniv
#' @docType data
NULL

#' @title Create a Regular SpatialPointsDataFrame
#' @name CreateGrid
#' @description This function creates a regular grid of SpatialPointDataFrame 
#' from the extent of a given Spatial*DataFrame and a given resolution.
#' @param w object of class sp (SpatialPointsDataFrame or 
#' SpatialPolygonsDataFrame). The spatial extent of this object is used to 
#' create the regular SpatialPointsDataFrame
#' @param resolution INTEGER, resolution of the output grid (in map units). 
#' @details The output of the function is a SpatialPointsDataFrame of regularly
#' spaced points with the same extent as \code{w}. 
#' @examples 
#' # Calculate potentials
#' data(spatData)
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
#' @description This function creates a distace matrix between two 
#' SpatialPointsDataFrame Objects
#' @param knownpts object of class SpatialPointsDataFrame. 
#' @param unknownpts object of class SpatialPointsDataFrame. 
#' @param longlat if FALSE, Euclidean distance, if TRUE Great Circle distance 
#' @details The function returns a full matrix of distances in the metric of the
#'  points if longlat=FALSE, or in kilometers if longlat=TRUE. This is a wrapper
#'   around the \code{\link{spDists}} function.
#' @examples 
#' # Calculate potentials
#' data(spatData)
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
          (Your computer is likely to collapse under the charge.) [y/n]" )
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
#' @description This function compute the potentials of spatial interaction as 
#' defined by J.Q. Stewart (1950).
#' @param knownpts object of class sp (SpatialPointsDataFrame or 
#' SpatialPolygonsDataFrame). It contains the known points to calculate 
#' potentials from.
#' @param unknownpts object of class sp (SpatialPointsDataFrame or 
#' SpatialPolygonsDataFrame). It contains points to calculate the potentials on. 
#' When used \code{resolution} is not used.
#' @param matdist NUMERIC MATRIX, a distance matrix. Row names match the first 
#' column of the knownpts object dataframe. Column names match the first column 
#' of the knownpts object dataframe. if used \code{longlat} is not used
#' @param varname CHARACTER, name of the variable (in the dataframe of the 
#' Spatial*DataFrame) from which potentials are computed.
#' Quantitative variable with no negative values. 
#' @param typefct CHARACTER, spatial interaction function. Options are "pareto" 
#' (default) or "exponential".
#' If "pareto" the interaction is defined as: (1 + alpha * mDistance) ^ (-beta).
#' If "exponential" the interaction is defined as: 
#' exp(- alpha * mDistance ^ beta).
#' The alpha parameter is computed from parameters given by the user 
#' (beta and reach).
#' @param span NUMERIC, distance where the density of probability of the spatial 
#' interaction function equals 0.5.
#' @param beta NUMERIC, impedance factor for the spatial interaction function.  
#' @param resolution INTEGER, resolution of the output SpatialPointsDataFrame
#'  (in map units). 
#' @param longlat BOOLEAN, type of distance. If TRUE longitudz and latitude are 
#' expected in decimal degrees, in WGS84 reference system.
#' @param mask object of class sp (SpatialPolygonsDataFrame) to clip the map of 
#' potentials.
#' @details The ouput value returned by the function is a SpatialPointsDataFrame 
#' of points with the computed potential in the "POTENTIALS" variable
#' @examples 
#' # Calculate potentials
#' data(spatData)
#' @references Stewart J. Q., Demographic gravitation: evidence and applications
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
#' @description This function compute the potentials of spatial interaction as 
#' defined by J.Q. Stewart (1950).
#' @param knownpts object of class sp (SpatialPointsDataFrame or 
#' SpatialPolygonsDataFrame). It contains the known points to calculate 
#' potentials from.
#' @param unknownpts object of class sp (SpatialPointsDataFrame or 
#' SpatialPolygonsDataFrame). It contains points to calculate the potentials on. 
#' When used \code{resolution} is not used.
#' @param matdist NUMERIC MATRIX, a distance matrix. Row names match the first 
#' column of the knownpts object dataframe. Column names match the first column 
#' of the knownpts object dataframe. if used \code{typedist} is not used
#' @param varname CHARACTER, name of the variable (in the dataframe of the 
#' Spatial*DataFrame) from which potentials are computed.
#' Quantitative variable with no negative values. 
#' @param typefct CHARACTER, spatial interaction function. Options are "pareto" 
#' (default) or "exponential".
#' If "pareto" the interaction is defined as: (1 + alpha * mDistance) ^ (-beta).
#' If "exponential" the interaction is defined as: 
#' exp(- alpha * mDistance ^ beta).
#' The alpha parameter is computed from parameters given by the user 
#' (beta and reach).
#' @param span NUMERIC, distance where the density of probability of the spatial 
#' interaction function equals 0.5.
#' @param beta NUMERIC, impedance factor for the spatial interaction function.  
#' @param resolution INTEGER, resolution of the output grid (in map units). 
#' @param longlat BOOLEAN, type of distance. If TRUE longitudz and latitude are 
#' expected in decimal degrees, in WGS84 reference system.
#' @param mask object of class sp (SpatialPolygonsDataFrame) to clip the map of 
#' potentials.
#' @details The ouput value returned by the function is a SpatialPointsDataFrame 
#' of points with the computed potential in the "POTENTIALS" variable
#' @examples 
#' # Calculate potentials
#' data(spatData)
#' @references Stewart J. Q., Demographic gravitation: evidence and applications
#' , Sociometry, 11(1-2), 31-58.
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
#' @description This function compute the potentials of spatial interaction as 
#' defined by J.Q. Stewart (1950).
#' @param knownpts object of class sp (SpatialPointsDataFrame or 
#' SpatialPolygonsDataFrame). It contains the known points to calculate 
#' potentials from.
#' @param unknownpts object of class sp (SpatialPointsDataFrame or 
#' SpatialPolygonsDataFrame). It contains points to calculate the potentials on. 
#' When used \code{resolution} is not used.
#' @param matdist NUMERIC MATRIX, a distance matrix. Row names match the first 
#' column of the knownpts object dataframe. Column names match the first column 
#' of the knownpts object dataframe. if used \code{typedist} is not used
#' @param varname CHARACTER, name of the variable (in the dataframe of the 
#' Spatial*DataFrame) from which potentials are computed.
#' Quantitative variable with no negative values. 
#' @param typefct CHARACTER, spatial interaction function. Options are "pareto" 
#' (default) or "exponential".
#' If "pareto" the interaction is defined as: (1 + alpha * mDistance) ^ (-beta).
#' If "exponential" the interaction is defined as: 
#' exp(- alpha * mDistance ^ beta).
#' The alpha parameter is computed from parameters given by the user 
#' (beta and reach).
#' @param span NUMERIC, distance where the density of probability of the spatial 
#' interaction function equals 0.5.
#' @param beta NUMERIC, impedance factor for the spatial interaction function.  
#' @param resolution INTEGER, resolution of the output grid (in map units). 
#' @param longlat BOOLEAN, type of distance. If TRUE longitudz and latitude are 
#' expected in decimal degrees, in WGS84 reference system.
#' @param mask object of class sp (SpatialPolygonsDataFrame) to clip the map of 
#' potentials.
#' @details The ouput value returned by the function is a SpatialPointsDataFrame 
#' of points with the computed reilly areas in the "OUTPUT" variable (row.names 
#' of the input Spatial*DataFrame).
#' @examples 
#' # Calculate potentials
#' data(spatData)
#' @references Stewart J. Q., Demographic gravitation: evidence and applications
#' , Sociometry, 11(1-2), 31-58.
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
#' @description This function create a raster layer from gridded 
#' SpatialPointsDataFrame layer outputed by the \code{\link{stewart}}. 
#' @param x object of class sp (SpatialPointsDataFrame) outputed by the
#'  stewart function.
#' @param mask object of class sp (SpatialPolygonsDataFrame) to clip the map of 
#' potentials.
#' @details The ouput of the function is a raster of the Stewart potentials
#' @examples 
#' # Calculate potentials
#' data(spatData)
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
#' @description This function create a raster layer from gridded 
#' SpatialPointsDataFrame layer outputed by the \code{\link{huff}}. 
#' @param x object of class sp (SpatialPointsDataFrame) outputed by the
#'  huff function.
#' @param mask object of class sp (SpatialPolygonsDataFrame) to clip the map of 
#' huff catchement areas.
#' @details The ouput of the function is a raster of the Huff values.
#' @examples 
#' # Calculate potentials
#' data(spatData)
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
#' @description This function create a raster layer from gridded 
#' SpatialPointsDataFrame layer outputed by the \code{\link{reilly}} function. 
#' @param x object of class sp (SpatialPointsDataFrame) outputed by the
#'  huff function.
#' @param mask object of class sp (SpatialPolygonsDataFrame) to clip the map of 
#' reilly catchement areas.
#' @details The ouput of the function is a raster of the Reilly values. 
#' The raster uses a RAT (\code{\link{ratify}}) that contains the 
#' correspondance between raster values and catchement areas values. Use \code{
#' unique(levels(rasterName)[[1]])} to see the correpondance table.
#' @examples 
#' # Calculate potentials
#' data(spatData)
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
#' @description This function plot a potential raster outputed by the 
#' \code{\link{rasterStewart}} function.
#' @param x object of class raster outputed by the 
#' \code{\link{rasterStewart}} function.
#' @param add BOOLEAN, classic add parameter 
#' @param breaks NUMERIC, a vector of break values to map. If used, 
#' this parameter override \code{typec} and \code{nclass} params 
#' @param typec CHARACTER, either "equal" or "quantile", how to discretize the 
#' map.
#' @param nclass NUMERIC, number of class to map.
#' @param legend.rnd NUMERIC, number of digits used to round the figures 
#' displayed in the legend  
#' @param col FUNCTION, colors to display the map like 
#' \code{\link{colorRampPalette}}. 
#' @details The ouput of the function is a plotted raster and the (invisible) 
#' list of break values
#' @examples 
#' # Calculate potentials
#' data(spatData)
#' 
#' @import sp
#' @import raster
#' @export
plotStewart <- function(x, add = TRUE, 
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
#' @description This function plot a potential raster outputed by the 
#' \code{\link{rasterHuff}} function.
#' @param x object of class raster outputed by the 
#' \code{\link{rasterHuff}} function.
#' @param add BOOLEAN, classic add parameter 
#' @param breaks NUMERIC, a vector of break values to map. If used, 
#' this parameter override \code{typec} and \code{nclass} params 
#' @param typec CHARACTER, either "equal" or "quantile", how to discretize the 
#' map.
#' @param nclass NUMERIC, number of class to map.
#' @param legend.rnd NUMERIC, number of digits used to round the figures 
#' displayed in the legend  
#' @param col FUNCTION, colors to display the map like 
#' \code{\link{colorRampPalette}}. 
#' @details The ouput of the function is a plotted raster and the (invisible) 
#' list of break values
#' @examples 
#' # Calculate potentials
#' data(spatData)
#' 
#' @import sp
#' @import raster
#' @export
plotHuff <- function(x, add = TRUE, 
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
#' @description This function plot a potential raster outputed by the 
#' \code{\link{rasterReilly}} function.
#' @param x object of class raster outputed by the 
#' \code{\link{rasterReilly}} function.
#' @param add BOOLEAN, classic add parameter 
#' @param col FUNCTION, colors to display the map, like \code{\link{rainbow}}. 
#' @details The ouput of the function is a plotted raster. 
#' @examples 
#' # Calculate potentials
#' data(spatData)
#' 
#' @import sp
#' @import raster
#' @export
plotReilly <- function(x, add = TRUE, 
                       col =  rainbow){
  nclass <- nrow(unique(levels(x)[[1]]))
  colorReilly <- col(n = nclass)
  plot(x, legend = FALSE, axes = FALSE,
       box = FALSE, col = colorReilly,  add = add)
}

#' @title Create a SpatialPolygonsDataFrame or a SpatialLinesDataFrame from the 
#' Stewart Raster
#' @name contourStewart 
#' @description This function create a SpatialPolygonsDataFrame or a 
#' SpatialLinesDataFrame from the stewart raster
#' @param x object of class raster outputed by the rasterStewart function
#' @param breaks NUMERIC, a vector of break values. 
#' @param type 'poly' or 'line'. WARNING: the poly option is experimental and 
#' can output errors
#' @details The ouput of the function is a SpatialPolygonsDataFrame or a 
#' SpatialLinesDataFrame
#' @examples 
#' # Calculate potentials
#' data(spatData)
#' @import sp
#' @import raster
#' @export
contourStewart <- function(x, breaks, type = "line"){
  if (type=="line"){
    return(rasterToContour(x = x, levels = breaks))
  } 
  if (type=="poly"){
    if (!requireNamespace("rgeos", quietly = T)) {
      stop("'rgeos' package needed for this function to work. Please install it.",
           call. = FALSE)
    }
    loadNamespace("rgeos")
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
        x <- Polygon(coords =  x, hole = F)
        x <- Polygons(srl = list(x), ID = j)
        Plist[j] <- x
      }  
      x <- SpatialPolygons(Srl = Plist)
      x <- rgeos::union(x = x)
      if (class(x) != "SpatialPolygonsDataFrame"){
        x <- SpatialPolygonsDataFrame(Sr = x, 
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
    
    x <- SpatialPolygons(Srl = SPlist)
    x <- SpatialPolygonsDataFrame(Sr = x, data = data.frame(levels = SPlevels))
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

