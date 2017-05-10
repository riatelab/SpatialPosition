#' @title Stewart Potentials Parallel
#' @name mcStewart
#' @description This function computes Stewart potentials using parallel 
#' computation. 
#' @param knownpts sp object (SpatialPointsDataFrame or SpatialPolygonsDataFrame);
#' this is the set of known observations to estimate the potentials from.
#' @param unknownpts sp object (SpatialPointsDataFrame or SpatialPolygonsDataFrame); 
#' this is the set of unknown units for which the function computes the estimates. 
#' Not used when \code{resolution} is set up. (optional)
#' @param varname character; name of the variable in the \code{knownpts} dataframe 
#' from which potentials are computed. Quantitative variable with no negative values. 
#' @param typefct character; spatial interaction function. Options are "pareto" 
#' (means power law) or "exponential".
#' If "pareto" the interaction is defined as: (1 + alpha * mDistance) ^ (-beta).
#' If "exponential" the interaction is defined as: 
#' exp(- alpha * mDistance ^ beta).
#' The alpha parameter is computed from parameters given by the user 
#' (\code{beta} and \code{span}).
#' @param span numeric; distance where the density of probability of the spatial 
#' interaction function equals 0.5.
#' @param beta numeric; impedance factor for the spatial interaction function.  
#' @param resolution numeric; resolution of the output SpatialPointsDataFrame
#'  (in map units). If resolution is not set, the grid will contain around 7250 
#'  points. (optional)
#' @param mask sp object; the spatial extent of this object is used to 
#' create the regularly spaced SpatialPointsDataFrame output. (optional)
#' @param cl numeric; number of clusters. By default cl is determined using 
#' \code{parallel::detectCores()}.
#' @param size numeric; mcStewart splits unknownpts in chunks, size indicates 
#' the size of each chunks.
#' @param longlat	logical; if FALSE, Euclidean distance, if TRUE Great Circle 
#' (WGS84 ellipsoid) distance.
#' @details The parallel implementation splits potentials computations along 
#' chunks of unknownpts (or chunks of the grid defined using resolution). It only 
#' uses Great Cercle distances (with \code{\link{CreateDistMatrix}}). 
#' @return A Spatial*DataFrame with the computed potentials in a new field 
#' named \code{OUTPUT} is returned. 
#' @seealso \link{stewart}.
#' @examples 
#' \dontrun{
#' if(require(cartography)){
#'   nuts3.spdf@data <- nuts3.df
#'   t1 <- system.time(
#'     s1 <- stewart(knownpts = nuts3.spdf,resolution = 40000, 
#'                   varname = "pop2008",
#'                   typefct = "exponential", span = 100000,
#'                   beta = 3, mask = nuts3.spdf)
#'   )
#'   
#'   t2 <- system.time(
#'     s2 <- mcStewart(knownpts = nuts3.spdf, resolution = 40000, 
#'                      varname = "pop2008",
#'                      typefct = "exponential", span = 100000,
#'                      beta = 3, mask = nuts3.spdf, cl = 3, size = 500)
#'   )
#'   identical(s1, s2)
#'   cat("Elapsed time\n", "stewart:", t1[3], "\n parStewart:",t2[3])
#'   
#'   r2 <- rasterStewart(s2)
#'   c2 <- rasterToContourPoly(r = r2, breaks = c(0,1000000,2000000, 5000000,
#'                                                10000000, 20000000, 200004342), 
#'                             mask = nuts3.spdf)
#'   # cartography
#'   opar <- par(mar = c(0,0,1.2,0))
#'   bks <- sort(unique(c(c2$min, c2$max)))
#'   choroLayer(spdf = c2, var = "center", breaks = bks, border = NA, 
#'              legend.title.txt = "pop")
#'   layoutLayer("potential population", "","", scale = NULL)
#'   par(opar)
#' }
#' }
#' @import sp
#' @import raster
#' @export
mcStewart <- function(knownpts, unknownpts = NULL,
                      varname,
                      typefct = "exponential",
                      span, beta, resolution = NULL,
                      mask = NULL, cl = NULL, size = 1000, 
                      longlat = TRUE){
  
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("'foreach' package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("foreach", quietly = TRUE)) {
    stop("'foreach' package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("'foreach' package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  if (is.null(unknownpts)){
    unknownpts <- CreateGrid(w = if(is.null(mask)){knownpts} else {mask},
                             resolution = resolution)
    unknownpts2 <- unknownpts
  }else{
    # SpatialPolygons to SpatialPoints
    if(methods::is(object = unknownpts, "SpatialPolygons")){
      unknownpts2 <- SpatialPointsDataFrame(coordinates(unknownpts),
                                            data = unknownpts@data,
                                            proj4string = unknownpts@proj4string)
    }
  }
  # SpatialPolygons to SpatialPoints
  if(methods::is(object = knownpts, "SpatialPolygons")){
    knownpts <- SpatialPointsDataFrame(coordinates(knownpts),
                                       data = knownpts@data,
                                       proj4string = knownpts@proj4string)
  }
  
  # force wgs84
  # if(sp::is.projected(knownpts)){
  #   knownpts <- sp::spTransform(knownpts,"+init=epsg:4326")
  #   unknownpts2 <- sp::spTransform(unknownpts2,"+init=epsg:4326")
  # }
  
  # launch multiple cores
  if (is.null(cl)){
    cl <- parallel::detectCores(all.tests = FALSE, logical = FALSE)
  }
  cl <- parallel::makeCluster(cl)
  doParallel::registerDoParallel(cl)
  
  # sequence to split unknowpts
  sequence <- unique(c(seq(1,nrow(unknownpts2), size),nrow(unknownpts2)+1))
  lseq <- length(sequence)-1
  
  # split unknownpts and put it on a list
  ml <- list()
  for  (i in 1:lseq){
    ml[[i]] <- unknownpts2[(sequence[i]):(sequence[i+1]-1),]
  }
  
  ls <- foreach::`%dopar%`(foreach::foreach(i = ml, 
                                            .packages = c('SpatialPosition'),
                                            .combine = rbind, .inorder = TRUE), 
                           {
                             mat <- CreateDistMatrix(knownpts = knownpts,
                                                     unknownpts = i,
                                                     bypassctrl = TRUE, 
                                                     longlat = longlat)
                             st <- stewart(knownpts = knownpts,
                                           unknownpts = i,
                                           matdist = mat,
                                           typefct = typefct,
                                           span = span, beta = beta,
                                           varname = varname)
                           })
  parallel::stopCluster(cl)
  unknownpts$OUTPUT <- ls$OUTPUT
  return(unknownpts)
}



