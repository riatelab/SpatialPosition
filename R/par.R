#' @title Stewart Potentials Parallel
#' @name stewartParallel
#' @description This function computes the potentials with parallel computing. 
#' EXPERIMENTAL use. Only Great Circle Distances are used. 
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
#' @param cl number of clusters
#' @param chunks chunks size of unknowpts
#' @return SpatialPointsDataFrame with the computed potentials in a new field 
#' named \code{OUTPUT}
#' @seealso \link{rasterStewart}, \link{plotStewart}, \link{quickStewart},
#' \link{rasterToContourPoly}, \link{CreateGrid}, \link{CreateDistMatrix}.
#' @examples 
#' \dontrun{
#' if(require(cartography)){
#' system.time(
#'   s1 <- stewart(knownpts = nuts3.spdf,resolution = 40000, 
#'                varname = "pop2008",
#'                typefct = "exponential", span = 100000,
#'                beta = 3, mask = nuts3.spdf)
#' )
#' 
#' system.time(
#'   s2 <- stewartParallel(knownpts = nuts3.spdf, resolution = 40000, 
#'                         varname = "pop2008",
#'                         typefct = "exponential", span = 100000,
#'                         beta = 3, mask = nuts3.spdf, cl = 4, chunks = 500)
#' )
#' 
#' identical(s1, s2)
#' 
#' r2 <- rasterStewart(s2)
#' c2 <- rasterToContourPoly(r2, breaks = c(0,10000,100000,1000000,2004000,
#'                                          5000000,10040000, 20000000200000000), 
#'                           mask = nuts3.spdf)
#' 
#'  par(mar = c(0,0,0,0))
#'  bks <- sort(unique(c(c2$min, c2$max)))
#'  choroLayer(spdf = c2,
#'             var = "center", legend.pos = "topleft",
#'             breaks = bks, border = NA)
#' }
#' }
#' @references 
#' STEWART J.Q. (1942) "Measure of the influence of a population at a distance", Sociometry, 5(1): 63-71.  
#' @import sp
#' @import raster
#' @import foreach
#' @import doParallel
#' @export
stewartParallel <- function(knownpts, unknownpts = NULL,
                            varname,
                            typefct = "exponential",
                            span, beta, resolution = NULL,
                            mask = NULL, cl = 2, chunks = 100){
  
  
  if (is.null(unknownpts)){
    unknownpts <- CreateGrid(w = if(is.null(mask)){knownpts} else {mask},
                             resolution = resolution)
    unknownpts2 <- unknownpts
  }else{
    if(methods::is(object = unknownpts, "SpatialPolygons")){
      unknownpts2 <- SpatialPointsDataFrame(coordinates(unknownpts),
                                            data = unknownpts@data,
                                            proj4string = unknownpts@proj4string)
    }
  }
  
  if(methods::is(object = knownpts, "SpatialPolygons")){
    knownpts <- SpatialPointsDataFrame(coordinates(knownpts),
                                       data = knownpts@data,
                                       proj4string = knownpts@proj4string)
  }
  
  
  if(sp::is.projected(knownpts)){
    knownpts <- sp::spTransform(knownpts,"+init=epsg:4326")
    unknownpts2 <- sp::spTransform(unknownpts2,"+init=epsg:4326")
  }
  
  if (is.null(cl)){
    cl <- parallel::detectCores(all.tests = FALSE, logical = FALSE)
  }
  cl <- parallel::makeCluster(cl)
  doParallel::registerDoParallel(cl)
  
  
  # sequence pour découper la grille
  sequence <- unique(c(seq(1,nrow(unknownpts2), chunks),nrow(unknownpts2)+1) )
  
  # nombre d'itération
  lseq <- length(sequence)-1
  
  # cut the unknownpts
  ml <- list()
  for  (i in 1:lseq){
    # cat(sequence[i], "-",sequence[i+1]-1, "\n")
    ml[[i]] <- unknownpts2[(sequence[i]):(sequence[i+1]-1),]
  }
  
  ls <- foreach(i = ml, .packages = c('SpatialPosition'),
                .combine = rbind, .inorder = TRUE) %dopar% {
                  mat <- CreateDistMatrix(knownpts = knownpts,
                                          unknownpts = i,
                                          bypassctrl = TRUE)
                  st <- stewart(knownpts = knownpts,
                                unknownpts = i,
                                matdist = mat,
                                typefct = typefct,
                                span = span, beta = beta,
                                varname = varname)
                }
  
  parallel::stopCluster(cl)
  unknownpts$OUTPUT <- ls$OUTPUT
  return(unknownpts)
}
