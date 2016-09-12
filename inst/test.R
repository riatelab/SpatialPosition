library(SpatialPosition)
library(foreach)
library(doParallel)
load("../LargeSpatialPostition/com.RData")
comS <- com

stewartParallel <- function(knownpts, unknownpts = NULL,
                            varname,
                            typefct = "exponential",
                            span, beta, resolution = NULL,
                            mask = NULL,cl = 2, chunks = 100){


  if (is.null(unknownpts)){
    unknownpts <- CreateGrid(w = if(is.null(mask)){knownpts} else {mask},
                              resolution = resolution)
    unknownpts2 <- unknownpts
  }else{
    if(is(object = unknownpts, "SpatialPolygons")){
      unknownpts2 <- SpatialPointsDataFrame(coordinates(unknownpts),
                                            data = unknownpts@data,
                                            proj4string = unknownpts@proj4string)
    }
  }

  if(is(object = knownpts, "SpatialPolygons")){
    knownpts <- SpatialPointsDataFrame(coordinates(knownpts),
                                       data = knownpts@data,
                                       proj4string = knownpts@proj4string)
  }


  if(sp::is.projected(knownpts)){
    knownpts <- sp::spTransform(knownpts,"+init=epsg:4326")
    unknownpts2 <- sp::spTransform(unknownpts2,"+init=epsg:4326")
  }

  if (is.null(cl)){
    cl <- detectCores(all.tests = FALSE, logical = FALSE)
  }
  cl <- makeCluster(cl)
  registerDoParallel(cl)


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

  stopCluster(cl)
  unknownpts$OUTPUT <- ls$OUTPUT
  return(unknownpts)
}


nuts3.spdf@data <- nuts3.df

system.time(
mystewart <- stewart(knownpts = nuts3.spdf, unknownpts = nuts3.spdf,
                     varname = "pop2008",
                     typefct = "exponential", span = 100000,
                     beta = 3, mask = nuts3.spdf)
)
plot(mystewart)
par(mar = c(0,0,0,0))
choroLayer(mystewart, var = "OUTPUT", border = NA, nclass=8)
stewart()

system.time(

s2 <- stewartParallel(knownpts = nuts3.spdf, unknownpts = nuts3.spdf,
                varname = "pop2008",
                typefct = "exponential", span = 100000,
                beta = 3, mask = nuts3.spdf, cl = 3, chunks = 100)

)

choroLayer(s2, var = "OUTPUT", border = NA, nclass=8)
identical(s2, mystewart)




system.time(
  s <- stewart(knownpts = nuts3.spdf,resolution = 10000,
                       varname = "pop2008",
                       typefct = "exponential", span = 100000,
                       beta = 3, mask = nuts3.spdf)
)
r1 <- rasterStewart(s)
c1 <- rasterToContourPoly(r1, breaks = c(0,10000,100000,1000000,2004000,5000000,10040000, 20000000200000000), mask = nuts3.spdf)
par(mar = c(0,0,0,0))
bks <- sort(unique(c(c1$min, c1$max)))
choroLayer(spdf = c1,
           var = "center", legend.pos = "topleft",
           breaks = bks, border = NA)



system.time(
  s2 <- stewartParallel(knownpts = nuts3.spdf, resolution = 10000,
                        varname = "pop2008",
                        typefct = "exponential", span = 100000,
                        beta = 3, mask = nuts3.spdf, cl = 4, chunks = 50)
)
r2 <- rasterStewart(s2)
c2 <- rasterToContourPoly(r2, breaks = c(0,10000,100000,1000000,2004000,
                                         5000000,10040000, 20000000200000000),
                          mask = nuts3.spdf)
par(mar = c(0,0,0,0))
bks <- sort(unique(c(c2$min, c2$max)))
choroLayer(spdf = c2,
           var = "center", legend.pos = "topleft",
           breaks = bks, border = NA)




summary(s2)
summary(mystewart)
s2@data==mystewart@data

# Les calculs de matrices de distances sont (beaucoup) plus rapides quand les
# objets spatiaux sont en WGS84 :
# Stocker la projection initiale
old.proj <- proj4string(com)
# passage en WGS84
com <- SpatialPointsDataFrame(coordinates(com), data = com@data, proj4string = com@proj4string)
com <- spTransform(com, "+init=epsg:4326")




nbcl <- detectCores(all.tests = FALSE, logical = FALSE)

#setup parallel backend to use 8 processors
cl <- makeCluster(nbcl)
registerDoParallel(cl)


# sequence pour découper la grille
sequence <- unique(c(seq(1,nrow(grid), 100),nrow(grid)+1) )

# nombre d'itération
lseq <- length(sequence)-1

# cut the unknownpts
ml <- list()
for  (i in 1:lseq){
  # cat(sequence[i], "-",sequence[i+1]-1, "\n")
  ml[[i]] <- grid[(sequence[i]):(sequence[i+1]-1),]
}

ls <- foreach(i = ml, .packages = c('SpatialPosition'),
              .combine = rbind, .inorder = TRUE) %dopar% {
                mat <- CreateDistMatrix(knownpts = com,
                                        unknownpts = i,
                                        bypassctrl = TRUE)
                st <- stewart(knownpts = com,
                              unknownpts = i,
                              matdist = mat,
                              typefct = "exponential",
                              span = 50000, beta = 3,
                              var = "POPULATION")
              }

stopCluster(cl)


strt1 <- Sys.time()
print(strt1 - strt0)
#
grid$OUTPUT <- ls$OUTPUT
rasSt <- rasterStewart(grid)
#
a <- rasterToContourPoly(r = rasSt,
                         breaks=c(10,20,50,100,200,500,1000,2000,5000,17588),
                         mask = comS)
#
#




#
# Cartographie
bks <- sort(unique(c(a$min, a$max)))
# Display the map
library(cartography)
dev.off()
opar <- par(mar = c(0,0,1.2,0))
choroLayer(spdf = a,
           var = "center", legend.pos = "topleft",
           breaks = bks, border = NA,
           legend.title.txt = "Population potentielle\ndans un voisinage de 20 km\n(en milliers)",
           legend.values.rnd = 0)
layoutLayer(title = "Pop",
            south = TRUE,
            sources = "", author = "")
par(opar)

#
#
#
#
#
#
#
#
#
#
#
# i3 <- iter(data.frame(x=1:3, y=10), by='row')
# nextElem(i3)
# nextElem(i3)
# nextElem(i3)




devtools::use_build_ignore('inst/test.R', escape = TRUE, pkg = ".")

