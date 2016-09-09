# library(SpatialPosition)
# 
# load("../LargeSpatialPostition/com.RData")
# coco <- com
# 
# # création d'une grille fine (5km)
# grid <- CreateGrid(w = com, resolution = 5000)
# 
# # Les calculs de matrices de distances sont (beaucoup) plus rapides quand les
# # objets spatiaux sont en WGS84 :
# # Stocker la projection initiale
# old.proj <- proj4string(com)
# # passage en WGS84
# com <- SpatialPointsDataFrame(coordinates(com), data = com@data, proj4string = com@proj4string)
# com <- spTransform(com, "+init=epsg:4326")
# grid <- spTransform(grid, "+init=epsg:4326")
# 
# knownpts = com
# unknownpts = grid
# varname = "POPULATION"
# typefct = "exponential"
# span = 20000
# beta = 2
# 
# library(foreach)
# library(doParallel)
# 
# 
# 
# #setup parallel backend to use 8 processors
# cl <- makeCluster(3)
# registerDoParallel(cl)
# 
# # sequence pour découper la grille
# sequence <- unique(c(seq(1,nrow(grid), 10),nrow(grid)+1) )
# 
# # nombre d'itération
# lseq <- length(sequence)-1
# 
# i <- lseq
# 
# ml <- list()
# for  (i in 1:lseq){
#   ml[[i]] <- unknownpts[sequence[i]:sequence[i+1]-1,]
# }
# 
# 
# 
# strt0 <- Sys.time()
# ls <- foreach(i = ml, .packages = c('SpatialPosition'),
#               .combine = rbind, .inorder = FALSE) %dopar% {
# 
#   # mat <- CreateDistMatrix(knownpts = knownpts,
#   #                         unknownpts = unknownpts[sequence[i]:sequence[i+1]-1,],
#   #                         bypassctrl = TRUE)
#   # st <- stewart(knownpts = knownpts,
#   #               unknownpts = unknownpts[sequence[i]:sequence[i+1]-1,],
#   #               matdist = mat,
#   #               typefct = "exponential",
#   #               span = 20000, beta = 3,
#   #               var = "POPULATION")
#   mat <- CreateDistMatrix(knownpts = knownpts,
#                           unknownpts = i,
#                           bypassctrl = TRUE)
#   st <- stewart(knownpts = knownpts,
#                 unknownpts = i,
#                 matdist = mat,
#                 typefct = "exponential",
#                 span = 20000, beta = 3,
#                 var = "POPULATION")
# }
# 
# strt1 <- Sys.time()
# print(strt1 - strt0)
# stopCluster(cl)
# 
# ls2 <- spTransform(ls, old.proj)
# rasSt <- rasterStewart(ls2)
# a <- rasterToContourPoly(r = rasSt,
#                          breaks=c(10,20,50,100,200,500,1000,2000,5000,17588),
#                          mask = coco)
# 
# # Cartographie
# bks <- sort(unique(c(a$min, a$max)))
# # Display the map
# library(cartography)
# dev.off()
# opar <- par(mar = c(0,0,1.2,0))
# choroLayer(spdf = a,
#            df = a@data,
#            var = "center", legend.pos = "topleft",
#            breaks = bks, border = NA,
#            lwd = 0.2,
#            legend.title.txt = "Population potentielle\ndans un voisinage de 20 km\n(en milliers)",
#            legend.values.rnd = 0)
# layoutLayer(title = "Pop",
#             south = TRUE,
#             sources = "", author = "")
# par(opar)
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
# 
# i3 <- iter(data.frame(x=1:3, y=10), by='row')
# nextElem(i3)
# nextElem(i3)
# nextElem(i3)
