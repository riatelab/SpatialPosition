## ---- fig.width=7, fig.height=5-----------------------------------------------
library(SpatialPosition)
library(sf)
data(hospital)

# Compute potentials (accessibility)
potentials <- stewart(
  knownpts = hospital, 
  varname = "capacity",
  typefct = "exponential", 
  span = 1000, 
  beta = 3,
  resolution = 50,
  mask = paris, 
  returnclass = "sf"
)

isopotentials <- isopoly(x = potentials, mask = paris, returnclass = "sf")

lab <- paste0(round(isopotentials$min,0),' to ', 
              round(isopotentials$max,0))
par(mar = c(4,2,2,1))
plot(st_geometry(isopotentials), col = heat.colors(8))
legend(x = "topright", legend = lab, fill = heat.colors(8), 
       cex = 0.7, title = "Potentials")
  
mtext("Global Accessibility to Public Hospitals", side = 3,cex = 1.5)
mtext(text = "Potential nb. of beds
      distance function: exponential, span = 1 km, beta = 3",
      side = 1, line = 1)   

## ---- fig.width=5, fig.height=5-----------------------------------------------
library(raster)
row.names(hospital)
catchReilly <- reilly(knownpts = hospital, varname = "capacity",
                      typefct = "exponential", span = 750, beta = 2,
                      resolution = 50, mask = paris)
# Create a raster
rasterCatch <- rasterReilly(x = catchReilly, mask = paris)

par(mar = c(4,2,2,1))
# Plot the raster and add the points
plotReilly(x = rasterCatch)
plot(st_geometry(hospital), pch = 20, add = TRUE)

mtext("Catchment Areas of Public Hospitals", side = 3,cex = 1.5)
mtext(text = "distance function: exponential, span = 0.75 km, beta = 2",
      side = 1, line = 0) 

## ---- fig.width=5, fig.height=5-----------------------------------------------
catchHuff <- huff(knownpts = hospital, varname = "capacity",
                  typefct = "exponential", span = 750, beta = 2,
                  resolution = 50, mask = paris)

# Create a raster
rasterCatch <- rasterHuff(x = catchHuff, mask = paris)

# Plot the raster and add the points
par(mar = c(4,2,2,1))
plotHuff(x = rasterCatch)
plot(st_geometry(hospital), pch = 20, col = "red", add = TRUE)

mtext("Probabilistic Catchment Areas \nof Public Hospitals", 
      side = 3,cex = 1.5, line=-1.5)
mtext(text = "distance function: exponential, span = 0.75 km, beta = 2",
      side = 1, line = 0) 

