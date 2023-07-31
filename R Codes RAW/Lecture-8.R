# Load packages

library(sp)
library(gstat)
library(rgdal)
library(RColorBrewer)
library(classInt)

# Load data

data(meuse)
coordinates(meuse) <- ~x+y
meuse$lzn <- log(meuse$zinc)
data(meuse.grid)
coordinates(meuse.grid) <- ~x+y
meuse.grid <- as(meuse.grid, "SpatialPixelsDataFrame")

# Inverse distance weighting

idw.out <- gstat::idw(log(zinc) ~ 1, meuse, meuse.grid, idp=2.5)
plot(idw.out, col=heat.colors(100))
head(as.data.frame(idw.out), n=2)

spplot(meuse, "lzn", main="Log-Zinc Measurements", 
       col.regions = heat.colors(4),key.space="right")


# Ordinary kriging

vgm_log_zinc <- variogram(log(zinc) ~ 1, meuse)
fit_v <- fit.variogram(vgm_log_zinc, vgm(1, "Sph", 800, 1))
lz.ok <- krige(log(zinc)~1, meuse, meuse.grid, fit_v)
plot(lz.ok, col=heat.colors(100))


# Kriging variance

plot(lz.ok['var1.var'])
plot(meuse, add=TRUE, col='grey50')


# Filtering

err.fit <- fit.variogram(vgm_log_zinc, vgm(1, "Sph", 800, Err=1))
krige(log(zinc) ~ 1, meuse, meuse[1,], err.fit)
krige(log(zinc) ~ 1, meuse, meuse[1,], fit_v)


# Universal kriging

vt.fit2 <- fit.variogram.gls(log(zinc) ~ sqrt(dist), data=meuse[1:40,], 
                             model= vgm(1, "Exp", 300, 1) ,trace=FALSE)
lz.uk.gls <- krige(log(zinc)~sqrt(dist), meuse, meuse.grid, vt.fit2)
plot(lz.uk.gls, col=heat.colors(100))

