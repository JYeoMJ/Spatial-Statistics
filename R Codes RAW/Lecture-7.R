# Load packages

library(sp)
library(gstat)
library(rgdal)
library(RColorBrewer)
library(classInt)

# Plot of Meuse river data

data(meuse)
coordinates(meuse) <- ~x+y
meuse$lzn <- log(meuse$zinc)
spplot(meuse, "lzn", main="Log-Zinc Measurements", 
       col.regions = heat.colors(5))


# Pairwise scatterplots with distance bins

hscat(log(zinc) ~ 1, meuse, (0:9)*100)


# Sample variogram 

vgm_log_zinc <- variogram(log(zinc) ~ 1, meuse)
plot(vgm_log_zinc)

# Sample variograms in different directions

plot(variogram(log(zinc) ~ 1, meuse, alpha=c(0, 45, 90, 135)))


# Plot of Smoky mountains data

smoky_data <- read.table("../../data/smoky/smoky.dat", header=TRUE, skip=1)
coordinates(smoky_data) <- ~easting + northing
smok <- readOGR("../../data/smoky", "smokypoly")
pal <- brewer.pal(4, "RdYlGn")
q4 <- classIntervals(smoky_data$ph, n=4, style="fisher")
spplot(smoky_data, "ph", cuts=q4$brks , col.regions=pal, colorkey=TRUE, 
       sp.layout=list("sp.polygons", smok), scales=list(draw=TRUE))


# Sample variogram

vgm_ph <- variogram(ph ~ 1, smoky_data, cutoff=120, width=10)
plot(vgm_ph)
head(vgm_ph, n=2)


# Show semivariograms

show.vgms(par.strip.text=list(cex=0.6))

head(vgm())


# Show a spherical semivariogram

sph1 <- vgm(psill=2, model='Sph', range=10, nugget=1)
class(sph1)
sph1
show.vgms(sill=2, range=10, models='Sph', nugget=1, max=13)


# Gaussian semivariogram

gauss1 <- vgm(psill=2, model='Gau', range=10, nugget=1)
gauss1


# Sum of semivariograms

gauss2 <- vgm(psill=2, model='Gau', range=10, nugget=0, add.to=sph1)
out <- variogramLine(gauss2, maxdist=15)
plot(out, type='l', main="Sum of Two Variograms")


# Fitting a parametric semivariogram model

fit_v <- fit.variogram(vgm_log_zinc, vgm(1, "Sph", 800, 1))
fit_v
attr(fit_v, "SSErr")
plot(vgm_log_zinc, model=fit_v)

