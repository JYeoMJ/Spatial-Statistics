
knitr::opts_chunk$set(echo=TRUE, message=FALSE, warnings=FALSE)


# Parts of this code are taken from:
# 
# Zammit-Mangion, A., and Cressie, N. (2021). FRK: An R Package for Spatial 
# and Spatio-Temporal Prediction with Large Datasets. 
# Journal of Statistical Software, 98(4), 1--48.
#

library("sp")
library("FRK")
library("gstat")
library("ggplot2")
library("gridExtra")
library("grid")
library("INLA")


data(meuse)
coordinates(meuse) = ~x+y # changes meuse into an sp object

# BAUs from hexagonal grid

HexPols_df <- auto_BAUs(manifold = plane(),
                        cellsize = 100,
                        type = "hex",
                        data = meuse,
                        nonconvex_hull = FALSE)
plot(HexPols_df)



G <- auto_basis(manifold = plane(),data=meuse,nres = 2,regular=2,prune=0.1,type = "bisquare")
show_basis(G,ggplot()) + geom_point(data=data.frame(meuse),aes(x,y))

# Initialise seed to ensure reproducibility
set.seed(1)

# formula for SRE model (just an intercept here)
f <- log(zinc) ~ 1

# Call FRK()
S <- FRK(f = f,                # formula
         data = list(meuse),   # list of data objects (just meuse in this case)
         basis=G,
         BAUs=HexPols_df)          

summary(S)   # Print out a summary of the returned (fitted) SRE model

# Get predictions (for the BAUs by default)

GridBAUsPred<- predict(S, obs_fs = FALSE)

# For plotting convert the SpatialPolygonsDataFrame GridBAUs1 to a data frame

Pred_df <- SpatialPolygonsDataFrame_to_df(sp_polys = GridBAUsPred,vars=c("mu","sd"))

# Set coordinate reference system for predictions
# 4326 is the EPSG code for WGS84

# Create plot for the means

g1<-ggplot() +
  geom_polygon(
    aes(x=x, y=y, fill=mu, group=id),
    data = Pred_df)

# Change x and y axis labels

g1<-g1 + xlab("Easting") + ylab("Northing")

# Plot for standard deviation

g2<-ggplot() +
  geom_polygon(
    aes(x=x, y=y, fill=sd, group=id),
    data = Pred_df)

g2<-g2 + xlab("Easting") + ylab("Northing")

# Combine the above two plots on one graph

grid.arrange(g1,g2,nrow=1)

