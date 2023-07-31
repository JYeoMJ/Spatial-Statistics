knitr::opts_chunk$set(echo=FALSE, warnings=FALSE)


# This code is modified from:
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
library("tidyr")
library("reshape2")
library("gtable")
library("dplyr")


# Predicting at the BAUs when the large scale trend has covariates

# load meuse.grid and convert to SpatialPixelsDataFrame

data("meuse", package="sp")
coordinates(meuse) = ~ x + y
data("meuse.grid", package = "sp")
coordinates(meuse.grid) = ~ x + y
gridded(meuse.grid) = TRUE

# Remove any common fields in the data and the BAUs (recall that in FRK all covariates
# need to be with the BAUs)

meuse$soil <- meuse$dist <- meuse$ffreq <- NULL

# formula for SRE
f <- log(zinc) ~ 1 + sqrt(dist)

# Run FRK
S <- FRK(f = f,                # formula
         data = list(meuse),    # data (just meuse)
         BAUs = meuse.grid,     # the BAUs
         regular = 0)           # irregularly placed basis functions

# Predict at all the BAUs
Pred <- predict(S, obs_fs = FALSE)  
Pred_df<-SpatialPolygonsDataFrame_to_df(sp_polys = Pred,vars=c("mu","sd"))



g1<-ggplot() +
  geom_polygon(aes(x=x, y=y, fill=mu,group=id),
               data=Pred_df)
g1<-g1 + xlab("Easting") + ylab("Northing")

g2<-ggplot() +
  geom_polygon(aes(x=x, y=y, fill=sd, group=id),
               data=Pred_df)
g2<-g2 + xlab("Easting") + ylab("Northing")

grid.arrange(g1,g2,nrow=1)


# Example where both observed data and prediction locations are
# supplied in one data frame, with response at prediction locations
# missing (NA)

# Reload meuse data and assume first 10 observations are missing
data("meuse", package = "sp")
meuse[1:10, "zinc"] <- NA

# Create a complete-data data frame
meuse2 <- subset(meuse, !is.na(zinc))
meuse2 <- meuse2[, c("x", "y", "zinc")]
coordinates(meuse2) <- ~x+y

# Create BAUs around all prediction and observation 
# points (after removing data from field)

meuse$zinc <- NULL
coordinates(meuse) <- c("x", "y")
meuse.grid2 <- BAUs_from_points(meuse)

# Run the function FRK with the data frame and created BAUs.  
f <- log(zinc) ~ 1 + sqrt(dist)
S <- FRK(f = f, 
         data = list(meuse2), 
         BAUs = meuse.grid2, 
         regular = 0)
Pred <- predict(S, obs_fs = FALSE)
data(meuse)
plot(Pred$mu[1:10],log(meuse$zinc[1:10]),xlab="Predicted",ylab="Actual",cex.lab=1.5)
abline(c(0,1))


# Construct model for meuse data with square grid, and consider
# prediction grid with different support

# Initialise seed to ensure reproducibility
set.seed(1)

# Load meuse data and convert to sp object
data("meuse", package = "sp")
coordinates(meuse) =  ~ x + y

# Construct the BAUs
GridBAUs1 <- auto_BAUs(manifold = plane(),    # we are on the plane
                       type = "grid",         # gridded BAUs (not hex)
                       cellsize = c(100, 100), # 100m x 100m
                       data = meuse,          # data for boundary
                       nonconvex_hull = TRUE, # nonconvex boundary
                       convex = -0.05)        # convex parameter (see
                                              # INLA::inla.nonconvex.hull)



G <- auto_basis(manifold = plane(),  # we are on the plane
                data = meuse,        # data around which to construct basis
                regular = 0,         # irregularly placed basis
                nres = 3,            # three resolutions
                type = "bisquare")   # bisquare functions

## Plot basis functions
dev.new()
print(show_basis(G) +            # main call
      xlab("Easting (m)") +      # x-label
      ylab("Northing (m)") +     # y-label
      coord_fixed())             # fixed asp. ratio



# formula (just intercept model)
f <- log(zinc) ~ 1


S <- FRK(f = f,                     
         data = list(meuse),        
         BAUs = GridBAUs1,          
         basis = G,                 
         average_in_BAU = FALSE)    



# Now look at prediction of averages on a 600m x 600m grid

Pred_regions <- auto_BAUs(manifold = plane(),     # we are on the plane
                          cellsize = c(600, 600),  # 600m x 600m
                          type = "grid",          # grid
                          data = meuse,           # meuse data for boundary
                          convex = -0.05)         # convex parameter (see
                                                  # INLA::inla.nonconvex.hull)

## Now predict over these regions
Pred <- predict(S, newdata = Pred_regions,   # prediction regions
                obs_fs = FALSE)              # Case 2 in paper (default)

## Plot results

## First converg SpatialPolygonsDataFrame to data frame
Pred_regions_df <- SpatialPolygonsDataFrame_to_df(sp_polys = Pred,      # the sp obj
                                                  vars = c("mu", "sd")) # the variables to put
                                                                        # into the data frame


g1<-ggplot() +
  geom_polygon(aes(x=x, y=y, fill=mu,group=id),
               data=Pred_regions_df)
g1<-g1 + xlab("Easting") + ylab("Northing")

g2<-ggplot() +
  geom_polygon(aes(x=x, y=y, fill=sd, group=id),
               data=Pred_regions_df)
g2<-g2 + xlab("Easting") + ylab("Northing")

grid.arrange(g1,g2,nrow=1)

## ----predict_areal_data1----------------------------------------------------------------------

# Now consider an example where the observed data do not have
# point support

set.seed(1)      # Fix seed
data("meuse.grid", package = "sp") # load meuse data grid
data("meuse", package = "sp")      # load meuse data

## Generate observations with large spatial support
meuse_pols <- NULL  # initialise
offset <- 150       # half-size of box (150 m)

## for each meuse data point
for (i in 1:nrow(meuse)) {
    this_meuse <- meuse[i, ]
    meuse_pols <- rbind(meuse_pols,  # create a box aroun it of sides 300 m in length
                        data.frame(x = c(this_meuse$x - offset,
                                         this_meuse$x + offset,
                                         this_meuse$x + offset,
                                         this_meuse$x - offset),
                                   y = c(this_meuse$y - offset,
                                         this_meuse$y - offset,
                                         this_meuse$y + offset,
                                         this_meuse$y + offset),
                                   id = i,
                                   zinc = this_meuse$zinc))
}

## Convert to SpatialPolygonsDataFrame
meuse_pols <- df_to_SpatialPolygons(meuse_pols,coords = c("x", "y"), keys = "id", proj = CRS())
meuse_pols <- SpatialPolygonsDataFrame(meuse_pols,
                                       data.frame(row.names = row.names(meuse_pols),
                                                  zinc = meuse$zinc))
coordnames(meuse_pols) <- c("x", "y")

## Convert meuse and meuse.grid to sp object
coordinates(meuse) <- ~x + y
coordinates(meuse.grid) <- ~x+y



## ----conditional_simulation----------------------------------------------------------------------
# Conditionally simulate field on the BAUs and aggregate to polygons 
# to generate realistic observations

# First construct the BAUs
GridBAUs <- auto_BAUs(manifold = plane(),      # 2D plane
                      cellsize = c(50, 50),    # BAU cellsize
                      type = "grid",           # grid (not hex)
                      data = meuse,            # data around which to create BAUs
                      convex = -0.05)          # border buffer factor
GridBAUs$fs <- 1

## Now fit variogram to meuse data and conditionally simulate
f <- log(zinc) ~ 1
lzn.vgm <- variogram(f, meuse)
lzn.fit <- fit.variogram(lzn.vgm, model = vgm(1, "Sph", 900, 1))
lzn.sim <- krige(formula = f, meuse,
                 GridBAUs,
                 model = lzn.fit,
                 nsim=1, nmax = 30)

# Observe field using large spatial support. In the below we run FRK to just generate
# the mapping matrix, which we then use to assign values to the meuse polygons
S <- FRK(f = f,  ## Just get out C Matrix
         data = list(meuse_pols),
         BAUs = GridBAUs,
         regular = 0,
         n_EM = 1, print_lik = FALSE)

meuse_pols$zinc <- as.numeric(S@Cmat %*% exp(as.numeric(lzn.sim$sim1)) +
                              0.1*rnorm(length(meuse_pols)))



## Generate the basis functions
G <- auto_basis(manifold = plane(),   # 2D plane
                data = meuse,         # meuse data
                nres = 3,             # number of resolutions
                type = "bisquare",    # type of basis function
                regular = 0)          # place irregularly in domain

## Construct and fit the SRE model

S <- FRK(f = f,
         data=list(meuse_pols),
         BAUs=GridBAUs,
         basis=G)
         
GridBAUs2 <- predict(S, obs_fs = FALSE)


## Create data frames for plotting
BAUs_df <- data.frame(GridBAUs2)
Obs_df <- SpatialPolygonsDataFrame_to_df(sp_polys = meuse_pols,   # BAUs to convert
                                         vars = c("zinc"))  # fields to extract

## Now plot the footprints
g1 <- ggplot() +
    geom_path(data = Obs_df,
              aes(x, y, group = id),
              colour = "black") +
    coord_fixed() +
    xlab("Easting (m)") + ylab("Northing (m)")


## Now plot the prediction
g2 <- ggplot() +                         # Use a plain theme
    geom_tile(data = BAUs_df,            # Draw BAUs
              aes(x, y, fill = mu),      # Colour <-> Mean
              colour = "light grey") +   # Border is light grey
    scale_fill_distiller(palette = "Spectral", name = "pred")  +  # Spectral palette
    coord_fixed() +                               # fix aspect ratio
    xlab("Easting (m)") + ylab("Northing (m)")    # axes labels

## Now plot the prediction error
g3 <- ggplot() +                          # Similar to above but with s.e.
    geom_tile(data = BAUs_df,
              aes(x, y, fill = sqrt(var)),
              colour = "light grey") +
    scale_fill_distiller(palette = "BrBG",
                         name = "s.e.") +
    coord_fixed() +
    xlab("Easting (m)") + ylab("Northing (m)")

dev.new()
plot(g1)


#library(gtable)
#p2 <- ggplotGrob(g2)
#p3 <- ggplotGrob(g3)
#g <- cbind(p2, p3, size="first")
#g$heights <- unit.pmax(p2$heights, p3$heights)
#grid.newpage()
#grid.draw(g)

grid.arrange(g2, g3, nrow = 1) # Plot

x.seq <- y.seq <- seq(0, 1, length.out = 100)
BAUs <- expand.grid(x = x.seq, y = y.seq)
coordinates(BAUs) = ~ x + y
gridded(BAUs) <- TRUE

BAUs_df <- coordinates(BAUs) %>% as.data.frame

## ----Poisson_generate_data----------------------------------------------------------------------

set.seed(2020)

# Define some complicated smooth trigonometric field

#f <- function(x, y, a = 0, b = 0, l = 1) exp(-l * sqrt((x - a)^2 + (y - b)^2))
logistic <- function(x, L = 1, k=1, x0 = 0) L / (1 + exp(-k * (x - x0)))
smooth_Y_process <- function(x, y) {
    a <- 4 + 2 * sin(5 * x) + 2 * cos(4 * y) + 2 * sin(3 * x * y)+
    sin(7 * x) + cos(9 * y) +
    sin(12 * x) + cos(17 * y) + sin(14 * x) * cos(16 * y)  +
    sin(23 * x) + cos(22 * y) + sin(24 * x) * cos(26 * y) +
    x + y + x^2 + y^2 + x * y + 
    sin(10 * pi * x * y) + cos(20 * x * y - 3) + sin(40 * x * y)
    
  logistic(a, x0 = 2, L = 6, k = 0.35)
}
  
BAUs_df <- BAUs_df %>% 
  dplyr::mutate(Y = smooth_Y_process(x, y) + rnorm(length(BAUs), mean = 0, sd = 0.2)) 

## Compute the true mean process, mu, over the BAUs
BAUs_df <- dplyr::mutate(BAUs_df, mu = exp(Y))

## Subsample n of the BAUs to act as observation locations, and simulate data
n <- 750
Poisson_simulated <- sample_n(BAUs_df, n) %>% 
  dplyr::mutate(Z = rpois(n, lambda = mu)) %>%
  dplyr::select(x, y, Z)  

## scalar matrix for fine scale variation, and converts to SpatialPixelsDF
BAUs$fs <- rep(1, length(BAUs)) 

## Convert Poisson_simulated to Spatial* object
coordinates(Poisson_simulated) <- ~ x + y



S <- FRK(f = Z ~ 1, data = list(Poisson_simulated), 
                      nres = 3, BAUs = BAUs, 
                      response = "poisson", 
                      link = "log", 
                      K_type = "precision", method = "TMB", est_error = FALSE) 
pred_list<- predict(S, type = c("link", "mean"))

## Predictions, uncertainty, and data
plot_list <- plot(S, pred_list$newdata, 
                  labels_from_coordnames = FALSE)

# Plots

g1<-ggplot(BAUs_df) + geom_tile(aes(x, y, fill = mu)) + 
  scale_fill_distiller(palette = "Spectral") + 
  labs(x = expression(s[1]), y = expression(s[2])) +
  theme_bw() + coord_fixed() + 
  scale_x_continuous(breaks=c(0.25, 0.75), expand = c(0, 0)) + 
  scale_y_continuous(breaks=c(0.25, 0.75), expand = c(0, 0)) + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 19), 
        legend.text = element_text(size = 16), 
        plot.title = element_text(hjust = 0.5, size = 19))
g1


g2<-plot_list$Z + labs(title = plot_list$Z$labels$fill) + 
  labs(fill = "", colour = "") + 
  scale_x_continuous(breaks=c(0.25, 0.75), expand = c(0, 0)) + 
  scale_y_continuous(breaks=c(0.25, 0.75), expand = c(0, 0)) + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 19), 
        legend.text = element_text(size = 16), 
        plot.title = element_text(hjust = 0.5, size = 19))
g2

g3<-plot_list$p_Y + labs(title = plot_list$p_Y$labels$fill) + 
  labs(fill = "", colour = "") + 
  scale_x_continuous(breaks=c(0.25, 0.75), expand = c(0, 0)) + 
  scale_y_continuous(breaks=c(0.25, 0.75), expand = c(0, 0)) + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 19), 
        legend.text = element_text(size = 16), 
        plot.title = element_text(hjust = 0.5, size = 19))
g3


