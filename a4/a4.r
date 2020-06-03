install.packages("rgdal")
install.packages("gstat")
install.packages("sp")
install.packages("spatstat")
install.packages("maptools")
install.packages("raster")
install.packages("tmap")

library(rgdal)
library(gstat)
library(sp)
library(spatstat)
library(maptools)
library(raster)
library(tmap)

setwd("/csc/geog418/a4")

#################################################
##Prepare Pollution Data

#DATASET 1
#Read the pollution csv dataset.
ozone = read.csv("OZONE_PICKDATA_2016-4-30.csv", header = T, sep = ",")

#DATASET 2
#Read the monitoring station spatial dataset as an OGR data object.
monitor = readOGR(dsn = ".", layer = "airmonitoringstations")
#Extract the monitoring stations for the South Coast (SC)
SC.monitor = monitor[monitor$AIRBASIN %in% c("South Coast"),]
#Reproject the data to a suitable projection. Here we use a UTM projection because of the scale of the analysis. 
SC.monitor.t = spTransform(SC.monitor, CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))

#DATASET 3
#Read the California Air Basin spatial dataset.
Ca.AirBasin = readOGR(dsn = ".", layer = "CaAirBasin")
#Extract the South Coast air basin from the spatial dataset. 
SC.AirBasin = Ca.AirBasin[Ca.AirBasin$NAME %in% c("South Coast"),] 
#Reproject the South Coast air basin spatial dataset to match the projeciton of the monitoring station dataset.  
SC.AirBasin.t = spTransform(SC.AirBasin, CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))


##################################################################
##Process Pollution Data
#You need to represent each location with a single value in order to perform statistical analyses.

#Examine the first several rows of the ozone dataset. 
head(ozone)

#Looking at the date and hour columns, you can see that we need to process the data
#to get summary statistics.

#Calculate the mean and max ozone level for each site for all readings.
mean.ozone = aggregate(value ~ site, ozone, mean)
max.ozone = aggregate(value ~ site, ozone, max)

#Join the mean and max ozone values to their respective monitoring stations. In doing so, you will need to rename the 
#first column of the monitoring data to site in order to have a unique name to match the two datasets.
names(SC.monitor.t)[1] ="site"  

#Merge the the monitoring station shapefile with the ozone data using the site column.  
mrg.tab.mean <- sp::merge(SC.monitor.t, mean.ozone, by = "site", all.x = FALSE) 
mrg.tab.max <- sp::merge(SC.monitor.t, max.ozone, by = "site", all.x = FALSE)

#Create a max and a mean spatialPointDataFrame. 
ozone.mean.spdf = na.omit(mrg.tab.mean)
ozone.max.spdf = na.omit(mrg.tab.max)

# Load and observe ozone data
tm_layout(basemaps = c('OpenStreetMap')) + tm_scale_bar()+ tm_compass()+ tm_shape(SC.AirBasin.t) + tm_polygons() +
  tm_shape(ozone.mean.spdf) +
  tm_dots(col="value", palette = "RdBu", auto.palette.mapping = FALSE,
          title="Sampled Ozone \n(in ppm)", size=0.7) + tm_legend(legend.outside=TRUE)

tm_compass()
#################################################
##Spatial Interpolation with Thiessen Polygons

# Create a tessellated surface
th  <-  as(dirichlet(as.ppp(ozone.mean.spdf)), "SpatialPolygons")

# The dirichlet function does not carry over projection information
# requiring that this information be added manually
proj4string(th) <- proj4string(ozone.mean.spdf)

# The tessellated surface does not store attribute information
# from the point data layer. We'll use the over() function (from the sp
# package) to join the point attributes to the tesselated surface via
# a spatial join. The over() function creates a dataframe that will need to
# be added to the `th` object thus creating a SpatialPolygonsDataFrame object
th.z     <- over(th, ozone.mean.spdf, fn=mean)
th.spdf  <-  SpatialPolygonsDataFrame(th, th.z)

# Finally, we'll clip the tessellated  surface to the Texas boundaries
th.clp   <- raster::intersect(SC.AirBasin.t,th.spdf)

# Map the data
tm_layout(basemaps = c('OpenStreetMap')) + tm_scale_bar()+ tm_shape(th.clp) + 
  tm_polygons(col="value", palette="RdBu", auto.palette.mapping=FALSE,
              title="Predicted Ozone \n(in ppm)") +
  tm_legend(legend.outside=TRUE)


#################################################
##Spatial Interpolation with Polynomial Trends
# Define the 1st order polynomial equation

f.1 <- as.formula(value ~ X + Y) 

# Add X and Y to P
ozone.mean.spdf$X <- coordinates(ozone.mean.spdf)[,1]
ozone.mean.spdf$Y <- coordinates(ozone.mean.spdf)[,2]

# Run the regression model
lm.1 <- lm( f.1, data=ozone.mean.spdf)

# Use the regression model output to interpolate the surface
dat.1st <- SpatialGridDataFrame(grd, data.frame(var1.pred = predict(lm.1, newdata=grd))) 

# Clip the interpolated raster to Texas
r   <- raster(dat.1st)
r.m <- mask(r, SC.AirBasin.t)

# Plot the map
tm_shape(r.m) + tm_layout(basemaps = c('OpenStreetMap')) +
  tm_raster(n=10, palette="RdBu", auto.palette.mapping=FALSE, 
            title="Predicted Ozone \n(in ppm)") +
  tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)



# Define the 2nd order polynomial equation
f.2 <- as.formula(value ~ X + Y + I(X*X)+I(Y*Y) + I(X*Y))

# Add X and Y to P
ozone.mean.spdf$X <- coordinates(ozone.mean.spdf)[,1]
ozone.mean.spdf$Y <- coordinates(ozone.mean.spdf)[,2]

# Run the regression model
lm.2 <- lm( f.2, data=ozone.mean.spdf)

# Use the regression model output to interpolate the surface
dat.2nd <- SpatialGridDataFrame(grd, data.frame(var1.pred = predict(lm.2, newdata=grd))) 

# Clip the interpolated raster to Texas
r   <- raster(dat.2nd)
r.m <- mask(r, SC.AirBasin.t)

# Plot the map
tm_shape(r.m) + tm_layout(basemaps = c('OpenStreetMap'))+
  tm_raster(n=10, palette="RdBu", auto.palette.mapping=FALSE,
            title="Predicted Ozone \n(in ppm)") +
  tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)


#################################################-------1
##Spatial Interpolation with IDW

# Create an empty grid where n is the total number of cells
grd <- as.data.frame(spsample(ozone.mean.spdf, "regular", n=50000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object

proj4string(grd) <- proj4string(SC.monitor.t)
P.idw <- gstat::idw(value ~ 1, ozone.mean.spdf, newdata=grd, idp=3.0)
r       <- raster(P.idw)
r.m     <- mask(r, SC.AirBasin.t)

tm_shape(r.m) + tm_layout(basemaps = c('OpenStreetMap')) + tm_scale_bar()+
  tm_raster(n=10,palette = "RdBu", auto.palette.mapping = FALSE,
            title="Predicted Ozone \n(in ppm)") + 
  tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)

#################################################
# Leave-one-out validation routine
IDW.out <- vector(length = length(ozone.mean.spdf))
for (i in 1:length(ozone.mean.spdf)) {
  IDW.out[i] <- idw(value ~ 1, ozone.mean.spdf[-i,], ozone.mean.spdf[i,], idp=3.0)$var1.pred
}

# Plot the differences
OP <- par(pty="s", mar=c(4,3,0,0))
plot(IDW.out ~ ozone.mean.spdf$value, asp=1, xlab="Observed", ylab="Predicted", pch=16,
     col=rgb(0,0,0,0.5))
abline(lm(IDW.out ~ ozone.mean.spdf$value), col="red", lw=2,lty=2)
abline(0,1)
par(OP)
sqrt( sum((IDW.out - ozone.mean.spdf$value)^2) / length(ozone.mean.spdf))


#################################################
# Implementation of a jackknife technique to estimate a confidence interval at each unsampled point.
# Create the interpolated surface
img <- gstat::idw(value~1, ozone.mean.spdf, newdata=grd, idp=3.0)
n   <- length(ozone.mean.spdf)
Zi  <- matrix(nrow = length(img$var1.pred), ncol = n)

# Remove a point then interpolate (do this n times for each point)
st <- stack()
for (i in 1:n){
  Z1 <- gstat::idw(value~1, ozone.mean.spdf[-i,], newdata=grd, idp=3.0)
  st <- addLayer(st,raster(Z1,layer=1))
  # Calculated pseudo-value Z at j
  Zi[,i] <- n * img$var1.pred - (n-1) * Z1$var1.pred
}

# Jackknife estimator of parameter Z at location j
Zj <- as.matrix(apply(Zi, 1, sum, na.rm=T) / n )

# Compute (Zi* - Zj)^2
c1 <- apply(Zi,2,'-',Zj)            # Compute the difference
c1 <- apply(c1^2, 1, sum, na.rm=T ) # Sum the square of the difference

# Compute the confidence interval
CI <- sqrt( 1/(n*(n-1)) * c1)

# Create (CI / interpolated value) raster
img.sig   <- img
img.sig$v <- CI /img$var1.pred 

# Clip the confidence raster to Texas
r <- raster(img.sig, layer="v")
r.m <- mask(r, SC.AirBasin.t)

# Plot the map
tm_shape(r.m) + tm_layout(basemaps = c('OpenStreetMap'))+ tm_raster(n=7,title="Confidence Interval / Zj") +
  tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)




#################################################
##Spatial Interpolation with Kriging

f.1 <- as.formula(value ~ X + Y) 
var.smpl <- variogram(f.1, ozone.mean.spdf, cloud = FALSE) #, cutoff=1000000, width=89900)
dat.fit  <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(psill=1.58e-05, model="Sph", range=13, nugget=0))
plot(var.smpl, dat.fit, main = "plot1")


var.smpl2 <- variogram(f.1, ozone.mean.spdf, cloud = FALSE) #, cutoff=1000000, width=89900)
dat.fit2  <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(psill=1.58e-05, model="Exp", range=13, nugget=0))
plot(var.smpl2, dat.fit2, main = "plot2")

var.smpl3 <- variogram(f.1, ozone.mean.spdf, cloud = FALSE) #, cutoff=1000000, width=89900)
dat.fit3  <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(psill=1.58e-05, model="Gau", range=13, nugget=0))
plot(var.smpl3, dat.fit3, main = "plot3")

var.smpl4 <- variogram(f.1, ozone.mean.spdf, cloud = FALSE) #, cutoff=1000000, width=89900)
dat.fit4  <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(psill=0.38e-05, model="Gau", range=13, nugget=0))
plot(var.smpl4, dat.fit4, main = "plot4")

var.smpl5 <- variogram(f.1, ozone.mean.spdf, cloud = FALSE) #, cutoff=1000000, width=89900)
dat.fit5  <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(psill=1.58e-05, model="Gau", range=15, nugget=0))
plot(var.smpl5, dat.fit5, main = "plot5")


# Define the trend model
f.1 <- as.formula(value ~ X + Y) 

# Perform the krige interpolation (note the use of the variogram model
# created in the earlier step)
st <- stack() 
for (i in 1:n){
  Z1 <- krige( f.1, ozone.mean.spdf, grd, dat.fit)
  st <- addLayer(st,raster(Z1,layer=1))
  Zi[,i] <- n*img$var1.pred - (n-1)* Z1$var1.pred
}
# Jackknife estimator of parameter Z at location j
Zj <- as.matrix(apply(Zi, 1, sum, na.rm=T) / n )

# Compute (Zi* - Zj)^2
c1 <- apply(Zi,2,'-',Zj)            # Compute the difference
c1 <- apply(c1^2, 1, sum, na.rm=T ) # Sum the square of the difference

# Compute the confidence interval
CI <- sqrt( 1/(n*(n-1)) * c1)

# Create (CI / interpolated value) raster
img.sig   <- img
img.sig$v <- CI /img$var1.pred 

# Clip the confidence raster to Texas
r <- raster(img.sig, layer="v")
r.m <- mask(r, SC.AirBasin.t)






dat.krg <- krige( f.1, ozone.mean.spdf, grd, dat.fit)
dat.krg2 <- krige( f.1, ozone.mean.spdf, grd, dat.fit2)
dat.krg3 <- krige( f.1, ozone.mean.spdf, grd, dat.fit3)
dat.krg4 <- krige( f.1, ozone.mean.spdf, grd, dat.fit4)
dat.krg5 <- krige( f.1, ozone.mean.spdf, grd, dat.fit5)


# Convert kriged surface to a raster object for clipping
r <- raster(dat.krg)
r2 <- raster(dat.krg2)
r3 <- raster(dat.krg3)
r4 <- raster(dat.krg4)
r5 <- raster(dat.krg5)

r.m <- mask(r, SC.AirBasin.t)
r.m2 <- mask(r2, SC.AirBasin.t)
r.m3 <- mask(r3, SC.AirBasin.t)
r.m4 <- mask(r4, SC.AirBasin.t)
r.m5 <- mask(r5, SC.AirBasin.t)


# Plot the map
m <- tm_shape(r.m) + tm_layout(basemaps = c('OpenStreetMap'))+
  tm_raster(n=10, palette="RdBu", auto.palette.mapping=FALSE, 
            title="Predicted Ozone \n(in ppm)") +
  tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)

m2 <- tm_shape(r.m2) + 
  tm_raster(n=10, palette="RdBu", auto.palette.mapping=FALSE, 
            title="Predicted Ozone \n(in ppm)") +
  tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)
m3 <- tm_shape(r.m3) + tm_layout(basemaps = c('OpenStreetMap'))+
  tm_raster(n=10, palette="RdBu", auto.palette.mapping=FALSE, 
            title="Predicted Ozone \n(in ppm)") +
  tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)
m4 <- tm_shape(r.m4) + 
  tm_raster(n=10, palette="RdBu", auto.palette.mapping=FALSE, 
            title="Predicted Ozone \n(in ppm)") +
  tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)
m5 <- tm_shape(r.m5) + 
  tm_raster(n=10, palette="RdBu", auto.palette.mapping=FALSE, 
            title="Predicted Ozone \n(in ppm)") +
  tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)

multiplot(m4,m5, cols= 2 )

img

r   <- raster(dat.krg, layer="var1.var")
r2   <- raster(dat.krg2, layer="var2.var")
r3   <- raster(dat.krg3, layer="var3.var")
r4   <- raster(dat.krg4, layer="var4.var")
r5   <- raster(dat.krg5, layer="var5.var")

r.m <- mask(r, SC.AirBasin.t)
r.m2 <- mask(r2, SC.AirBasin.t)
r.m3 <- mask(r3, SC.AirBasin.t)
r.m4 <- mask(r4, SC.AirBasin.t)
r.m5 <- mask(r5, SC.AirBasin.t)

qq <- tm_shape(r.m) + 
  tm_raster(n=7, palette ="Reds",
            title="Variance map \n") +tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)
qq2 <- tm_shape(r.m2) + 
  tm_raster(n=7, palette ="Reds",
            title="Variance map \n") +tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)
qq3 <- tm_shape(r.m3) + 
  tm_raster(n=7, palette ="Reds",
            title="Variance map \n") +tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)
qq4 <- tm_shape(r.m4) + 
  tm_raster(n=7, palette ="Reds",
            title="Variance map \n") +tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)
qq5 <- tm_shape(r.m5) + 
  tm_raster(n=7, palette ="Reds",
            title="Variance map \n") +tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)

multiplot(qq,qq4,qq2,qq5,qq3, cols= 3 )

r   <- sqrt(raster(dat.krg, layer="var1.var")) * 1.96
r2   <- sqrt(raster(dat.krg2, layer="var1.var")) * 1.96
r3   <- sqrt(raster(dat.krg3, layer="var1.var")) * 1.96
r4   <- sqrt(raster(dat.krg4, layer="var1.var")) * 1.96
r5   <- sqrt(raster(dat.krg5, layer="var1.var")) * 1.96

r.m <- mask(r, SC.AirBasin.t)
r.m2 <- mask(r2, SC.AirBasin.t)
r.m3 <- mask(r3, SC.AirBasin.t)
r.m4 <- mask(r4, SC.AirBasin.t)
r.m5 <- mask(r5, SC.AirBasin.t)

pp <- tm_shape(r.m) + 
  tm_raster(n=7, palette ="Reds",
            title="95% CI map") +tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)
pp2 <- tm_shape(r.m2) + 
  tm_raster(n=7, palette ="Reds",
            title="95% CI map") +tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)
pp3 <- tm_shape(r.m3) + 
  tm_raster(n=7, palette ="Reds",
            title="95% CI map") +tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)
pp4 <- tm_shape(r.m4) + 
  tm_raster(n=7, palette ="Reds",
            title="95% CI map") +tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)
pp5 <- tm_shape(r.m5) + 
  tm_raster(n=7, palette ="Reds",
            title="95% CI map") +tm_shape(ozone.mean.spdf) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)

multiplot(pp4,pp5, cols= 2 )