###PROJECTING AND FORMATTING DATA
#install packages 

install.packages("spatstat")
install.packages("rgdal")
install.packages("maptools")
install.packages("raster")
install.packages("sp")
install.packages("plyr")

#load required libraries
library(spatstat)
library(rgdal)
library(maptools)
library(raster)
library(sp)
library(plyr)

###ACQUIRE AND FORMAT DATA ###

#get data from here
#https://www.google.com/maps/d/viewer?mid=18u0QER64-OR_Kacg_EoKQpDUU5g&hl=en&ll=49.16837282831341%2C-122.63991487265628&z=10
setwd("/csc/geog418/a2")
#read in homicide events as shapefile, uses rgdal readOGR
km <- readOGR("v_kill.kml")
#clean up the columns
km$Name <- as.character(km$Name)
#create a year column
km$year <- as.numeric(substr(km$Name, nchar(km$Name)-4, nchar(km$Name)))
#remove one observation without a year
km <- km[complete.cases(km$year),]

#project to bc albers
kma <- spTransform(km, CRS("+init=epsg:3005"))
#kma <- spTransform(km, CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 
#                           +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

#add coordinates to the data
kma$x <- coordinates(kma)[,1]
kma$y <- coordinates(kma)[,2]

#check for and remove duplicated points
#check for duplicated points
#finds zero distance among points
zd <- zerodist(kma)
zd
#remove duplicates
kma <- remove.duplicates(kma)


#create an "extent" object which can be used to create the observation window for spatstat
kma.ext <- as.matrix(extent(kma)) 

#observation window
window <- as.owin(list(xrange = kma.ext[1,], yrange = kma.ext[2,]))

#create ppp oject from spatstat
kma.ppp <- ppp(x = kma$x, y = kma$y, window = window)

###KERNEL DENSITY ESTIMATION
#2D (gaussian) kernel, compare how bandwidth (sigma) selection influences the point density estimates
#since data are projected, sigma is represented in metres
#eps is the width and height of the pixels (1000m X 1000m)
#coerce to a SpatialGridDataFrame for plotting
kde.100 <- density(kma.ppp, sigma = 100, at = "pixels", eps = c(1000, 1000))
kde.SG <- as(kde.100, "SpatialGridDataFrame")
kde.500 <- density(kma.ppp, sigma = 500, at = "pixels", eps = c(1000, 1000))
kde.SG <- cbind(kde.SG, as(kde.500, "SpatialGridDataFrame"))
kde.1k <- density(kma.ppp, sigma = 1000, at = "pixels", eps = c(1000, 1000)) 
kde.SG <- cbind(kde.SG, as(kde.1k, "SpatialGridDataFrame"))
kde.5k <- density(kma.ppp, sigma = 5000, at = "pixels", eps = c(1000, 1000))
kde.SG <- cbind(kde.SG, as(kde.5k, "SpatialGridDataFrame"))

names(kde.SG) <- c("kde.100m", "kde.500m", "kde.1km", "kde.5km")
#plot
x11() #opens a new plot window
spplot(kde.SG)

#can see how the bandwidth selection influences the density estimates
summary(kde.SG)

#use cross-validation to get the bandwidth that minimizes MSE
bw.d <- bw.diggle(kma.ppp)
#plot the "optimal" bandwidth
plot(bw.d, ylim=c(-10, 10), main="Cross validation for murder events")

#density using the cross-validation bandwidth
kde.bwo <- density(kma.ppp, sigma = bw.d, at = "pixels", eps = c(1000, 1000))
plot(kde.bwo)

###NEAREST NEIGHBOUR
nearestNeighbour <- nndist(kma.ppp)
##Convert the nearestNeighbor object into a dataframe.
nearestNeighbour=as.data.frame(as.numeric(nearestNeighbour))
##Change the column name to "Distance"
colnames(nearestNeighbour) = "Distance"


##Calculate the nearest neighbor statistic to test for a random spatial distribution.
#mean nearest neighbour
nnd = sum(nearestNeighbour$Distance)/nrow(nearestNeighbour)


#mean nearest neighbour for random spatial distribution
r.nnd = 1/(2*sqrt(nrow(nearestNeighbour)/3884100000))
  
d.nnd = 1.07453/sqrt(nrow(nearestNeighbour)/3884100000)
  
r = nnd/r.nnd
sd.nnd = 0.26136/sqrt((nrow(nearestNeighbour)*nrow(nearestNeighbour))/3884100000)
z = (nnd-r.nnd)/sd.nnd
  

###K-FUNCTION 
#basic k-function
k.fun <- Kest(kma.ppp, correction = "Ripley")
plot(k.fun)

#use simulation to test the point pattern against CSR
k.fun.e <- envelope(kma.ppp, Kest, nsim = 99, correction = "Ripley")
plot(k.fun.e)

###QUADRAT ANALYSIS

##First, determine the number of qusdrats 
quads <- 10

qcount = quadratcount(kma.ppp, nx = quads, ny = quads)

plot(kma.ppp, pch = "+", cex = 0.5)
plot(qcount, add = T, col = "red")

qcount.df <- as.data.frame(qcount)

##Second, count the number of quadrats with a distinct number of points.
qcount.df = plyr::count(qcount.df,'Freq')
##Change the column names so that x=number of points and f=frequency of quadrats with x cells.
colnames(qcount.df) = c("x","f")

##Third, create new columns for total number of points and for fx^2.
qcount.df$TotPoints <- qcount.df$x * qcount.df$f
qcount.df$fx2 = (qcount.df$x)^2 * qcount.df$f
qcount.df$xfx2 = qcount.df$fx2 * qcount.df$f #adjusted for the count 

##Fourth, calculate the sum of each column, which you will use as inputs into the 
##formula for VMR.
f.sum = sum(qcount.df$f)
TotPoints.sum = sum(qcount.df$TotPoints) 
fx2.sum = sum(qcount.df$fx2) 
l

##Fifth, calculate VAR, MEAN, and VMR. ### OF WHICH VARIABLES? f.sum, TotPoints.Sum, fx2.sum?
VAR = sum(qcount.df$xfx2)/(sum(qcount.df$f)-1)
mean.points = TotPoints.sum/(quads*quads) 
VMR = VAR/mean.points

#Finally, perform the test statistic to test for the existence of a random spatial pattern.
chi.square = VMR*(100-1)

#KDE

kde.100 <- density(pm.income.ppp, sigma = 100, at = "pixels", eps = c(1000, 1000))
kde.SG <- as(kde.100, "SpatialGridDataFrame")

