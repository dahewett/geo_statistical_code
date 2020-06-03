

install.packages("plyr")
install.packages("dplyr")
install.packages("spdep")
install.packages("GISTools")
install.packages("raster")
install.packages("maptools")
install.packages("rgdal")
install.packages("spatstat")
install.packages("sp")
install.packages("spatstat")
install.packages("spgwr")
install.packages("tmap")
install.packages("ggmap")

#Geog 418/518 Final Project
library(plyr)
library(dplyr)
library(spdep)
library(GISTools)
library(raster)
library(maptools)
library(rgdal)
library(spatstat)
library(sp)
library(spgwr)
library(tmap)
library(ggmap)


#Set working directory
setwd("/csc/geog418/project")


#Reading in particulate matter dataset
pm25 <- read.csv("PM25.csv") #Read in PM2.5 data
#Select only columns 1 and 2
pm25 <- pm25[,1:2]
#Change the column names 
colnames(pm25) <- c("POSTALCODE", "PM25")

#Reading in postal code shapefile
postalcodes <- shapefile("BC_Postal_Codes") #Read in related postal code data
#Join PM2.5 data with postal code data using the POSTALCODE column
pm25.spatial <- merge(postalcodes,pm25,by = "POSTALCODE")
#Plot the points on a map
plot(pm25.spatial)
#Examine the first several rows
head(pm25.spatial)
#You should notice that the dataset contains NA's, so these need to be removed.
pm25.spatial <- pm25.spatial[!is.na(pm25.spatial$PM25),]


#Reading in the income dataset
income <- read.csv("Income.csv") #Read in census income data  
#Change the column names
colnames(income) <- c("DAUID", "Income") #Select only ID and Income columns

#Read in the dissemination tract shapefile
census.tracts <- shapefile("BC_DA.shp") 


#Merge the income dataset and the DA shapefile
income.tracts <- merge(census.tracts,income, by = "DAUID") 
#Remove any NA's from the merged dataset
income.tracts <- income.tracts[!is.na(income.tracts$Income),]

#Study area map
tmap_mode("view")
data("World")
## tmap mode set to interactive viewing
tm_basemap("OpenStreetMap.Mapnik") +tm_shape(income.tracts)+tm_borders("purple", lwd = .4) +tm_scale_bar()+tm_shape(pm.income)+tm_symbols(size=0.01,alpha = 0,border.col = "Green", border.alpha = 0.7) +tm_style("classic", legend.only = TRUE)


memory.limit(size=50000)
#Create choropleth map of income
med.income <- income.tracts$Income
shades <- auto.shading(med.income, n=6, cols = brewer.pal(6, "Blues"))

choropleth(income.tracts, med.income, shades, border="transparent")  #map the data with associated colours
  choro.legend(-122.7, 49.495, shades) #add a legend (you might need to change the location)


#Perform a spatial intersection on the PM2.5 and Income data
pm.income <- intersect(pm25.spatial,income.tracts)
#Observe the result
head(pm.income)
pm.income$DAUID <- as.numeric(pm.income$DAUID)
pm.income$PM25 <- as.numeric(pm.income$PM25)
#Aggregate the the multiple PM2.5 values for each DA. Here the mean function is used.
pm.income <- aggregate(pm.income$PM25~pm.income$DAUID, FUN=mean)
#Change the column names
colnames(pm.income) <- c("DAUID", "PM25")
#Remove any NA's
pm.income <- na.omit(pm.income)

#Seeing as how the datasets are not properly merged, perform another merge to have PM and income together
pm.income.poly <- merge(income.tracts,pm.income,by = "DAUID")
pm.income.poly <- na.omit(pm.income.poly)
#Remove unwanted columns
pm.income.poly <- pm.income.poly[,-(2:23)]
#Observe the result. Are there still NA's? If so, apply following line to get rid of them.
pm.income.poly <- pm.income.poly[!is.na(pm.income.poly$PM25),]

#Create choropleth map of PM25
avg.pm <- pm.income.poly$PM25
shades <- auto.shading(avg.pm, n=6, cols = brewer.pal(6, 'Reds'))
choropleth(income.tracts, avg.pm, shades, border ="transparent") #map the data with associated colours
choro.legend(-122.7, 49.496, shades) #add a legend (you might need to change the location)


# Testing Spatial autocorrelation of income

#Create spatial neigh weight matrix
income.tracts <- income.tracts[,-(2:23)]
income.nb <- poly2nb(income.tracts)

plot(income.tracts)
plot(income.nb, coordinates(income.tracts), add = TRUE, col = "red")

income.lw <- nb2listw(income.nb, zero.policy = TRUE, style = "W")
print.listw(income.lw, zero.policy = TRUE)

income.a <- income.tracts$Income

income.lagged.means = lag.listw(income.lw, income.a, zero.policy = TRUE)
shades2 <- auto.shading(income.lagged.means, n=6, cols = brewer.pal(6, 'Oranges'))
choropleth(income.tracts, income.lagged.means,shades2, border="transparent")
choro.legend(-122.7, 49.49, shades2)

#Global Moran's I
mi <- moran.test(income.a, income.lw, zero.policy = TRUE)
mi

#To contextualize your Moran's I value, retrieve range of potential Moran's I values.
moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(income.lw)

#Perform the Z-test
#You can get the necessary values from your mi object resulting from your Moran's I test above.
#For example, the Moran's I value is the first value in the output mi, so you call mi$estimate[1] to get the value.
z=((mi$estimate[1]-mi$estimate[2])/(mi$estimate[3]))
z


#Local Moran's I
lisa.test <- localmoran(income.a, income.lw)
lisa.test
lisa.test <- na.omit(lisa.test)
#Create a choropleth map of the LISA values.
lisa.shades <- auto.shading(c(lisa.test[,1],-lisa.test[,1]),cols=brewer.pal(5,"PRGn"))
choropleth(income.tracts, lisa.test[,1],shading=lisa.shades, border="transparent")
choro.legend(-122.7, 49.49,lisa.shades,fmt="%6.2f")


#Create a Moran's I scatterplot
moran.plot(income.a, income.lw, zero.policy=NULL, spChk=NULL, labels=NULL, xlab="Median Income", 
           ylab="Spatially Lagged Median Income", quiet=NULL)


tm_layout(basemaps = c('OpenStreetMap')) + tm_scale_bar()+ tm_compass()





######Linear Regression##########

head(pm.income.poly)
#Plot income and PM2.5 from the pm.income.poly dataset you created
plot(pm.income.poly$Income~pm.income.poly$PM25)
#Notice that there are a lot of 0's in this dataset. If you decide to remove them, use the following line:
pm.income.poly <-  pm.income.poly[pm.income.poly$PM25 != 1, ]
#Now plot the data again
plot(pm.income.poly$Income~pm.income.poly$PM25)

#Perform a linear regression on the two variables. You should decide which one is dependent.
lm.model <- lm(pm.income.poly$Income~pm.income.poly$PM25)
#Add the regression model to the plot you created
abline(lm.model)
#Get the summary of the results
summary(lm.model)

#You want to determine if the model residuals are spatially clustered. 
#First obtain the residuals from the model
model.resids <- as.data.frame(residuals.lm(lm.model))
#Then add the residuals to your spatialpolygon dataframe
pm.income.poly$residuals <- residuals.lm(lm.model)
#Observe the result to make sure it looks correct
head(pm.income.poly)

#Now, create choropleth map of residuals
resids <- pm.income.poly$residuals
shades <- auto.shading(resids, n=6, cols = brewer.pal(6, 'Greens'))
choropleth(income.tracts, resids, shades, border ="transparent") #map the data with associated colours
choro.legend(-122.7, 49.49, shades) #add a legend (you might need to change the location)


## Global Moran's I 

# Create Neigbourhood Weights Matrix
# Queen's Neigbour
head(pm.income.poly)
tracts <- poly2nb(pm.income.poly)
# Create the spatial weighted neighbour list with queens
tracts.lw <- nb2listw(tracts, zero.policy = TRUE, style = "W")

mi <- moran.test(pm.income.poly$residuals, tracts.lw, zero.policy = TRUE)
mi

#To contextualize your Moran's I value, retrieve range of potential Moran's I values.
moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(tracts.lw)

#Perform the Z-test
#You can get the necessary values from your mi object resulting from your Moran's I test above.
#For example, the Moran's I value is the first value in the output mi, so you call mi$estimate[1] to get the value.
z=((mi$estimate[1]-mi$estimate[2])/(mi$estimate[3]))
z

#Local Moran's I
lisa.test <- localmoran(pm.income.poly$residuals, tracts.lw)
lisa.test
lisa.test <- na.omit(lisa.test)
#Create a choropleth map of the LISA values.
lisa.shades <- auto.shading(c(lisa.test[,1],-lisa.test[,1]),cols=brewer.pal(5,"PRGn"))
choropleth(pm.income.poly, lisa.test[,1],shading=lisa.shades, border ="transparent")
choro.legend(-122.7, 49.49,lisa.shades,fmt="%6.2f")


tracts.lw <- na.omit(tracts.lw)
pm.income.poly$residuals <- na.omit(pm.income.poly$residuals)
#Create a Moran's I scatterplot
moran.plot(pm.income.poly$residuals, tracts.lw, zero.policy=NULL, spChk=NULL, labels=NULL, xlab="Residuals", 
           ylab="Spatially Lagged Residuals", quiet=NULL)



####Geographically Weighted Regression
#The first thing you need to do is to add the polygon coordinates to the spatialpolygondataframe.
#You can obtain the coordinates using the "coordinates" function from the sp library
pm.income.poly.coords <- sp::coordinates(pm.income.poly)
#Observe the result
head(pm.income.poly.coords)
#Now add the coordinates back to the spatialpolygondataframe
pm.income.poly$X <- pm.income.poly.coords[,1]
pm.income.poly$Y <- pm.income.poly.coords[,2]
head(pm.income.poly)

###Determine the bandwidth for GWR: this will take a while
GWRbandwidth <- gwr.sel(pm.income.poly$Income~pm.income.poly$PM25, 
                        data=pm.income, coords=cbind(pm.income.poly$X,pm.income.poly$Y),adapt=T) 

###Perform GWR on the two variables with the bandwidth determined above
###This will take a looooooong while
gwr.model = gwr(pm.income.poly$Income~pm.income.poly$PM25, 
                data=pm.income.poly, coords=cbind(pm.income.poly$X,pm.income.poly$Y), 
                adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 

#Print the results of the model
gwr.model

#Look at the results in detail
results<-as.data.frame(gwr.model$SDF)
head(results)

#Now for the magic. Let's add our local r-square values to the map
pm.income.poly$localr <- results$localR2

#Create choropleth map of r-square values
local.r.square <- pm.income.poly$localr
shades <- auto.shading(local.r.square, n=6, cols = brewer.pal(6, 'Oranges'))
choropleth(income.tracts, local.r.square, shades, border = "transparent") #map the data with associated colours
choro.legend(-122.7, 49.49,shades) #add a legend (you might need to change the location)

#Time for more magic. Let's map the coefficients
pm.income.poly$coeff <- results$pm.income.poly.PM25

#Create choropleth map of the coefficients
local.coefficient <- pm.income.poly$coeff
shades <- auto.shading(local.coefficient, n=6, cols = brewer.pal(6, 'Oranges'))
choropleth(income.tracts, local.coefficient, shades,border="transparent") #map the data with associated colours
choro.legend(3864000, 1965000, shades) #add a legend (you might need to change the location)

plot(income.tracts)
plot(pm.income.poly)

pm.income.poly.ext <- as.matrix(extent(pm.income.poly))

window <- as.owin(list(xrange=pm.income.poly.ext[1,], yrange = pm.income.poly.ext[2,]))
pm.income.ppp <- ppp(x=pm.income.poly$X, y = pm.income.poly$Y, window = window)

nearestNeighbour <- nndist(pm.income.ppp)
##Convert the nearestNeighbor object into a dataframe.
nearestNeighbour=as.data.frame(as.numeric(nearestNeighbour))
##Change the column name to "Distance"
colnames(nearestNeighbour) = "Distance"

##Calculate the nearest neighbor statistic to test for a random spatial distribution.
#mean nearest neighbour
nnd = sum((nearestNeighbour$Distance))
nnd = nnd/nrow(nearestNeighbour)

#mean nearest neighbour for random spatial distribution
r.nnd = 1/(2*sqrt(nrow(nearestNeighbour)/0.44502))

d.nnd = 1.07453/sqrt(nrow(nearestNeighbour)/0.44502)

r = nnd/r.nnd
sd.nnd = 0.26136/sqrt((nrow(nearestNeighbour)*nrow(nearestNeighbour))/0.44502)
z = (nnd-r.nnd)/sd.nnd
z

# k function
k.fun <- Kest(pm.income.ppp, correction = "Ripley")
plot(k.fun)

#use simulation to test the point pattern against CSR
k.fun.e <- envelope(pm.income.ppp, Kest, nsim = 99, correction = "Ripley")
plot(k.fun.e)

#quadrat analysis

quads <- 10

qcount = quadratcount(pm.income.ppp, nx = quads, ny = quads)

plot(pm.income.ppp, pch = "+", cex = 0.5)
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