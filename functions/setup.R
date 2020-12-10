### Setup

library(raster)
library(rgeos)
library(sp)
library(gdalUtils)
library(rgdal)
library(geosphere)
library(maptools)
  data("wrld_simpl")
library(ncdf4)
library(RNCEP)
library(data.table)
library(pracma)
library(parallel)
library(sf)
library(SearchTrees)
library(RPostgres)
library(rpostgis)
library(DBI)
library(foreach)

deg2rad = pi/180
rad2deg = 180/pi

set.seed(120)

#################
######### Map ----

e <- file.exists("~/Dropbox/Data/GeoDat/NaturalEarth/50m_physical/ne_50m_land/ne_50m_land.shp")

if(e){
  land <- suppressWarnings(readOGR("~/Dropbox/Data/GeoDat/NaturalEarth/50m_physical/ne_50m_land/ne_50m_land.shp"))
  # proj4string(land) <- proj4string(wrld_simpl)
} else {
  land <- wrld_simpl
}

xlim = c( 0, 100)
ylim = c(-5, 72)
ext  <- extent(c(xlim, ylim))

landC <- gIntersection(land, as(extent(c(xlim, ylim)), "SpatialPolygons"), byid = T)
# plot(landC)


proj <- "+proj=laea"
xlim <- extent(landC)[1:2]
ylim <- extent(landC)[3:4]

xcentre <- round(xlim[1] + diff(xlim)/2)
ycentre <- round(ylim[1] + diff(ylim)/2)

proj <- sprintf("%s +lon_0=%f +lat_0=%f +ellps=sphere", proj, xcentre, ycentre)

landP <- spTransform(land, CRS(proj))

circ <- gBuffer(spTransform(SpatialPoints(matrix(c(xcentre, ycentre), ncol = 2), proj4string = CRS(proj4string(landC))), CRS(proj)), width = 4050000, quadsegs = 40)

landCirc <- gIntersection(landP, circ)

# plot(landCirc)


#################
## Locations ----

poly              <- as(extent(landC), "SpatialPolygons")
proj4string(poly) <- proj4string(wrld_simpl)
xcentre           <- round(extent(poly)[1] + diff(extent(poly)[1:2])/2)
ycentre           <- round(extent(poly)[3] + diff(extent(poly)[3:4])/2)
proj              <- sprintf("%s +lon_0=%f +lat_0=%f +ellps=sphere", "+proj=laea", xcentre, ycentre)
polyP             <- spTransform(poly, CRS(proj))

HexPts  <-spsample(polyP, type="hexagonal", cellsize = 150000)
HexPols <- HexPoints2SpatialPolygons(HexPts)


# plot(landC)
pts <- project(coordinates(HexPts), proj = proj, inv = TRUE)
dm  <- distm(pts)
points(pts)


load("FlywayAnalysis/Results/StartEnd.RData")
breed   <- stend[[2]]
winter  <- stend[[1]]
