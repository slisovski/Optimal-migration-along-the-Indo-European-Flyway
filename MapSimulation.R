### Map Simulations
source("functions/setup.R")
source("functions/WindFnc_db.R")
source("functions/NPPFnc_db.R")

wnPath  <- "~/Wind/" ## path to ERA interim wind data
fls  <- list.files(wnPath, pattern = "*.nc")
flsT <- data.frame(name = fls, Year = NA, Month = NA)

for(i in 1:length(fls)) {
  nf <- nc_open(paste0(wnPath, fls[i]))
  tm <- as.POSIXct("1900-01-01", tz = "GMT") + nf$dim$time$vals*60*60
  flsT[i,2] <- as.numeric(format(tm[1], "%Y"))
  flsT[i,3] <- as.numeric(format(tm[1], "%m"))
  nc_close(nf)
}
flsWT <- flsT[order(flsT$Year, flsT$Month),]

npPath  <- "~/MODIS/VHP_SM_SMN/" # Path to MODIS data
npFiles <- list.files(npPath, pattern = "tif")
year    <- as.numeric(ifelse(nchar(npFiles)==34, substring(npFiles, 17, 20), substring(npFiles, 18, 21)))
week    <- as.numeric(ifelse(nchar(npFiles)==34, substring(npFiles, 21, 23), substring(npFiles, 22, 24)))

date0 <- cbind(year, week)
npDate <- as.POSIXct(apply(date0, 1, function(x) {
  tm <- seq(as.POSIXct(paste0(x[1], "-01-01")), as.POSIXct(paste0(x[1], "-12-31")), by = "day")
  w  <- which(x[2]==as.numeric(format(tm, "%U")))
  mean(tm[w])
}), origin = "1970-01-01", tz = "GMT")


##############################
################ Cluster -----

cl <- parallel::makeForkCluster(2)
doParallel::registerDoParallel(cl)

##############################
############### simPaths  ----

tm0  <- seq(as.POSIXct("2012-01-01 15:00:00", tz = "GMT"), as.POSIXct("2018-01-01 15:00:00", tz = "GMT"), by = "2 days")
tmS  <- tm0[as.numeric(format(tm0, "%j"))%in%c(75:154)]
tmA  <- tm0[as.numeric(format(tm0, "%j"))%in%c(197:276)]

tracks <- array(dim = c(100, 3, nrow(stend[[2]]), length(tmS), nrow(stend[[1]]), 2, 2))

for(seas in 1:2) {
  
  if(seas==1) {
    startP <- stend[[2]]
    endP   <- stend[[1]]
  } else {
    startP <- stend[[1]]
    endP   <- stend[[2]]
  }
  
  for(t in 1:ifelse(seas==1, length(tmS), length(tmA))) {
    
      if(seas==1) {
        startT <- as.POSIXct(tmS[t], origin = "1970-01-01")
        indT   <- which(tmS==startT)
      } else {
        startT <- as.POSIXct(tmA[t], origin = "1970-01-01")
        indT   <- which(tmA==startT)
      }
    
      cat(sprintf('\rSeas %d - time %d of %d',
                 seas, indT, length(tmA)))
    
      if(all(is.na(tracks[, , ,indT, ,seas, 2]))) {  
      
            #### Wind Lookup ----
            flT <- flsWT[flsWT$Year==as.numeric(format(startT, "%Y")) & flsWT$Month==as.numeric(format(startT, "%m")),"name"]

                u0_1 <- crop(brick(paste0(wnPath, flT), varname = "u", level = 1), ext)
                u0_2 <- crop(brick(paste0(wnPath, flT), varname = "u", level = 2), ext)
                u0_3 <- crop(brick(paste0(wnPath, flT), varname = "u", level = 3), ext)

                v0_1 <- crop(brick(paste0(wnPath, flT), varname = "v", level = 1), ext)
                v0_2 <- crop(brick(paste0(wnPath, flT), varname = "v", level = 2), ext)
                v0_3 <- crop(brick(paste0(wnPath, flT), varname = "v", level = 3), ext)

                ind  <- .bincode(as.numeric(startT), as.numeric(as.POSIXct(substr(names(u0_1), 2, nchar(names(u0_1)[1])), "%Y.%m.%d.%H.%M", tz = "GMT")))

              geom   <- SpatialPointsDataFrame(coordinates(u0_1),
                                               data = data.frame(U1 = u0_1[[ind]][], U2 = u0_2[[ind]][], U3 = u0_3[[ind]][],
                                                                 V1 = v0_1[[ind]][], V2 = v0_2[[ind]][], V3 = v0_3[[ind]][]),
                                               proj4string = sp::CRS("+init=epsg:4326"))
            ####################

            grdL <- foreach::foreach(i = 1:nrow(startP), .packages = c("SearchTrees")) %dopar% {
              tree <- createTree(sp::coordinates(geom))
              wlk <- function(p) {
                matrix(geom[c(knnLookup(tree, newdat = p, k=1)),]@data[,-1], ncol = 2, byrow = FALSE)
              }
              mkWindMap(pts = pts, startPts = startP[i,], lookup = wlk, airspeed = 10, tmstep = 3)
            }

            trkL <- lapply(grdL, function(x) {
              apply(endP, 1, function(y) {
                windTrack(x, y)
              })
            })

            for(s in 1:nrow(startP)) {
              for(e in 1:nrow(endP)) {
                tracks[1:nrow(trkL[[s]][[e]]), , s, indT, e, seas, 1]  <- as.matrix(trkL[[s]][[e]][,1:3])
              }
            }
            
            #### MODIS Lookup ----
            ext <- extent(c(range(pts[,1])+c(-1,1), range(pts[,2])+c(-1, 1)))
            
            if(!exists("r0") || (any(names(r0)!=c(substring(npFiles[max(which(npDate<=startT))], 1, nchar(npFiles[max(which(npDate<=startT))])-4),
                                    substring(npFiles[min(which(npDate>startT))], 1, nchar(npFiles[min(which(npDate>startT))])-4))))) {
            r0  <- aggregate(
                        crop(
                          brick(raster(paste0(npPath, npFiles[max(which(npDate<=startT))])), 
                                raster(paste0(npPath, npFiles[min(which(npDate>startT))]))), 
                            ext),
                              fact = 25, fun = max, na.rm = T)
            }
            
            r   <- calc(r0, function(x) approx(npDate[c(max(which(npDate<=startT)), min(which(npDate>startT)))], x, xout = startT)$y)

            survParm <- c(20, 0.3)
            npp      <- seq(0, 1, length.out = 101)
            survY    <- sigmoid(npp, a = survParm[1], b = survParm[2])
              # plot(npp, survY, type = "o")
            survC    <- approxfun(x = npp, y = survY)
            
            # r[] <- survC(ifelse(r[]<0, NA, r[]))
            #   plot(r)
            surv     <- survC(extract(r, pts))
            surv[is.na(surv)] <- 0
            ####################
            
            grdL <- foreach::foreach(i = 1:nrow(startP)) %dopar% {
              mkNPPMap(pts, startPts = startP[i,], surv, bout = 550000)
            }
            
            trkL <- lapply(grdL, function(x) {
              apply(endP, 1, function(y) {
                tryCatch(NPPTrack(x, y), error = function(e) NULL)
              })
            })
            
            for(s in 1:nrow(startP)) {
              for(e in 1:nrow(endP)) {
                if(!is.null(trkL[[s]][[e]])) tracks[1:nrow(trkL[[s]][[e]]), , s, indT, e, seas, 2]  <- as.matrix(trkL[[s]][[e]][,c(1,2,4)])
              }
            }
      
      save(tracks, file = "MapSimulation.RData")      
      }

  } ## end t
  
} ## end seas  

stopCluster(cl)