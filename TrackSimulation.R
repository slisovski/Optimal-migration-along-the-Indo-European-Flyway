### Track Simulations
source('functions/setup.R')
source('functions/WindFnc.R')
source('functions/NPPFnc.R')

#############################
############### getTracks  ----

gTab <- read.csv("/SummaryLocations.csv") ## Summary locations of geoloctor data (available upon request from Movebank)
ids <- unique(as.character(gTab[,1]))

flsSummary <- list.files("/RosefinchMigration/Geolocation/Results", pattern = "*Group_fit.RData") ## Fit data of geolocation analysis
flsID <- unlist(lapply(strsplit(flsSummary, "_"), function(x) x[[1]]))

SimTab <- data.frame(ID = ids, fls = flsSummary[match(ids, flsID)], dep01 = NA, arr01 = NA, dep02 = NA, arr02 = NA, 
                     depInd01 = NA, arrInd01 = NA, depInd02 = NA, arrInd02 = NA,
                     lon01 = NA, lat01 = NA, lon02 = NA, lat02 = NA, lon03 = NA, lat03 = NA)

for(i in 1:nrow(SimTab)) {
  tmp <- gTab[gTab[,1]==ids[i],]
  load(SimTab$fls[i])
  dts <- as.numeric(fit$model$time)
  
  if(any(tmp$Site==1 & tmp$Type==0)) {
    SimTab$dep01[i]    <- as.POSIXct(as.character(tmp$EndTime)[tmp$Site==1 & tmp$Type==0], tz = "GMT")
    SimTab$depInd01[i] <- min(which(dts>=SimTab$dep01[i]))
    
    SimTab$lon01[i]    <- tmp$Lon.50.[tmp$Site==1 & tmp$Type==0]
    SimTab$lat01[i]    <- tmp$Lat.50.[tmp$Site==1 & tmp$Type==0]
  }
  if(any(tmp$Type==2)) {
    SimTab$arr01[i]    <- as.POSIXct(as.character(tmp$StartTime)[min(which(tmp$Type==2))], tz = "GMT")
    SimTab$arrInd01[i] <- min(which(dts>=SimTab$arr01[i]))
    
    SimTab$lon02[i]    <- tmp$Lon.50.[min(which(tmp$Type==2))]
    SimTab$lat02[i]    <- tmp$Lat.50.[min(which(tmp$Type==2))]
  }
  if(any(tmp$Type==3)) {
    SimTab$dep02[i]    <- as.POSIXct(as.character(tmp$EndTime)[max(which(tmp$Type==2))], tz = "GMT")
    SimTab$depInd02[i] <- min(which(dts>=SimTab$dep02[i]))
    
    SimTab$lon03[i]    <- tmp$Lon.50.[max(which(tmp$Type==2))]
    SimTab$lat03[i]    <- tmp$Lat.50.[max(which(tmp$Type==2))]
  }
  if(any(tmp$Site>2 & tmp$Type==0)) {
    SimTab$arr02[i]    <- as.POSIXct(as.character(tmp$StartTime)[max(which(tmp$Site>2 & tmp$Type==0))], tz = "GMT")
    SimTab$arrInd02[i] <- min(which(dts>=SimTab$arr02[i]))
  }
}

SimTab$dep01 <- as.POSIXct(SimTab$dep01, origin = "1970-01-01", tz = "GMT")
SimTab$arr01 <- as.POSIXct(SimTab$arr01, origin = "1970-01-01", tz = "GMT")
SimTab$dep02 <- as.POSIXct(SimTab$dep02, origin = "1970-01-01", tz = "GMT")
SimTab$arr02 <- as.POSIXct(SimTab$arr02, origin = "1970-01-01", tz = "GMT")

SimTab$Pop   <- c(rep("FIN", 9), rep("GER", 4), rep("CZE",2), rep("BUL", 3), rep("SWE", 3)) 


#############################
############### simPaths  ----

trks <- data.frame(ID = NA, type = NA, season = NA, tm = NA, lon = NA, lat = NA, val = NA)

cl <- parallel::makeCluster(detectCores())
clusterEvalQ(cl, {
  library("data.table")
  library("geosphere")
  library("RNCEP")
  library("raster")
  NULL
})

clusterExport(cl, varlist = c("pts", "mkWindTable", "flsT"))

for(i in 1:nrow(SimTab)) {
  
  ### Autumn Migration ----
  startT   <- SimTab$dep01[i] - 2*24*60*60
  
  days     <- round(((as.numeric(SimTab$arr01[i]) - as.numeric(startT))/24/60/60)) + 20
  tmstep   <- 3
  airspeed <- 10
  
  clusterExport(cl, varlist= c("tmstep", "airspeed"))
  
  startP <- SimTab[i, c("lon01", "lat01")]
  endP   <- SimTab[i, c("lon02", "lat02")]
  
  clusterExport(cl, varlist = c("days", "startT"))
  
  null <- clusterEvalQ(cl, {
    lookupWnd  <- mkWindTable(as.POSIXct(startT, origin = "1970-01-01"), days, pts, flsT)
  })
  
  lookupNPP  <- mkNPPTable(as.POSIXct(startT, origin = "1970-01-01"), days, pts, flsNPPTab)
  
  tmSeq      <- seq(SimTab$dep01[i], SimTab$arr01[i], by = 48*60*60)
  
  
  for(t in as.numeric(tmSeq)) {
    
    cat(sprintf('\rInd %d of %d - Autumn: time %d of %d \b\b\b\b',
                i, nrow(SimTab), which(tmSeq==t), length(tmSeq)))
    
    ########### Wind Map
    startTm <- as.POSIXct(t, origin = "1970-01-01", tz = "GMT")
    clusterExport(cl, varlist= "startTm")
    grd <- mkWindMap(pts = pts, startP = startP, airspeed, tmstep)

    ########## Wind Track
    windTrk <- windTrack(grd, endP)
    trks <- rbind(trks, data.frame(ID = SimTab$ID[i], type = 1, season = 1, tm = rep(t, nrow(windTrk)), 
                                   lon = windTrk[,1], lat = windTrk[,2], val = max(windTrk$time)))
    
    ############## NP Map
    grdNPP  <- mkNPPMap(pts, lookupNPP, startP, t)
    NPPtrk  <- NPPTrack(grdNPP, endP)
    trks <- rbind(trks, data.frame(ID = SimTab$ID[i], type = 2, season = 1, tm = rep(t, nrow(NPPtrk)), 
                                   lon = NPPtrk[,1], lat = NPPtrk[,2], val = NPPtrk[nrow(NPPtrk),4]))
    
    gc()
  }

  ### Spring Migration ----
  if(!is.na(SimTab$dep02[i])) {
    
    startT <- SimTab$dep02[i] - 2*24*60*60
    
    days   <- round(((as.numeric(SimTab$arr02[i]) - as.numeric(startT))/24/60/60)) + 20
    days <- ifelse(is.na(days), 60, days)
    
    startP <- SimTab[i, c("lon03", "lat03")]
    endP   <- SimTab[i, c("lon01", "lat01")]
    
    clusterExport(cl, varlist = c("days", "startT"))

    null <- clusterEvalQ(cl, {
      lookupWnd  <- mkWindTable(as.POSIXct(startT, origin = "1970-01-01"), days, pts, flsT)
    })
    lookupNPP  <- mkNPPTable(as.POSIXct(startT, origin = "1970-01-01"), days, pts, flsNPPTab)
    
    tmSeq      <- seq(SimTab$dep02[i], as.POSIXct(ifelse(is.na(SimTab$arr02[i]), SimTab$dep02[i] + 45*24*60*60, SimTab$arr02[i]), 
                                                  origin = "1970-01-01", tz = "GMT"), by = 48*60*60) 
    
    for(t in as.numeric(tmSeq)) {
      
      cat(sprintf('\rInd %d of %d - Sring: time %d of %d \b\b\b\b\b\b',
                  i, nrow(SimTab), which(tmSeq==t), length(tmSeq)))
      
      ########### Wind Map
      startTm <- as.POSIXct(t, origin = "1970-01-01", tz = "GMT")
      clusterExport(cl, varlist= "startTm")
      grd <- mkWindMap(pts = pts, startP = startP, airspeed, tmstep)
      
      ########## Wind Track
      windTrk <- windTrack(grd, endP)
      trks <- rbind(trks, data.frame(ID = SimTab$ID[i], type = 1, season = 2, tm = rep(t, nrow(windTrk)), lon = windTrk[,1], lat = windTrk[,2], val = max(windTrk$time)))
      
      ############## NP Map
      grdNPP  <- mkNPPMap(pts, lookupNPP, startP, t)
      NPPtrk  <- NPPTrack(grdNPP, endP)

      trks <- rbind(trks, data.frame(ID = SimTab$ID[i], type = 2, season = 2, tm = rep(t, nrow(NPPtrk)), lon = NPPtrk[,1], lat = NPPtrk[,2], val = NPPtrk[nrow(NPPtrk),4]))
      
      gc()
    }
    
  }
  ######
  
}

stopCluster(cl)

save(trks, file = "TrackSimulation.RData")