###############################
########### Optimal route: NPP

##############################
############### mkNPPMap  ----
mkNPPMap <- function(pts, startPts, surv, bout = 550000) {
  
  # if(!exists("dm") | nrow(dm)!=nrow(pts)) dm  <- distm(pts)
  startInd   <- which.min(distVincentySphere(startPts, pts))
  
  grd         <- data.frame(lon = pts[,1], lat =  pts[,2])
  grd$reached <- NA
  grd$reached[startInd] <- 0
  
  grd$surv <- surv
  grd$surv[is.na(grd$surv)] <- 0
  
  grd$survToGetThere <- NA
  grd$survToGetThere[!is.na(grd$reached)] <- 1
  
  grd$done <- F
  
  s = 1
  
  # plot(landC, col = "grey90", border = "grey70")
  # points(grd[,1:2], cex = 1.5, pch = 16, col = rev(terrain.colors(200, alpha = 0.5))[approx(c(0,1), c(1, 200), xout = grd$surv)$y])
  
  while(any(is.na(grd$reached)))
  {
    
    simFrom           <- which.max(grd$survToGetThere*!grd$done*!is.na(grd$reached))
    grd$done[simFrom] <- T
    
    neighbours        <- which(dm[simFrom,] < bout & dm[simFrom,] > 0)
    
    grd$surv[neighbours] <- surv[neighbours]
    tmp1                 <- (grd$surv[neighbours]*grd$survToGetThere[simFrom]>grd$survToGetThere[neighbours] | is.na(grd$reached[neighbours])) & !grd$done[neighbours]
    
    grd$reached[neighbours[tmp1]] <- simFrom
    grd$survToGetThere[neighbours[tmp1]] <- grd$surv[neighbours[tmp1]]*grd$survToGetThere[simFrom]
    
    # with(grd[neighbours[tmp1],], arrows(pts[reached,1], pts[reached,2], lon, lat, length = 0.1))
    
    if(all(surv[!grd$done]==0) | sum(grd$survToGetThere*!grd$done*!is.na(grd$reached), na.rm = T)==0) break
    
    s = s + 1
    
  }
  
  grd
  
}

# Test
# grd <- mkNPPMap(pts = pts, startP = startP, surv = surv)
# 
# plot(r0)
# points(pts[,1], pts[,2], pch = 16, cex = 1.5, col = viridis::viridis(200)[approx(range(grd$survToGetThere, na.rm = T), c(1,200),
#                                                                                  xout = grd$survToGetThere)$y])


##############################
############### NPPTrack  ----

NPPTrack   <- function(grd, endPts) {
  
  endInd   <- which.min(distVincentySphere(endPts, grd[,1:2]))
  
  trk <- grd[endInd,c(1,2,3,5)]
  repeat {
    if(trk$reached[nrow(trk)]==0) break
    trk <- rbind(trk, grd[trk$reached[nrow(trk)],c(1,2,3,5)])
  }
  trk[nrow(trk):1,]
}

# Test
# trk <- NPPTrack(grd, endP[1,])
# lines(trk$lon, trk$lat, lwd = 3, col = "orange")
