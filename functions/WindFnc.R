###############################
########### Optimal route: Wind


##############################
############## mkWindMap  ----


mkWindMap   <- function(pts, startPts, lookup = wlk, airspeed = 10, tmstep = 3) {
  
  # if(!exists("dm") | nrow(dm)!=nrow(pts)) dm  <- distm(pts)
  startInd   <- which.min(distVincentySphere(startPts, pts))
  
  grd                       <- data.frame(lon = pts[,1], lat =  pts[,2])
  grd$reachedFrom           <- NA
  grd$reachedAt             <- NA
  grd$reachedFrom[startInd] <- 0
  grd$reachedAt[startInd]   <- 0
  
  
  # plot(landC, col = "grey90", border = "grey90")
  # with(grd[is.na(grd$reachedFrom) || grd$reachedFrom>=0,], points(lon, lat, pch = 21, bg = "white", col = "grey80"))

  s = 0
  
  while(any(is.na(grd$reachedFrom))) {
    
    sitesTG <- as.numeric(which(!is.na(grd$reachedFrom) & grd$reachedAt == s))
    
      # points(pts[sitesTG,1], pts[sitesTG,2])
    
    if(length(sitesTG)>0) {
      
      kt = s  
      
      tab  <-  rbindlist(apply(cbind(c(sitesTG)), 1,  function(x) {
        data.table(simFrom    = x,
                   neighbours = as.numeric(which(dm[x,] < 280000 & dm[x,] > 0)),
                   lon        = as.numeric(pts[x,1]),
                   lat        = as.numeric(pts[x,2]),
                   time0      = as.numeric(0),
                   time       = as.numeric(NA),
                   reached    = FALSE)}))
      
      tab <- tab[is.na(grd[tab$neighbours,"reachedAt"]) | grd[tab$neighbours,"reachedAt"] > kt+1,]
      
      # points(pts[tab$neighbours,1], pts[tab$neighbours,2], pch = 16, cex = 0.5)
      
      if(nrow(tab)>0) {
        
        while(any(!tab$reached)) {
        
          tab[!tab$reached, c('lon', 'lat', 'time0') := as.data.table(t(apply(tab[!tab$reached,], 1, function(x) {
            
            d  <- bearingRhumb(pts[unlist(x)[1],], pts[unlist(x)[2],])
            Vg <- max(apply(lookup(matrix(c(unlist(x)[3], unlist(x)[4]), ncol = 2)), 1, function(w) NCEP.M.Groundspeed(unlist(w)[1], unlist(w)[2], d, airspeed)$groundspeed), na.rm = T)
            
            if(is.infinite(abs(Vg)) || is.na(Vg) || Vg<1) Vg <- 1
            
            Tt <- (distVincentySphere(unlist(x)[3:4], pts[unlist(x)[2],])/1000)/(Vg*3.6)
            Dp <- destPoint(unlist(x)[3:4], d, Vg*tmstep*60*60)
            
            c(Dp[,1], Dp[,2], time0 = Tt)
            
          })))]
          
          tab[tab$time0<=tmstep & !tab$reached, time    := kt+1]
          tab[tab$time0<=tmstep & !tab$reached, reached := TRUE]
          
          kt <- kt + 1
          
          if(kt>(s+10)) {
            tab[!tab$reached, time    := kt+1]
            tab[!tab$reached, reached := TRUE]
            break
          }
          
        }
        
        tabOut <- as.data.frame(rbindlist(lapply(split(tab, tab$neighbours), function(x) x[which.min(x$time),])))        
        
        grd[tabOut[,2], "reachedFrom"] <- ifelse(!is.na(grd[tabOut[,2],"reachedAt"]) & tabOut$time>grd[tabOut[,2],"reachedAt"], grd[tabOut[,2],"reachedFrom"], tabOut$simFrom)
        grd[tabOut[,2], "reachedAt"]   <- ifelse(!is.na(grd[tabOut[,2],"reachedAt"]) & tabOut$time>grd[tabOut[,2],"reachedAt"], grd[tabOut[,2],"reachedAt"], tabOut$time)
        
        # points(pts[!is.na(grd$reachedAt),], pch = 16, cex = 1, col = rev(heat.colors(max(grd$reachedAt, na.rm = T)))[grd$reachedAt[!is.na(grd$reachedAt)]], lwd = 2)
      }      
    }
    
    s = s+1
    
  }
  
  grd
}


# # Test
# grd <- mkWindMap(pts, startP, startT = as.POSIXct("2016-05-01"), 10, 3)
# 
# stopCluster(cl)
# 
# plot(landC, col = "grey90", border = "grey90")
# with(grd[is.na(grd$reachedFrom) || grd$reachedFrom>=0,], points(lon, lat, pch = 21, cex = 1.5, col = "grey80",
#                                                                 bg = viridis::viridis(200)[approx(c(0, max(grd$reachedAt, na.rm = T)), c(1,200),
#                                                                                                     xout = grd$reachedAt)$y]))
# points(startP[1], startP[2], pch = 21, bg = "white", cex = 2)


##############################
############## windTrack  ----

windTrack   <- function(grd, endP) {
  
  endInd   <- which.min(distVincentySphere(endP, grd[,1:2]))
  
  trk <- data.frame(grd[endInd,1:2], time = grd$reachedAt[endInd], cell = grd$reachedFrom[endInd])
  repeat {
    trk <- rbind(trk, data.frame(grd[trk[nrow(trk),4],1:2], time = grd$reachedAt[trk[nrow(trk),4]], cell = grd$reachedFrom[trk[nrow(trk),4]]))
    if(trk[nrow(trk),3]==0) break
  }
  trk[nrow(trk):1,]
}

## test
# trk <- windTrack(grd, endP[5,])
# lines(trk$lon, trk$lat, lwd = 3, col = "orange")


