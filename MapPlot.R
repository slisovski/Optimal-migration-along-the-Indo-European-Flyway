### Track Simulations
load("MapSimulation.RData")
npp <- tracks

source("functions/setup.R")

tm0  <- seq(as.POSIXct("2012-01-01 15:00:00", tz = "GMT"), as.POSIXct("2017-01-01 15:00:00", tz = "GMT"), by = "2 days")
tmS  <- tm0[as.numeric(format(tm0, "%j"))%in%c(75:154)]
tmA  <- tm0[as.numeric(format(tm0, "%j"))%in%c(197:276)]

mat   <- array(0, dim = c(nrow(pts), nrow(pts), 2,2))
sList <- list(t1 = list(s1 = list(pop1 = list(dest1= list(), dest2 = list(), dest3 = list(), dest4 = list(), dest5 = list()), 
                                  pop2 = list(dest1= list(), dest2 = list(), dest3 = list(), dest4 = list(), dest5 = list()), 
                                  pop3 = list(dest1= list(), dest2 = list(), dest3 = list(), dest4 = list(), dest5 = list()), 
                                  pop4 = list(dest1= list(), dest2 = list(), dest3 = list(), dest4 = list(), dest5 = list()), 
                                  pop5 = list(dest1= list(), dest2 = list(), dest3 = list(), dest4 = list(), dest5 = list())),
                        s2 = list(pop1 = list(dest1= list(), dest2 = list(), dest3 = list(), dest4 = list(), dest5 = list()), 
                                  pop2 = list(dest1= list(), dest2 = list(), dest3 = list(), dest4 = list(), dest5 = list()), 
                                  pop3 = list(dest1= list(), dest2 = list(), dest3 = list(), dest4 = list(), dest5 = list()), 
                                  pop4 = list(dest1= list(), dest2 = list(), dest3 = list(), dest4 = list(), dest5 = list()), 
                                  pop5 = list(dest1= list(), dest2 = list(), dest3 = list(), dest4 = list(), dest5 = list()))), 
              t2 = list(s1 = list(pop1 = list(dest1= list(), dest2 = list(), dest3 = list(), dest4 = list(), dest5 = list()), 
                                  pop2 = list(dest1= list(), dest2 = list(), dest3 = list(), dest4 = list(), dest5 = list()), 
                                  pop3 = list(dest1= list(), dest2 = list(), dest3 = list(), dest4 = list(), dest5 = list()), 
                                  pop4 = list(dest1= list(), dest2 = list(), dest3 = list(), dest4 = list(), dest5 = list()), 
                                  pop5 = list(dest1= list(), dest2 = list(), dest3 = list(), dest4 = list(), dest5 = list())),
                        s2 = list(pop1 = list(dest1= list(), dest2 = list(), dest3 = list(), dest4 = list(), dest5 = list()), 
                                  pop2 = list(dest1= list(), dest2 = list(), dest3 = list(), dest4 = list(), dest5 = list()), 
                                  pop3 = list(dest1= list(), dest2 = list(), dest3 = list(), dest4 = list(), dest5 = list()), 
                                  pop4 = list(dest1= list(), dest2 = list(), dest3 = list(), dest4 = list(), dest5 = list()), 
                                  pop5 = list(dest1= list(), dest2 = list(), dest3 = list(), dest4 = list(), dest5 = list()))))

for(t in 1:2) {
  for(s in 1:2) {

    if(s==1) dates <- tmS else dates   <- tmA
    years   <- as.numeric(format(dates, "%Y"))

    trackS <- tracks[,,,,,s,t]

    for(pop in 1:dim(trackS)[5]) {
      
      if(t==1) {
        
        if(s==1) {
          tm    <- lapply(1:5, function(y) t(apply(trackS[,3,,,y], 2:3, max, na.rm = T)))
          surv  <- lapply(tm,  function(y) apply(y, 2, function(x) 1 - (x^2 / (rep(aggregate(x, by = list(years), 
                                                            function(y) median(y^2))$x,length(unique(years))) + x^2))))
        } else {
          tm    <- lapply(1:5, function(y) t(apply(trackS[,3,y,,], 3:2, max, na.rm = T)))
          surv  <- lapply(tm,  function(y) apply(y, 2, function(x) 1 - (x^2 / (rep(aggregate(x, by = list(years), 
                                                                                             function(y) median(y^2))$x,length(unique(years))) + x^2))))
        }    

        for(dest in 1:dim(trackS)[5]) { ### breeding pops
            
            tmp01 <- trackS[,1:2,dest, surv[[pop]][,dest]>=quantile(surv[[pop]][,dest], probs = 0.75) ,pop]
            sList[[t]][[s]][[pop]][[dest]][[1]] <- (apply(trackS[,3,dest, surv[[pop]][,dest]>=quantile(surv[[pop]][,dest], probs = 0   ) ,pop], 2, max, na.rm = T)*3)/24
            sList[[t]][[s]][[pop]][[dest]][[2]] <- (apply(trackS[,3,dest, surv[[pop]][,dest]>=quantile(surv[[pop]][,dest], probs = 0.75) ,pop], 2, max, na.rm = T)*3)/24
            
            trks <- lapply(1:dim(tmp01)[3], function(x) {
                            st_as_sf(SpatialLines(list(Lines(list(Line(project(tmp01[!is.na(tmp01[,1,x]),1:2,x],
                                proj = proj4string(HexPols)))), ID = 1)),
                                    CRS(proj4string(HexPols)))) })

            inds <- lapply(trks, function(x)  unlist(c(st_intersects(x, st_as_sf(HexPols)))))

            for(i in 1:length(inds)) {
              diag(mat[,,s,t])[inds[[i]]] <- diag(mat[,,s,t])[inds[[i]]] + 1
            }
            
        }
        
          
      
        } else {
        
        if(s==1) {
          tm    <- lapply(1:5, function(y) t(apply(trackS[,3,,,y], 2:3, function(x) ifelse(all(is.na(x)), 0, min(x, na.rm = T)))))
          surv  <- lapply(tm,  function(y) apply(y, 2, function(x) (x^2 / (rep(aggregate(x, by = list(years), 
                                                                       function(y) median(y^2))$x,length(unique(years))) + x^2))))
        } else {
          tm    <- lapply(1:5, function(y) t(apply(trackS[,3,y,,], 3:2, function(x) ifelse(all(is.na(x)), 0, min(x, na.rm = T)))))
          surv  <- lapply(tm,  function(y) apply(y, 2, function(x) (x^2 / (rep(aggregate(x, by = list(years), 
                                                                                         function(y) median(y^2))$x,length(unique(years))) + x^2))))
        }
        
        
        for(dest in 1:dim(trackS)[5]) { ### breeding pops
          
          tmp01 <- trackS[,1:2,dest, surv[[pop]][,dest]>=quantile(surv[[pop]][,dest], probs = 0.75) ,pop]
          sList[[t]][[s]][[pop]][[dest]][[1]] <- apply(trackS[, 3, dest, surv[[pop]][,dest]>=quantile(surv[[pop]][,dest], probs = 0   ) ,pop], 2, 
                                                function(x) ifelse(all(is.na(x)), 0, min(x, na.rm = T)))
          sList[[t]][[s]][[pop]][[dest]][[2]] <- apply(trackS[, 3, dest, surv[[pop]][,dest]>=quantile(surv[[pop]][,dest], probs = 0.75) ,pop], 2, 
                                                function(x) ifelse(all(is.na(x)), 0, min(x, na.rm = T)))
          
          trks <- lapply(1:dim(tmp01)[3], function(x) {
            st_as_sf(SpatialLines(list(Lines(list(Line(project(tmp01[!is.na(tmp01[,1,x]),1:2,x],
                                                               proj = proj4string(HexPols)))), ID = 1)),
                                  CRS(proj4string(HexPols)))) })

          inds <- lapply(trks, function(x)  unlist(c(st_intersects(x, st_as_sf(HexPols)))))

          for(i in 1:length(inds)) {
            diag(mat[,,s,t])[inds[[i]]] <- diag(mat[,,s,t])[inds[[i]]] + 1
          }
          
        }
      }
          
        }
  
  }
}


cls <- rev(wesanderson::wes_palette("Darjeeling1", n = 5, type = "discrete")[c(1,3,2,4,5)])

end0   <- stend[[2]]
start0 <- stend[[1]]

opar <- par(mfrow = c(2,2), mar = c(0,0,0,0)) ###----

test     <- diag(mat[,,2,1])
cols     <- c("transparent", rev(viridis::inferno(max(mat[,,,1], na.rm = T))))
plot(landCirc, border = NA, col = "grey80")
plot(spTransform(HexPols, CRS(proj)), col = ifelse(test==0, "transparent", cols[test+1]), border = NA, add = T)
plot(landCirc, add = T, lwd = 0.6)

points(project(start0, proj = proj)[,1], project(start0, proj = proj)[,2], pch = 21, bg = cls, cex = 2.5, col = "white", lwd = 2)
points(project(end0, proj = proj)[,1], project(end0, proj = proj)[,2], pch = 23, bg = "white", cex = 2)

test     <- diag(mat[,,1,1])
plot(landCirc, border = NA, col = "grey80")
plot(spTransform(HexPols, CRS(proj)), col = ifelse(test==0, "transparent", cols[test+1]), border = NA, add = T)
plot(landCirc, add = T, lwd = 0.6)

points(project(start0, proj = proj)[,1], project(start0, proj = proj)[,2], pch = 21, bg = cls, cex = 2.5, col = "white", lwd = 2)
points(project(end0, proj = proj)[,1], project(end0, proj = proj)[,2], pch = 23, bg = "white", cex = 2)

test     <- diag(mat[,,2,2])
cols     <- c("transparent", rev(viridis::viridis(max(mat[,,,2], na.rm = T))))
plot(landCirc, border = NA, col = "grey80")
plot(spTransform(HexPols, CRS(proj)), col = ifelse(test==0, "transparent", cols[test+1]), border = NA, add = T)
plot(landCirc, add = T, lwd = 0.6)

points(project(start0, proj = proj)[,1], project(start0, proj = proj)[,2], pch = 21, bg = cls, cex = 2.5, col = "white", lwd = 2)
points(project(end0, proj = proj)[,1], project(end0, proj = proj)[,2], pch = 23, bg = "white", cex = 2)

test     <- diag(mat[,,1,2])
plot(landCirc, border = NA, col = "grey80")
plot(spTransform(HexPols, CRS(proj)), col = ifelse(test==0, "transparent", cols[test+1]), border = NA, add = T)
plot(landCirc, add = T, lwd = 0.6)

points(project(start0, proj = proj)[,1], project(start0, proj = proj)[,2], pch = 21, bg = cls, cex = 2.5, col = "white", lwd = 2)
points(project(end0, proj = proj)[,1], project(end0, proj = proj)[,2], pch = 23, bg = "white", cex = 2)

par(opar) 


###----

dsts <- lapply(1:5, function(x) distVincentySphere(stend[[1]][x,], stend[[2]])/1000)

    dates   <- tmS
    years   <- as.numeric(format(dates, "%Y"))

SOut <- list(t1 = list(s1 = list(), s2 = list()), t2 = list(s1 = list(), s2 = list()))
       
prob = 0.75
 
for(t in 1:2) {
 
   for(s in 1:2) {
    
    trackS <- tracks[,,,,,s,t] # [1] 100   3   5 200   5
    
    if(t==1) {
      if(s==1) {
        surv  <- lapply(1:5, function(y) apply(t(apply(trackS[,3,,,y], 2:3, max, na.rm = T)), 2, 
                             function(x) 1 - (x^2 / (rep(aggregate(x, by = list(years), function(y) median(y^2))$x,length(unique(years))) + x^2))))
        
        ftim  <- lapply(1:5, function(x) t(apply(trackS[,3,,,x], 2:3, function(y) (max(y, na.rm =T)*3)/24)))
        ntime <- lapply(1:5, function(x) t(apply(t(apply(trackS[,3,,,x], 2:3, function(y) (max(y, na.rm = T)*3)/24)), 1, function(z) dsts[[x]]/z)))
     
     } else {
       surv  <- lapply(1:5, function(y) apply(t(apply(trackS[,3,y,,], 3:2, max, na.rm = T)), 2, 
                                              function(x) 1 - (x^2 / (rep(aggregate(x, by = list(years), function(y) median(y^2))$x,length(unique(years))) + x^2))))
       
       ftim  <- lapply(1:5, function(x) t(apply(trackS[,3,x,,], 3:2, function(y) (max(y, na.rm =T)*3)/24)))
       ntime <- lapply(1:5, function(x) t(apply(t(apply(trackS[,3,x,,], 3:2, function(y) (max(y, na.rm = T)*3)/24)), 1, function(z) dsts[[x]]/z)))
     }
      
    SOut[[1]][[s]]$ftS <-  cbind(do.call("cbind", lapply(1:5, function(x) quantile(do.call("c", lapply(1:5, 
                                        function(y) ftim[[x]][surv[[x]][,y]>=quantile(surv[[x]][,y], probs = prob),y])), probs = c(0.025, 0.5, 0.975)))),
                                              quantile(unlist(unlist(lapply(1:5, function(x) lapply(1:5, 
                                                function(y) ftim[[x]][surv[[x]][,y]>=quantile(surv[[x]][,y], probs = prob),y])))),  probs = c(0.025, 0.5, 0.975)))

    SOut[[1]][[s]]$ntS1 <-  cbind(do.call("cbind", lapply(1:5, function(x) quantile(do.call("c", lapply(1:5, 
                                        function(y) ntime[[x]][surv[[x]][,y]>=quantile(surv[[x]][,y], probs = prob),y])), probs = c(0.025, 0.5, 0.975)))),
                                             quantile(unlist(unlist(lapply(1:5, function(x) lapply(1:5, 
                                                 function(y) ntime[[x]][surv[[x]][,y]>=quantile(surv[[x]][,y], probs = prob),y])))), probs = c(0.025, 0.5, 0.975)))
    
    SOut[[1]][[s]]$ntS2 <-  cbind(do.call("cbind", lapply(1:5, function(x) quantile(do.call("c", lapply(1:5, 
                               function(y) ntime[[x]][surv[[x]][,y]>=quantile(surv[[x]][,y], probs = prob),y])), probs = c(0.2, 0.5, 0.8)))),
                                 quantile(unlist(unlist(lapply(1:5, function(x) lapply(1:5, 
                                    function(y) ntime[[x]][surv[[x]][,y]>=quantile(surv[[x]][,y], probs = prob),y])))), probs = c(0.2, 0.5, 0.8)))
    
    } else {
      
      if(s==1) {
        surv  <- lapply(1:5, function(y) apply(t(apply(trackS[,3,,,y], 2:3, function(z) suppressWarnings(ifelse(is.infinite(min(z, na.rm = T)), 0, min(z, na.rm = T))))), 2, 
                                               function(x) (x^2 / (rep(aggregate(x, by = list(years), 
                                                                                 function(y) median(y^2))$x,length(unique(years))) + x^2))))
        
        
        
      } else {
        surv  <- lapply(1:5, function(y) apply(t(apply(trackS[,3,y,,], 3:2, function(z) suppressWarnings(ifelse(is.infinite(min(z, na.rm = T)), 0, min(z, na.rm = T))))), 2, 
                                               function(x) (x^2 / (rep(aggregate(x, by = list(years), function(y) median(y^2))$x,length(unique(years))) + x^2))))
        
      }
      
      SOut[[2]][[s]]$survS1 <-  cbind(do.call("cbind", lapply(1:5, function(x) quantile(do.call("c", lapply(1:5, 
                                      function(y) surv[[x]][surv[[x]][,y]>=quantile(surv[[x]][,y], probs = c(0,1)),y])), probs = c(0.025, 0.5, 0.975)))),
                                        quantile(unlist(unlist(lapply(1:5, function(x) lapply(1:5, 
                                          function(y) surv[[x]][surv[[x]][,y]>=quantile(surv[[x]][,y], probs = c(0,1)),y])))), probs = c(0.025, 0.5, 0.975)))
      
      SOut[[2]][[s]]$survS2 <-  cbind(do.call("cbind", lapply(1:5, function(x) quantile(do.call("c", lapply(1:5, 
                                      function(y) surv[[x]][surv[[x]][,y]>=quantile(surv[[x]][,y], probs = prob),y])), probs = c(0.2, 0.5, 0.8)))),
                                        quantile(unlist(unlist(lapply(1:5, function(x) lapply(1:5, 
                                          function(y) surv[[x]][surv[[x]][,y]>=quantile(surv[[x]][,y], probs = prob),y])))), probs = c(0.2, 0.5, 0.8)))
      
      
    }
    
   }
  
}
      
    
opar <- par(mfrow = c(2,1), mar = c(2,2,2,2), oma = c(2,2,2,2), las = 1, bty = "n")

  
  WS1 <- cbind(do.call("cbind", lapply(1:length(sList[[t=1]][[s=1]]), function(x) quantile(as.numeric(do.call("c", lapply(1:length(x), 
                                                                        function(y) dsts[[x]][y]/sList[[t=1]][[s=1]][[x]][[y]][[1]]))), prob = c(0, 0.5, 1)))),
               quantile(do.call("c", lapply(1:length(sList[[t=1]][[s=1]]), function(x) as.numeric(do.call("c", lapply(1:length(x), 
                                                                              function(y) dsts[[x]][y]/sList[[t=1]][[s=1]][[x]][[y]][[1]]))))), prob = c(0, 0.5, 1)))
  
  WS2 <- cbind(do.call("cbind", lapply(1:length(sList[[t=1]][[s=2]]), function(x) quantile(as.numeric(do.call("c", lapply(1:length(x), 
                                                                        function(y) dsts[[x]][y]/sList[[t=1]][[s=2]][[x]][[y]][[1]]))), prob = c(0, 0.5, 1)))),
               quantile(do.call("c", lapply(1:length(sList[[t=1]][[s=2]]), function(x) as.numeric(do.call("c", lapply(1:length(x), 
                                                                              function(y) dsts[[x]][y]/sList[[t=1]][[s=2]][[x]][[y]][[1]]))))), prob = c(0, 0.5, 1)))
  
  plot(NA, xlim = c(0.5, 6.5), ylim = range(WS1), xaxt = "n")
  
  arrows(c(1:6) +0.15, WS1[1,], c(1:6) +0.15, WS1[3,], length = 0, lwd = 2, col = c(cls, "black"))  
  points(c(1:6) +0.15, WS1[2,], pch = 22, lwd = 4, cex = 1.5, bg = "white", col = c(cls, "black"))
  
  arrows(c(1:6) -0.15, WS2[1,], c(1:6) -0.15, WS2[3,], length = 0, lwd = 2, lty = 2, col = c(cls, "black"))
  points(c(1:6) -0.15, WS2[2,], pch = 23, lwd = 4, cex = 1.5, bg = c(cls, "black"), col = c(cls, "black"))
  
  text(x = 1:5, 1250, as.character(round(unlist(lapply(dsts, median)))), xpd = T)
  # axis(1, at = c(1:6), labels = NA)
  
  NP1 <- cbind(do.call("cbind", lapply(1:length(sList[[t=2]][[s=1]]), function(x) quantile(as.numeric(do.call("c", lapply(1:length(x), 
                                                                        function(y) sList[[t=2]][[s=1]][[x]][[y]][[1]]))), prob = c(0, 0.5, 1)))),
               quantile(do.call("c", lapply(1:length(sList[[t=2]][[s=1]]), function(x) as.numeric(do.call("c", lapply(1:length(x), 
                                                                              function(y) sList[[t=2]][[s=1]][[x]][[y]][[1]]))))), prob = c(0, 0.5, 1)))
  
  NP2 <- cbind(do.call("cbind", lapply(1:length(sList[[t=2]][[s=2]]), function(x) quantile(as.numeric(do.call("c", lapply(1:length(x), 
                                                                        function(y) sList[[t=2]][[s=2]][[x]][[y]][[1]]))), prob = c(0, 0.5, 1)))),
               quantile(do.call("c", lapply(1:length(sList[[t=2]][[s=2]]), function(x) as.numeric(do.call("c", lapply(1:length(x), 
                                                                             function(y) sList[[t=2]][[s=2]][[x]][[y]][[1]]))))), prob = c(0, 0.5, 1)))
  
  plot(NA, xlim = c(0.5, 6.5), ylim = c(0,1), xaxt = "n")
  
  arrows(c(1:6) +0.15, NP1[1,], c(1:6) +0.15, NP1[3,], length = 0, lwd = 2, col = c(cls, "black"))  
  points(c(1:6) +0.15, NP1[2,], pch = 22, lwd = 4, cex = 1.5, bg = "white", col = c(cls, "black"))
  
  arrows(c(1:6) -0.15, NP2[1,], c(1:6) -0.15, NP2[3,], length = 0, lwd = 2, lty = 2, col = c(cls, "black"))
  points(c(1:6) -0.15, NP2[2,], pch = 23, lwd = 4, cex = 1.5, bg = c(cls, "black"), col = c(cls, "black"))
  
  text(x = 1:5, 1250, as.character(round(unlist(lapply(dsts, median)))), xpd = T)
  axis(1, at = c(1:6), labels = c("FIN", "SWE", "GER", "CZE", "BUL", "All"))

par(opar)
  
  