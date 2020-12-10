### Indio-European Flyway

source("FlywayAnalysis/setup.R")

########### Tracks                   -----

gTab <- read.csv("SummaryLocations.csv") ## Summary locations of geoloctor data (available upon request from Movebank)
  ids <- unique(as.character(gTab[,1]))
gTab$StartTime <- as.POSIXct(gTab$StartTime, format = "%Y-%d-%m %H:%M", tz = "GMT")
gTab$EndTime <- as.POSIXct(gTab$EndTime, format = "%Y-%d-%m %H:%M", tz = "GMT")
  
flsSummary <- list.files("/RosefinchMigration/Geolocation/Results", pattern = "*Group_fit.RData") ## Fit data of geolocation analysis
  flsID    <- unlist(lapply(strsplit(flsSummary, "_"), function(x) x[[1]]))

SimTab <- data.frame(ID = ids, fls = flsSummary[match(ids, flsID)], dep01 = NA, arr01 = NA, dep02 = NA, arr02 = NA, 
                                                                    depInd01 = NA, arrInd01 = NA, depInd02 = NA, arrInd02 = NA,
                     lon01 = NA, lat01 = NA, lon02 = NA, lat02 = NA, lon03 = NA, lat03 = NA)

for(i in 1:nrow(SimTab)) {
  tmp <- gTab[gTab[,1]==ids[i],]
  load(paste0("/RosefinchMigration/Geolocation/Results/", SimTab$fls[i]))
  dts <- as.numeric(fit$model$time)
  
  if(any(tmp$Type==1)) {
    SimTab$dep01[i]    <- as.POSIXct(as.character(tmp$EndTime)[tmp$Type==1], tz = "GMT")
    SimTab$depInd01[i] <- min(which(dts>=SimTab$dep01[i]))
    
    SimTab$lon01[i]    <- tmp$Lon.50.[tmp$Type==1]
    SimTab$lat01[i]    <- tmp$Lat.50.[tmp$Type==1]
  }
  if(any(tmp$Type==3)) {
    SimTab$arr01[i]    <- as.POSIXct(as.character(tmp$StartTime)[min(which(tmp$Type==3))], tz = "GMT")
    SimTab$arrInd01[i] <- min(which(dts>=SimTab$arr01[i]))
    
    SimTab$lon02[i]    <- tmp$Lon.50.[min(which(tmp$Type==3))]
    SimTab$lat02[i]    <- tmp$Lat.50.[min(which(tmp$Type==3))]
  }
  if(any(tmp$Type==4)) {
    SimTab$dep02[i]    <- as.POSIXct(as.character(tmp$EndTime)[max(which(tmp$Type==3))], tz = "GMT")
    SimTab$depInd02[i] <- min(which(dts>=SimTab$dep02[i]))
    
    SimTab$lon03[i]    <- tmp$Lon.50.[max(which(tmp$Type==3))]
    SimTab$lat03[i]    <- tmp$Lat.50.[max(which(tmp$Type==3))]
  }
  if(any(tmp$Type==5)) {
    SimTab$arr02[i]    <- as.POSIXct(as.character(tmp$StartTime)[max(which(tmp$Type==5))], tz = "GMT")
    SimTab$arrInd02[i] <- min(which(dts>=SimTab$arr02[i]))
  }
}

SimTab$dep01 <- as.POSIXct(SimTab$dep01, origin = "1970-01-01", tz = "GMT")
SimTab$arr01 <- as.POSIXct(SimTab$arr01, origin = "1970-01-01", tz = "GMT")
SimTab$dep02 <- as.POSIXct(SimTab$dep02, origin = "1970-01-01", tz = "GMT")
SimTab$arr02 <- as.POSIXct(SimTab$arr02, origin = "1970-01-01", tz = "GMT")

cnt <- aggregate(gTab$Country, by = list(gTab$ID), function(x) unique(x))
SimTab$Pop   <- cnt[match(SimTab$ID, cnt[,1]),2]

plot(landC, col = "grey80", border = "grey80")
segments(SimTab$lon01, SimTab$lat01, SimTab$lon02, SimTab$lat02, col = "orange")
segments(SimTab$lon03, SimTab$lat03, SimTab$lon01, SimTab$lat01, col = "blue")


########### Comparison               ----

library(SGAT)
library(abind)
library(mgcv)

getGeom <- function(xy, start=0, sep = 1){ 
  
  dx <- c(0,diff(xy[,1])) 
  dy <- c(0,diff(xy[,2])) 
  dseg   <- sqrt(dx^2+dy^2) 
  dtotal <- cumsum(dseg) 
  
  linelength = sum(dseg) 
  
  pos = seq(start, linelength, by = sep) 
  
  whichseg = unlist(lapply(pos, function(x){sum(dtotal<=x)})) 
  
  pos=data.frame(pos=pos,whichseg=whichseg, 
                 x0=xy[whichseg,1], 
                 y0=xy[whichseg,2], 
                 dseg = dseg[whichseg+1], 
                 dtotal = dtotal[whichseg], 
                 x1=xy[whichseg+1,1], 
                 y1=xy[whichseg+1,2] 
  ) 
  
  pos$further =  pos$pos - pos$dtotal 
  pos$f = pos$further/pos$dseg 
  pos$x = pos$x0 + pos$f * (pos$x1-pos$x0) 
  pos$y = pos$y0 + pos$f * (pos$y1-pos$y0) 
  
  pos$theta = atan2(pos$y0-pos$y1,pos$x0-pos$x1) 
  
  return(pos[,c("x","y","x0","y0","x1","y1","theta")]) 
  
} 
transect <- function( tpts, tlen){ 
  
  tpts$thetaT = tpts$theta+pi/2 
  dx = tlen*cos(tpts$thetaT) 
  dy = tlen*sin(tpts$thetaT) 
  return( 
    data.frame(x0 = tpts$x + dx, 
               y0 = tpts$y + dy, 
               x1 = tpts$x - dx, 
               y1 = tpts$y -dy) 
  ) 
  
} 



for(i in 1:nrow(SimTab)) {
  
  for(seas in 1:2) {
  
  load(paste0("/RosefinchMigration/Geolocation/Results/", SimTab$fls[i]))
  
  if(seas==1) {  
  startID <- SimTab$depInd01[i]
  endID   <- SimTab$arrInd01[i]
  } else {
    startID <- SimTab$depInd02[i]
    endID   <- ifelse(!is.na(SimTab$arrInd02[i]), SimTab$arrInd02[i], dim(fit$x[[1]])[1])
  }
  
  if(!is.na(startID)) {
    
    sm <- abind(fit$x[[1]][startID:endID,,1501:2000], along = 3)
    
    out <- do.call("rbind", apply(sm, 3, function(t) {
    
    smSp  <- SpatialLines(list(Lines(list(Line(t)), ID = 1)), proj4string = CRS(proj4string(wrld_simpl)))
    
    if(seas==1) {
    ds  <- distVincentySphere(SimTab[i,c("lon01", "lat01")], SimTab[i,c("lon02", "lat02")])/1000
    gcL <- gcIntermediate(SimTab[i,c("lon01", "lat01")], SimTab[i,c("lon02", "lat02")], n = round(ds/100), addStartEnd = T, sp = F)
    } else {
      ds  <- distVincentySphere(SimTab[i,c("lon03", "lat03")], SimTab[i,c("lon01", "lat01")])/1000
      gcL <- gcIntermediate(SimTab[i,c("lon03", "lat03")], SimTab[i,c("lon01", "lat01")], n = round(ds/100), addStartEnd = T, sp = F)
    }
    tspts   <- getGeom(gcL, 0, 1)   
    tslines <- cbind(transect(tspts, 30), tspts[,1:2])
    tslines[tslines[,4]>90, 4] <- 85  
    gcL     <- SpatialLines(list(Lines(list(Line(gcL)), ID = 1)), proj4string = CRS(proj4string(wrld_simpl)))
    
    simT0  <- trks[trks$ID==SimTab$ID[i] & trks$season==seas,]
    select <- aggregate(simT0$val, by = list(simT0$type, simT0$tm), function(x) unique(x)) 
      simT    <- simT0[which((simT0$type==1 & simT0$tm%in%select[select$Group.1==1,][order(select[select$Group.1==1,3]),2][1:5]) | 
                             (simT0$type==2 & simT0$tm%in%select[select$Group.1==2,][order(select[select$Group.1==2,3], decreasing = T),2][1:5])),]
    simL    <- split(simT, list(simT$tm, simT$type))[lapply(split(simT, list(simT$tm, simT$type)), nrow)>0]
    
       # n <- lapply(simL, function(l) lines(l$lon, l$lat, lwd = 2, col = adjustcolor(ifelse(l$type==1, "blue", "darkgreen"), alpha.f = 0.1)))
    
    data.frame(ID = SimTab$ID[i], Pop = SimTab$Pop[i], do.call("rbind", lapply(split(tslines, f = 1:nrow(tslines)), function(x) {
      
      transL <- SpatialLines(list(Lines(list(Line(matrix(as.numeric(x[,1:4]), ncol = 2, byrow = T))), ID = 1)), proj4string = CRS(proj4string(wrld_simpl)))
      # plot(landC, col = "grey90", border = "grey90")
      # plot(transL, add = T, col = "red")
      intTrk <- gIntersection(transL, smSp)
      intGc  <- gIntersection(transL, gcL)
      
      cbind(x[1,5], do.call("rbind", lapply(simL, function(y) {
        
          spL    <- SpatialLines(list(Lines(list(Line(matrix(unlist(c(y[,c("lon", "lat")])), ncol = 2))), ID = 1)), proj4string = CRS(proj4string(wrld_simpl)))
          intSim <- gIntersection(transL, spL)
          
          if(all(!is.null(c(intSim, intGc, intTrk)))) {
          if(!is.null(intSim) & !is.null(intGc)) {
            dst1 <- apply(cbind(coordinates(intSim)), 1, function(d) distVincentySphere(d, coordinates(intGc))/1000)
            dst1 <- dst1[which.min(dst1)] * ifelse(coordinates(intGc)[,2]>=coordinates(intSim)[which.min(dst1),2], -1, 1)
          } else dst1 <- NA
          if(!is.null(intTrk) & !is.null(intGc)) {
            dst2 <- apply(cbind(coordinates(intTrk)), 1, function(d) distVincentySphere(d, coordinates(intGc))/1000)
            dst2 <- dst2[which.min(dst2)] * ifelse(coordinates(intGc)[,2]>=coordinates(intTrk)[which.min(dst2),2], -1, 1)
          } else dst2 <- NA
          if(!is.null(intTrk) & !is.null(intSim)) {
            dst3 <- apply(cbind(coordinates(intTrk)), 1, function(d) distVincentySphere(d, coordinates(intSim))/1000)
            dst3 <- dst3[which.min(dst3)]
          } else dst3 <- NA
            
            matrix(c(y$type[1], y$season[1], dst1, dst2, dst3), nrow = 1)
          } else matrix(c(y$type[1], y$season[1], NA, NA, NA), nrow = 1)

      })))
      
    }))) 
    
    }))
    
    if(!exists("outDist")) outDist <- out else outDist <- rbind(outDist, out)
    
    
  } ### end if !is.na(startID)    

 
  } # end seas loop

} # end ind loop

# names(outDist) <- c("ID", "Pop", "Lon", "Type", "Season", "Dist_sim", "Dist_track", "Dist_sim_track")  


########### Plot 1                   ----

# cairo_pdf(paste0(getwd(), "/_plot1.pdf"), width = 5, height = 11)
layout(matrix(c(1:10), ncol = 2, byrow = T), widths = c(4,2))
opar <- par(mar = c(3,3,1,1))
for(pop in c("FIN", "SWE", "GER", "CZE", "BUL")) {
  
    
    seqL <- seq(15, 80, length.out = 30)
    
    tab1    <- subset(outDist, Season==1 & Type==1 & Pop==pop)
    cutInd1 <- seqL[cut(tab1$Lon, breaks = seqL, labels = FALSE)]
    aggr1   <- aggregate(tab1$Dist_sim, by = list(cutInd1), FUN = function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = T))
    aggr1   <- merge(data.frame(Group.1 = seqL), as.data.frame(aggr1), all.x = T)

    tab2 <- subset(outDist, Season==1 & Type==2 & Pop==pop)
    cutInd2 <- seqL[cut(tab2$Lon, breaks =seqL, labels = FALSE)]
    aggr2  <- aggregate(tab2$Dist_sim, by = list(cutInd2), FUN = function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = T))
    aggr2   <- merge(data.frame(Group.1 = seqL), as.data.frame(aggr2), all.x = T)
    
    bp <- barplot(t(aggr2[,-1]), plot = FALSE)
    
    tab3 <- subset(outDist, Season==1 & Pop==pop & !duplicated(Dist_track))
    cutInd3 <- seqL[cut(tab3$Lon, breaks = seqL, labels = FALSE)]
    aggr3  <- aggregate(tab3$Dist_track, by = list(cutInd3, tab3$ID), FUN = function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = T))
    aggr3L <- split(data.frame(Lon = as.numeric(aggr3$Group.1), ID = as.numeric(aggr3$Group.2), lower = aggr3$x[,1], median = aggr3$x[,2], upper = aggr3$x[,3]), aggr3$Group.2)
    aggr3L <- lapply(aggr3L, function(x) if(nrow(x)>0) {
      cbind(data.frame(Lon = bp[apply(x, 1, function(m) which(seqL==m[1]))]), x[-1])
      }
      )
    
    bnd  <- cbind(aggregate(tab1$Dist_sim, by = list(cutInd1), FUN = function(x) quantile(x, probs = c(0, 1), na.rm = T)),
                  aggregate(tab2$Dist_sim, by = list(cutInd2), FUN = function(x) quantile(x, probs = c(0, 1), na.rm = T))$x)
    
    aggr4 <- do.call("rbind", lapply(split(tab3, cutInd3), function(x) {
      lns <- seqL[cut(x$Lon, breaks = seqL, labels = FALSE)]
      t(apply(cbind(lns, x$Dist_track), 1, function(y) cbind(y[1], ifelse(y[2]>=bnd[which(bnd$Group.1==y[1]),2][1] & y[2]<=bnd[which(bnd$Group.1==y[1]),2][2] & 
                                                                          y[2]>=bnd[which(bnd$Group.1==y[1]),3]    & y[2]<=bnd[which(bnd$Group.1==y[1]),4], 2, 
                                                                   ifelse( (y[2]>=bnd[which(bnd$Group.1==y[1]),2][1] & y[2]<=bnd[which(bnd$Group.1==y[1]),2][2]) &
                                                                          !(y[2]>=bnd[which(bnd$Group.1==y[1]),3]    & y[2]<=bnd[which(bnd$Group.1==y[1]),4]), 1, ifelse(
                                                                          !(y[2]>=bnd[which(bnd$Group.1==y[1]),2][1] & y[2]<=bnd[which(bnd$Group.1==y[1]),2][2]) &
                                                                           (y[2]>=bnd[which(bnd$Group.1==y[1]),3]    & y[2]<=bnd[which(bnd$Group.1==y[1]),4]), 3, 0))))))
    }))
    aggr5 <- aggregate(aggr4[,2], by = list(aggr4[,1]), FUN = function(x) cbind(sum(x==0), sum(x==1), sum(x==2), sum(x==3))/length(x))
    aggr6 <- merge(data.frame(Group.1 = seqL), as.data.frame(aggr5), all.x = T)
    
    
    bp <- barplot(t(aggr6[,-1]), col = c("white", "cornflowerblue", "brown", "darkolivegreen4"), ylim = c(-0.15, 6), xaxt = "n", yaxt = "n", xlab = "", ylab = "", lwd = 0.25)
    
    # apply(aggr6[aggr6[,1]>50,], 2, mean, na.rm = T)
    # apply(aggr6[aggr6[,1]>50,], 2, sd, na.rm = T)
    
    # apply(aggr6[aggr6[,1]<50,], 2, mean, na.rm = T)
    # apply(aggr6[aggr6[,1]<50,], 2, sd, na.rm = T)
    
    # apply(aggr6, 2, mean, na.rm = T)
    # apply(aggr6, 2, sd, na.rm = T)
    
    par(new = T)
    plot(NA, xlim = range(bp), ylim = c(-3000, 1850), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
    axis(2, at = c(-1500, 0, 1500))
    abline(h = 0, lty = 3)
    axis(1, approx(seqL, bp, xout = seq(20, 80, by = 10))$y, labels = seq(20, 80, by = 10))
    
    
    polygon(c(bp[!is.na(aggr1$x[,1])], rev(bp[!is.na(aggr1$x[,1])])), c(aggr1$x[!is.na(aggr1$x[,1]),1], rev(aggr1$x[!is.na(aggr1$x[,1]),3])), border = NA, col = adjustcolor("cornflowerblue", alpha.f = 0.4))
    polygon(c(bp[!is.na(aggr2$x[,1])], rev(bp[!is.na(aggr2$x[,1])])), c(aggr2$x[!is.na(aggr2$x[,1]),1], rev(aggr2$x[!is.na(aggr2$x[,1]),3])), border = NA, col = adjustcolor("darkolivegreen4", alpha.f = 0.4))
    lapply(aggr3L, function(x) {
      if(!is.null(x)) apply(x[,c(1,3,5)], 1, function(y) segments(y[1], y[2], y[1], y[3], col = adjustcolor("grey20", alpha.f = 0.25)))
      points(x[,1], x[,4], pch = 16, col = "grey40")
    })
    
    # aggr7 <- aggregate(outDist[outDist$Season==1 & !is.na(outDist$Dist_sim_track) & outDist$Pop==pop,"Dist_sim_track"], by = list(outDist[outDist$Season==1 & !is.na(outDist$Dist_sim_track) & outDist$Pop==pop,"Type"]),
                       # FUN = function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = T))
    
    plot(NA, ylim =  c(-3000, 1850), xlim = c(0.5, 3.5), yaxt = "n", xlab = "", ylab = "", xaxt = "n")
      abline(h = 0, lty = 3)
    vioplot(add = TRUE, tab1$Dist_sim[tab1$Dist_sim>=quantile(tab1$Dist_sim, prob = 0.025, na.rm = T) & tab1$Dist_sim<=quantile(tab1$Dist_sim, prob = 0.975, na.rm = T)], 
                        tab2$Dist_sim[tab2$Dist_sim>=quantile(tab2$Dist_sim, prob = 0.025, na.rm = T) & tab2$Dist_sim<=quantile(tab2$Dist_sim, prob = 0.975, na.rm = T)],
                        tab3$Dist_track[tab3$Dist_track>=quantile(tab3$Dist_track, prob = 0.025, na.rm = T) & tab3$Dist_track<=quantile(tab3$Dist_track, prob = 0.975, na.rm = T)], 
            col = c('cornflowerblue', 'darkolivegreen4', 'grey80'))

    # cbind(mean(tab3$Dist_track[tab3$Dist_track>=quantile(tab3$Dist_track, prob = 0.025, na.rm = T) & tab3$Dist_track<=quantile(tab3$Dist_track, prob = 0.975, na.rm = T) & tab3$Lon>50], na.rm = TRUE),
    #       sd(tab3$Dist_track[tab3$Dist_track>=quantile(tab3$Dist_track, prob = 0.025, na.rm = T) & tab3$Dist_track<=quantile(tab3$Dist_track, prob = 0.975, na.rm = T) & tab3$Lon>50], na.rm = TRUE))

    # cbind(mean(tab1$Dist_sim[tab1$Dist_sim>=quantile(tab1$Dist_sim, prob = 0.025, na.rm = T) & tab1$Dist_sim<=quantile(tab1$Dist_sim, prob = 0.975, na.rm = T) & tab1$Lon<50], na.rm = TRUE),
    #         sd(tab1$Dist_sim[tab1$Dist_sim>=quantile(tab1$Dist_sim, prob = 0.025, na.rm = T) & tab1$Dist_sim<=quantile(tab1$Dist_sim, prob = 0.975, na.rm = T) & tab1$Lon<50], na.rm = TRUE))

    # cbind(mean(tab1$Dist_sim[tab1$Dist_sim>=quantile(tab1$Dist_sim, prob = 0.025, na.rm = T) & tab1$Dist_sim<=quantile(tab1$Dist_sim, prob = 0.975, na.rm = T)], na.rm = TRUE),
          # sd(tab1$Dist_sim[tab1$Dist_sim>=quantile(tab1$Dist_sim, prob = 0.025, na.rm = T) & tab1$Dist_sim<=quantile(tab1$Dist_sim, prob = 0.975, na.rm = T)], na.rm = TRUE))
    
    # cbind(mean(tab2$Dist_sim[tab2$Dist_sim>=quantile(tab2$Dist_sim, prob = 0.025, na.rm = T) & tab2$Dist_sim<=quantile(tab2$Dist_sim, prob = 0.975, na.rm = T)], na.rm = TRUE),
          # sd(tab2$Dist_sim[tab2$Dist_sim>=quantile(tab2$Dist_sim, prob = 0.025, na.rm = T) & tab2$Dist_sim<=quantile(tab1$Dist_sim, prob = 0.975, na.rm = T)], na.rm = TRUE))
    
    
    # darkolivegreen4
    # 
    # lines(c(1,1), quantile(tab1$Dist_sim, probs = c(0.0275, 0.975), na.rm = T))
    # points(1, median(tab1$Dist_sim, na.rm = T))
    # 
    # lines(c(2,2), quantile(tab2$Dist_sim, probs = c(0.0275, 0.975), na.rm = T))
    # points(2, median(tab2$Dist_sim, na.rm = T))
    # 
    # lines(c(3,3), quantile(tab3$Dist_track, probs = c(0.0275, 0.975), na.rm = T))
    # points(3, median(tab3$Dist_track, na.rm = T))
    # abline(h = 0, lty = 3)
    
    # segments(c(1,2,3), aggr7$x[,1], c(1,2), aggr7$x[,3], col = c("cornflowerblue", "darkolivegreen4"), lwd = 3)
    # points(c(1,2), aggr7$x[,2], pch = 16, cex = 4, col = c("cornflowerblue", "darkolivegreen4"))
    # axis(4)
    # axis(1, at = c(1,2), labels = NA)
    # mod <- lme4::lmer(outDist[outDist$Season==1 & !is.na(outDist$Dist_sim_track) & outDist$Pop==pop,"Dist_sim_track"] ~ (outDist[outDist$Season==1 & !is.na(outDist$Dist_sim_track) & outDist$Pop==pop,"Type"])
                 # + (1|outDist[outDist$Season==1 & !is.na(outDist$Dist_sim_track) & outDist$Pop==pop,"ID"]))
    # summary(mod)
    # plot(outDist[outDist$Season==1 & !is.na(outDist$Dist_sim_track) & outDist$Pop==pop,"Dist_sim_track"]~as.factor(outDist[outDist$Season==1 & !is.na(outDist$Dist_sim_track) & outDist$Pop==pop,"Type"]),
    #      xlab = "", ylab = "", ylim = c(0, max(outDist$Dist_sim_track, na.rm = T)))
    
}
par(opar)
# dev.off()

# cairo_pdf(paste0(getwd(), "/_plot2.pdf"), width = 5, height = 11)
layout(matrix(c(1:10), ncol = 2, byrow = T), widths = c(4,2))
opar <- par(mar = c(3,3,1,1))
for(pop in c("FIN", "SWE", "GER", "CZE", "BUL")) {
  
  if(pop!="CZE") {
  seqL <- seq(15, 80, length.out = 30)
  # labs <- seqL[-length(seqL)] + diff(seqL)/2
  
  tab1 <- subset(outDist, Season==2 & Type==1 & Pop==pop)
  cutInd1 <- seqL[cut(tab1$Lon, breaks = seqL, labels = FALSE)]
  aggr1   <- aggregate(tab1$Dist_sim, by = list(cutInd1), FUN = function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = T))
  aggr1   <- merge(data.frame(Group.1 = seqL), as.data.frame(aggr1), all.x = T)
  
  tab2 <- subset(outDist, Season==2 & Type==2 & Pop==pop)
  cutInd2 <- seqL[cut(tab2$Lon, breaks =seqL, labels = FALSE)]
  aggr2  <- aggregate(tab2$Dist_sim, by = list(cutInd2), FUN = function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = T))
  aggr2   <- merge(data.frame(Group.1 = seqL), as.data.frame(aggr2), all.x = T)
  
  bp <- barplot(t(aggr2[,-1]), col = c("white", "cornflowerblue", "brown", "darkolivegreen4"), plot = FALSE)
  
  tab3 <- subset(outDist, Season==2 & Pop==pop & !duplicated(Dist_track))
  cutInd3 <- seqL[cut(tab3$Lon, breaks = seqL, labels = FALSE)]
  aggr3  <- aggregate(tab3$Dist_track, by = list(cutInd3, tab3$ID), FUN = function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = T))
  aggr3L <- split(data.frame(Lon = as.numeric(aggr3$Group.1), ID = as.numeric(aggr3$Group.2), lower = aggr3$x[,1], median = aggr3$x[,2], upper = aggr3$x[,3]), aggr3$Group.2)
  aggr3L <- lapply(aggr3L, function(x) if(nrow(x)>0) {
    cbind(data.frame(Lon = bp[apply(x, 1, function(m) which(seqL==m[1]))]), x[-1])
  }
  )
  
  bnd  <- cbind(aggregate(tab1$Dist_sim, by = list(cutInd1), FUN = function(x) quantile(x, probs = c(0, 1), na.rm = T)),
                aggregate(tab2$Dist_sim, by = list(cutInd2), FUN = function(x) quantile(x, probs = c(0, 1), na.rm = T))$x)
  
  aggr4 <- do.call("rbind", lapply(split(tab3, cutInd3), function(x) {
    lns <- seqL[cut(x$Lon, breaks = seqL, labels = FALSE)]
    t(apply(cbind(lns, x$Dist_track), 1, function(y) cbind(y[1], ifelse(y[2]>=bnd[which(bnd$Group.1==y[1]),2][1] & y[2]<=bnd[which(bnd$Group.1==y[1]),2][2] & 
                                                                          y[2]>=bnd[which(bnd$Group.1==y[1]),3]    & y[2]<=bnd[which(bnd$Group.1==y[1]),4], 2, 
                                                                        ifelse( (y[2]>=bnd[which(bnd$Group.1==y[1]),2][1] & y[2]<=bnd[which(bnd$Group.1==y[1]),2][2]) &
                                                                                  !(y[2]>=bnd[which(bnd$Group.1==y[1]),3]    & y[2]<=bnd[which(bnd$Group.1==y[1]),4]), 1, ifelse(
                                                                                    !(y[2]>=bnd[which(bnd$Group.1==y[1]),2][1] & y[2]<=bnd[which(bnd$Group.1==y[1]),2][2]) &
                                                                                      (y[2]>=bnd[which(bnd$Group.1==y[1]),3]    & y[2]<=bnd[which(bnd$Group.1==y[1]),4]), 3, 0))))))
  }))
  aggr5 <- aggregate(aggr4[,2], by = list(aggr4[,1]), FUN = function(x) cbind(sum(x==0), sum(x==1), sum(x==2), sum(x==3))/length(x))
  aggr6 <- merge(data.frame(Group.1 = seqL), as.data.frame(aggr5), all.x = T)
  
  
  apply(aggr6, 2, mean, na.rm = T)
  apply(aggr6, 2, sd, na.rm = T)
  
  
  bp <- barplot(t(aggr6[,-1]), col = c("white", "cornflowerblue", "brown", "darkolivegreen4"), ylim = c(-0.15, 6), xaxt = "n", yaxt = "n", xlab = "", ylab = "", lwd = 0.25)
  
  par(new = T)
  plot(NA, xlim = range(bp), ylim = c(-3000, 1850), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  axis(2, at = c(-1500, 0, 1500))
  abline(h = 0, lty = 3)
  axis(1, approx(seqL, bp, xout = seq(20, 80, by = 10))$y, labels = seq(20, 80, by = 10))
  
  
  polygon(c(bp[!is.na(aggr1$x[,1])], rev(bp[!is.na(aggr1$x[,1])])), c(aggr1$x[!is.na(aggr1$x[,1]),1], rev(aggr1$x[!is.na(aggr1$x[,1]),3])), border = adjustcolor("cornflowerblue", alpha.f = 0.4), col = adjustcolor("cornflowerblue", alpha.f = 0.4))
  polygon(c(bp[!is.na(aggr2$x[,1])], rev(bp[!is.na(aggr2$x[,1])])), c(aggr2$x[!is.na(aggr2$x[,1]),1], rev(aggr2$x[!is.na(aggr2$x[,1]),3])), border = adjustcolor("darkolivegreen4", alpha.f = 0.4), col = adjustcolor("darkolivegreen4", alpha.f = 0.4))
  lapply(aggr3L, function(x) {
    if(!is.null(x)) apply(x[,c(1,3,5)], 1, function(y) segments(y[1], y[2], y[1], y[3], col = adjustcolor("grey20", alpha.f = 0.25)))
    points(x[,1], x[,4], pch = 16, col = "grey40")
  })
  
  # aggr7 <- aggregate(outDist[outDist$Season==2 & !is.na(outDist$Dist_sim_track) & outDist$Pop==pop,"Dist_sim_track"], by = list(outDist[outDist$Season==2 & !is.na(outDist$Dist_sim_track) & outDist$Pop==pop,"Type"]),
  #                    FUN = function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = T))
  # 
  
  plot(NA, ylim =  c(-3000, 1850), xlim = c(0.5, 3.5), yaxt = "n", xlab = "", ylab = "", xaxt = "n")
  abline(h = 0, lty = 3)
  vioplot(add = TRUE, tab1$Dist_sim[tab1$Dist_sim>=quantile(tab1$Dist_sim, prob = 0.025, na.rm = T) & tab1$Dist_sim<=quantile(tab1$Dist_sim, prob = 0.975, na.rm = T)], 
                      tab2$Dist_sim[tab2$Dist_sim>=quantile(tab2$Dist_sim, prob = 0.025, na.rm = T) & tab2$Dist_sim<=quantile(tab2$Dist_sim, prob = 0.975, na.rm = T)],
                      tab3$Dist_track[tab3$Dist_track>=quantile(tab3$Dist_track, prob = 0.025, na.rm = T) & tab3$Dist_track<=quantile(tab3$Dist_track, prob = 0.975, na.rm = T)], 
          col = c('cornflowerblue', 'darkolivegreen4', 'grey80'))
  
  # cbind(mean(tab3$Dist_track[tab3$Dist_track>=quantile(tab3$Dist_track, prob = 0.025, na.rm = T) & tab3$Dist_track<=quantile(tab3$Dist_track, prob = 0.975, na.rm = T)], na.rm = TRUE),
  #         sd(tab3$Dist_track[tab3$Dist_track>=quantile(tab3$Dist_track, prob = 0.025, na.rm = T) & tab3$Dist_track<=quantile(tab3$Dist_track, prob = 0.975, na.rm = T)], na.rm = TRUE))

  # cbind(mean(tab1$Dist_sim[tab1$Dist_sim>=quantile(tab1$Dist_sim, prob = 0.025, na.rm = T) & tab1$Dist_sim<=quantile(tab1$Dist_sim, prob = 0.975, na.rm = T) & tab1$Lon<50], na.rm = TRUE),
  #         sd(tab1$Dist_sim[tab1$Dist_sim>=quantile(tab1$Dist_sim, prob = 0.025, na.rm = T) & tab1$Dist_sim<=quantile(tab1$Dist_sim, prob = 0.975, na.rm = T) & tab1$Lon<50], na.rm = TRUE))
  
  # cbind(mean(tab1$Dist_sim[tab1$Dist_sim>=quantile(tab1$Dist_sim, prob = 0.025, na.rm = T) & tab1$Dist_sim<=quantile(tab1$Dist_sim, prob = 0.975, na.rm = T)], na.rm = TRUE),
  # sd(tab1$Dist_sim[tab1$Dist_sim>=quantile(tab1$Dist_sim, prob = 0.025, na.rm = T) & tab1$Dist_sim<=quantile(tab1$Dist_sim, prob = 0.975, na.rm = T)], na.rm = TRUE))
  
  # cbind(mean(tab2$Dist_sim[tab2$Dist_sim>=quantile(tab2$Dist_sim, prob = 0.025, na.rm = T) & tab2$Dist_sim<=quantile(tab2$Dist_sim, prob = 0.975, na.rm = T)], na.rm = TRUE),
  # sd(tab2$Dist_sim[tab2$Dist_sim>=quantile(tab2$Dist_sim, prob = 0.025, na.rm = T) & tab2$Dist_sim<=quantile(tab1$Dist_sim, prob = 0.975, na.rm = T)], na.rm = TRUE))
  
  
  # plot(NA, ylim = c(0, 3500), xlim = c(0.5, 2.5), yaxt = "n", xlab = "", ylab = "", xaxt = "n")
  # segments(c(1,2), aggr7$x[,1], c(1,2), aggr7$x[,3], col = c("cornflowerblue", "darkolivegreen4"), lwd = 3)
  # points(c(1,2), aggr7$x[,2], pch = 16, cex = 4, col = c("cornflowerblue", "darkolivegreen4"))
  # axis(4)
  # axis(1, at = c(1,2), labels = NA)
  } else {
    plot(NA, xlim = c(1,1), ylim = c(1,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    plot(NA, xlim = c(1,1), ylim = c(1,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  }
  
}
par(opar)
dev.off()


