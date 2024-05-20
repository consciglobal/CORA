library(ctmm)
library(move)
library(amt)
library(dplyr)
library(sp)
source("outliehelper.R")
load("/home/jess/Documents/CSG/CORA/data/hr.RData")
cora1 =as.data.frame(hr[1,2]$data[[1]])
test2 <- move(x=cora1$x_,y=cora1$y_,time=cora1$t_, proj=CRS("+proj=longlat +datum=WGS84")) # proj=CRS("+proj=longlat +ellps=WGS84"
test2$HDOP <- cora1$HDOP
test2$VDOP <- cora1$VDOP
test <- as.telemetry(test2)
outliers <- outlie(as.telemetry(test2))

BESSEL_LIMIT <- 2^16
by='d'
if(class(test)[1]=="list")
{ return( lapply(1:length(test), function(i){OUT <- outlie(test[[i]],plot=plot,by=by); if(plot){ graphics::title(names(test)[i]) }; OUT} ) ) }

UERE <- uere(test)
if(!DOP.LIST$horizontal$VAR %in% names(test)) { uere(test) <- UERE } # adds VAR guesstimate columns if missing (DOF==0)
error <- UERE$UERE[,'horizontal']
names(error) <- rownames(UERE$UERE) # R drops dimnames
error <- ctmm(error=error,axes=c('x','y'))
error <- get.error(test,error,calibrate=TRUE)

# time resolution
DT <- diff(test$t_)
time.res <- time_res(DT)
ZERO <- DT==0
if(any(ZERO)) { DT[ZERO] <- time.res[2] }

Vs <- assign_speeds(test,UERE=error,DT=DT,axes=c('x','y'))
v <- Vs$v.t
VAR.v <- Vs$VAR.t

mu <- median.telemetry(test) #where is this function??
d <- get.telemetry(test,axes=c("x","y"))
mu <- get.telemetry(mu,axes=c("x","y"))
mu <- c(mu)

# detrend median
d <- t(d) - mu
# distances from median
if(length(dim(error))==3)
{ d <- t(d) } else {
  d <- colSums(d^2)
  d <- sqrt(d)
}
D <- distanceMLE(d,error,return.VAR=TRUE)
d <- D[,1]
VAR.d <- D[,2]
rm(D)

if('z' %in% names(test))
{
  error <- UERE$UERE[,'vertical']
  names(error) <- rownames(UERE$UERE) # R drops dimnames
  error <- ctmm(error=error,axes=c('z'))
  error <- get.error(test,error,calibrate=TRUE) # VAR
  
  Vz <- assign_speeds(test,UERE=error,DT=DT,axes='z')
  vz <- Vz$v.t
  VAR.vz <- Vz$VAR.t
  
  dz <- get.telemetry(test,axes=c("z"))
  dz <- dz - stats::median(test$z)
  dz <- abs(dz)
  DZ <- distanceMLE(dz,error,axes='z',return.VAR=TRUE)
  dz <- DZ[,1]
  VAR.dz <- DZ[,2]
  rm(DZ)
}

#if(plot)
#{
  # bounding box
  new.plot(test)
  # convert units on telemetry object to match base plot
  test <- unit.telemetry(test,length=get("x.scale",pos=plot.env))
  
  lwd <- Vs$v.dt
  lwd <- if(diff(range(lwd))) { lwd/max(lwd) } else { 0 }
  col <- grDevices::rgb(0,0,lwd,lwd)
  lwd <- 2*lwd
  lwd_indicator <- c(0,lwd) #JESS ADDED TO FLAG SPEED
  #n <- length(data$t)
  n <- length(test$t)
  graphics::segments(x0=test$x[-n],y0=test$y[-n],x1=test$x[-1],y1=test$y[-1],col=col,lwd=lwd,asp=1)
  
  by <- get(by)
  cex <- if(diff(range(by))) { by/max(by) } else { 0 }
  cex_indicator <- unname(cex) #JESS ADDED TO FLAG
  col <- grDevices::rgb(cex,0,0,cex)
  graphics::points(test$x,test$y,col=col,cex=cex,pch=20)
#}

VAR.v <- nant(VAR.v,0)
VAR.d <- nant(VAR.d,0)
if('z' %in% names(test))
{
  VAR.vz <- nant(VAR.vz,0)
  VAR.dz <- nant(VAR.dz,0)
}

R <- data.frame(t=test$t,distance=d,VAR.distance=VAR.d,speed=v,VAR.speed=VAR.v)
if('z' %in% names(test))
{
  R$vertical.distance <- dz
  R$VAR.vertical.distance <- VAR.dz
  R$vertical.speed <- vz
  R$VAR.vertical.speed <- VAR.vz
}
R <- new.outlie(R)
