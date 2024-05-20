NAMES.CI <- c("low","est","high")

# this is stuff that needs to be run first (and in the right order) for S4 crap to work

methods::setOldClass("UERE")
new.UERE <- methods::setClass("UERE",contains="list",representation=methods::representation(info="list"),
                              prototype=methods::prototype(list(UERE=cbind(0),DOF=cbind(0),AICc=Inf,Zsq=Inf,VAR.Zsq=Inf,N=0),info=list()))
#DOF="matrix",AICc="numeric",Zsq="numeric",VAR.Zsq="numeric",N="numeric"

methods::setOldClass("telemetry")
new.telemetry <- methods::setClass("telemetry",contains="data.frame",representation=methods::representation(info="list",UERE="UERE"),
                                   prototype=methods::prototype(data.frame(),info=list(),UERE=new.UERE()) )

methods::setOldClass("ctmm")
new.ctmm <- methods::setClass("ctmm",contains="list",representation=methods::representation(info="list"),
                              prototype=methods::prototype(list(),info=list()))

methods::setOldClass("UD")
new.UD <- methods::setClass("UD",contains="list",representation=methods::representation(info="list",type="character",variable="character",CTMM="ctmm"),
                            prototype=methods::prototype(list(),info=list(),type=character(),variable=character(),CTMM=new.ctmm()) )

methods::setOldClass("variogram")
new.variogram <- methods::setClass("variogram",representation=methods::representation("data.frame",info="list",UERE="UERE"),
                                   prototype=methods::prototype(data.frame(),info=list(),UERE=new.UERE()) )

methods::setOldClass("outlie")
new.outlie <- methods::setClass("outlie",representation=methods::representation("data.frame"),prototype=methods::prototype(data.frame()))

# R drop is very annoying and yet this doesn't do anything despite the error on det(1)
#setMethod('determinant', signature(x='numeric'), identity)

# existing functions -> S4 generics
# this doesn't work
#methods::setGeneric("SpatialPoints",package="sp",signature=signature("coords",...))
#methods::setGeneric("SpatialPolygonsDataFrame",package="sp",signature="Sr")

# existing funtions -> S3 generics
# this works but is masked if you load sp
#SpatialPoints <- function(object,...) UseMethod("SpatialPoints")
#SpatialPoints.matrix <- function(object,...) sp::SpatialPoints(coords=object,...)
#SpatialPoints.data.frame <- function(object,...) sp::SpatialPoints(coords=object,...)

#SpatialPolygonsDataFrame <- function(object,...) UseMethod("SpatialPolygonsDataFrame")
#SpatialPolygonsDataFrame.SpatialPolygons <- function(object,...) sp::SpatialPolygonsDataFrame(Sr=object,...)

# existing S4 generic functions
methods::setGeneric("projection", getGeneric("projection", package="raster"))
#methods::setGeneric("projection<-", getGeneric("projection<-", package="raster"))
methods::setGeneric("raster", getGeneric("raster", package="raster"))
methods::setGeneric("writeVector", getGeneric("writeVector", package="terra"))


# new S3 generic functions
emulate <- function(object,...) UseMethod("emulate")
AICc <- function(object,...) UseMethod("AICc")
speed <- function(object,...) UseMethod("speed")
speeds <- function(object,...) UseMethod("speeds")
mag <- function(x,...) UseMethod("mag")
#modes <- function(object,...) UseMethod("modes") # CRAN is annoying: do not declare generics for non-exported functions
ridges <- function(object,...) UseMethod("ridges")

# internal S3 generic function
pars <- function(...) { UseMethod("pars") }

# speeds assigned by blame
# dt[1] is recording interval
# dt[2] is minimum time between fixes, which can be smaller than dt[1]
# UERE is misnamed and is actually the per-time error covariance
assign_speeds <- function(data,DT=NULL,UERE=0,method="max",axes=c('x','y'))
{
  method <- match.arg(method,c("max","min"))
  
  if(is.null(DT))
  {
    DT <- diff(data$t)
    dt <- time_res(DT)
  }
  
  # inner speed estimates
  v.dt <- speedMLE(data,UERE=UERE,DT=DT,axes=axes)
  VAR.dt <- v.dt$VAR
  v.dt <- v.dt$X
  if(length(v.dt)==1)
  {
    v <- c(v.dt,v.dt)
    VAR <- c(VAR.dt,VAR.dt)
    return(list(v.t=v,VAR.t=VAR,v.dt=v.dt,VAR.dt=VAR.dt))
  }
  
  N <- length(data$t)
  if(length(UERE)==1)
  {
    # end point contingency estimates
    v1 <- speedMLE(data[c(1,3),],DT=sum(DT[1:2]),UERE=UERE)
    v2 <- speedMLE(data[c(N-2,N),],DT=sum(DT[(N-1):N]),UERE=UERE)
  }
  else if(length(dim(UERE))==3) # pull out correct errors for calculation if fed all errors
  {
    v1 <- speedMLE(data[c(1,3),],DT=sum(DT[1:2]),UERE=UERE[c(1,3),,])
    v2 <- speedMLE(data[c(N-2,N),],DT=sum(DT[(N-2):(N-1)]),UERE=UERE[c(N-2,N),,])
  }
  else # pull out correct errors for calculation if fed all errors (circular)
  {
    v1 <- speedMLE(data[c(1,3),],DT=sum(DT[1:2]),UERE=UERE[c(1,3)])
    v2 <- speedMLE(data[c(N-2,N),],DT=sum(DT[(N-2):(N-1)]),UERE=UERE[c(N-2,N)])
  }
  VAR1 <- v1$VAR; v1 <- v1$X
  VAR2 <- v2$VAR; v2 <- v2$X
  
  # left and right estimates - n estimates, 1 lag apart
  v1 <- c(v1,v.dt)
  v2 <- c(v.dt,v2)
  
  VAR1 <- c(VAR1,VAR.dt)
  VAR2 <- c(VAR.dt,VAR2)
  
  if(method=="max")
  {
    # n-1 estimates, 2 lags apart
    v1 <- v1[-N]
    v2 <- v2[-1]
    
    # which side of lag index looks smaller
    LESS <- v1 < v2
    
    vs <- sort(v.dt,index.return=TRUE,decreasing=TRUE)
    is <- vs$ix
    vs <- vs$x
    
    v <- numeric(length(LESS))
    VAR <- numeric(length(LESS))
    # assign blame for smallest speeds, from greatest to least, last/least taking precedence in R
    LESS <- LESS[is]
    v[is+!LESS] <- vs # v.dt[is]
    VAR[is+!LESS] <- VAR.dt[is]
    
    # assign blame for largest speeds, from least to greatest, last/greatest taking precedence in R
    LESS <- rev(LESS)
    is <- rev(is)
    vs <- rev(vs)
    v[is+LESS] <- vs # v.dt[is]
    VAR[is+LESS] <- VAR.dt[is]
    
    rm(is,vs,LESS)
  }
  else if(method=="min")
  {
    # v <- pmin(v1,v2)
    
    v <- cbind(v1,v2)
    is <- apply(v,1,which.min)
    v <- vapply(1:nrow(v),function(i){v[i,is[i]]},v[,1])
    VAR <- vapply(1:nrow(VAR),function(i){VAR[i,is[i]]},v)
  }
  
  return(list(v.t=v,VAR.t=VAR,v.dt=v.dt,VAR.dt=VAR.dt))
}


# estimate the speed between the two rows of data with error UERE & temporal resolution dt
# dt[1] is recording interval
# dt[2] is minimum time between fixes, which can be smaller than dt[1]
speedMLE <- function(data,DT=NULL,UERE=0,axes=c('x','y'),CTMM=ctmm(error=UERE,axes=axes))
{
  AXES <- length(axes)
  
  ######################
  # 1/diff(t) estimate
  
  if(is.null(DT))
  {
    DT <- diff(data$t)
    dt <- time_res(DT)
  }
  f <- 1/DT
  
  ###################
  if(length(UERE)==1) # 1-time errors
  { error <- get.error(data,ctmm(error=UERE,axes=axes)) }
  else # full error array was passed
  { error <- UERE }
  
  ######################
  # distance estimate
  
  # measured distances
  if(AXES==2)
  {
    dx <- diff(data$x)
    dy <- diff(data$y)
    
    if(length(dim(error))==3) # error ellipse
    { dr <- cbind(dx,dy) } # need full vector
    else # error circle or interval
    { dr <- sqrt(dx^2+dy^2) }
    rm(dx,dy)
  }
  else if(AXES==1)
  {
    dr <- diff(data$z)
    dr <- abs(dr)
  }
  
  # point estimate of distance with error>0
  if(length(UERE)>1 || UERE)
  {
    # 2-time errors
    if(length(dim(error))==3) # error ellipse
    { error <- error[-1,,,drop=FALSE] + error[-nrow(error),,,drop=FALSE] }
    else # error circle or interval
    { error <- error[-1] + error[-length(error)] }
    
    DR <- distanceMLE(dr,error,axes=axes,return.VAR=TRUE)
    VAR <- DR[,2]
    DR <- DR[,1]
  }
  else
  {
    DR <- dr
    VAR <- numeric(length(dr))
  }
  
  RETURN <- data.frame(X=DR*f,VAR=VAR*f^2)
  return(RETURN)
}


####################
distanceMLE <- function(dr,error,axes=c('x','y'),return.VAR=FALSE)
{
  if(length(dim(error))==3) # error ellipse
  {
    dr <- sapply(1:nrow(dr),function(i){ abs_bivar(dr[i,],error[i,,],return.VAR=TRUE) })
    dr <- t(dr) #
  }
  else # error circle or interval
  {
    AXES <- length(axes)
    DR <- dr
    SUB <- dr>0 & error>0
    
    if(any(SUB))
    {
      if(AXES==2) # error circle
      {
        # coefficient in transcendental Bessel equation
        # x I0(x) == y I1(x)
        y <- DR[SUB]^2/error[SUB]
        x <- BesselSolver(y)
        # x = DR*dR/error
        
        # fixed for Inf error
        DR[SUB] <- ifelse(error[SUB]<Inf,error[SUB]/DR[SUB]*x,0)
      }
      else if(AXES==1) # error interval
      {
        error[SUB] <- sqrt(error[SUB]) # now standard deviation
        
        # x = y tanh(x y)
        y <- DR[SUB]/error[SUB]
        x <- TanhSolver(y)
        
        DR[SUB] <- ifelse(error[SUB]<Inf,error[SUB]*x,0)
      } # end error interval
    } # end SUB
    
    VAR <- error/(AXES-(dr^2-DR^2)/error)
    dr <- cbind(DR,VAR)
  } # end error circle or error interval
  
  if(!return.VAR) { dr <- dr[,1] }
  return(dr)
}


# x I0(x) == y I1(x)
BesselSolver <- function(y)
{
  # solution storage
  x <- numeric(length(y))
  
  # critical point, below which all point estimates are zero
  # SUB1 <- (y<=2)
  # if(any(SUB1)) { dr[SUB1] <- 0 }
  # x[SUB1] <- 0 # this is now default
  
  # perturbation of sqrt(EQ) of both sides (from x=0) and solving
  SUB2 <- (2<y) & (y<3.6)
  if(any(SUB2))
  {
    x.SUB <- sqrt(2*y[SUB2])
    x[SUB2] <- 4 * sqrt( (x.SUB-2)/(4-x.SUB) )
  }
  
  # expansion EQ/exp (from x=y) and solving
  SUB3 <- (3.6<=y) & (y<=BESSEL_LIMIT)
  if(any(SUB3))
  {
    y.SUB <- y[SUB3]
    
    BI0 <- besselI(y.SUB,0,expon.scaled=TRUE)
    BI1 <- besselI(y.SUB,1,expon.scaled=TRUE)
    
    x[SUB3] <- 2*y.SUB * ( y.SUB*BI0 - (y.SUB+1)*BI1 ) / ( (2*y.SUB-1)*BI0 - (2*y.SUB+1)*BI1 )
  }
  
  # expansion from x=y, solving, and then rational expansion from y=Inf
  SUB4 <- BESSEL_LIMIT<y
  if(any(SUB4))
  {
    y.SUB <- y[SUB4]
    
    x[SUB4] <- 2*y.SUB*(y.SUB-1)/(2*y.SUB-1)
  }
  
  # iterative to solution by expanding EQ/exp from x=x.old and solving
  SUB <- SUB2 & SUB3
  if(any(SUB))
  {
    x.SUB <- x[SUB]
    y.SUB <- y[SUB]
    
    BI0 <- besselI(x.SUB,0,expon.scaled=TRUE)
    BI1 <- besselI(x.SUB,1,expon.scaled=TRUE)
    
    x[SUB] <- x.SUB * ( x.SUB*(x.SUB+y.SUB)*BI0 - (x.SUB^2+(x.SUB+2)*y.SUB)*BI1 ) / ( x.SUB*(x.SUB+y.SUB-1)*BI0 - (x.SUB^2+(x.SUB+1)*y.SUB)*BI1 )
  }
  # one iteration is good enough to get everything below 1% relative error
  
  # this is now the point estimate, including telemetry error
  # SUB <- !SUB1
  
  return(x)
}


# x == y tanh(x y)
TanhSolver <- function(y)
{
  x <- y
  
  SUB <- y<1
  if(any(SUB)) { x[SUB] <- 0 }
  
  # Taylor expansion from x=0
  SUB <- 1<y & y<=1.1106677241426588
  if(any(SUB))
  {
    x[SUB] <- sqrt((5*y[SUB]-sqrt(120-95*y[SUB]^2))/y[SUB]^3)/2
  }
  
  # Taylor expansion from x=y
  SUB <- 1.1106677241426588<=y & y<Inf
  if(any(SUB))
  {
    TEMP <- 2*y[SUB]^2
    x[SUB] <- nant((sinh(TEMP)-TEMP)/(cosh(TEMP)-TEMP+1),1) * y[SUB]
  }
  
  # One Newton-Raphson iteration and we are within 0.5% error max
  SUB <- 1<y & y<Inf
  if(any(SUB))
  {
    TANH <- tanh(x[SUB]*y[SUB])
    GRAD <- (1-TANH^2) * y[SUB]
    # x + dx == y (TANH + GRAD*dx) + O(dx^2)
    # (1 - y*GRAD) dx == y*TANH - x + O(dx^2)
    dx <- (y[SUB]*TANH-x[SUB])/(1-y[SUB]*GRAD)
    x[SUB] <- x[SUB] + dx
  }
  
  return(x)
}


# estimate temporal resolution from data
# dt[1] = recording interval (truncation)
# dt[2] = fix time
time_res <- function(DT)
{
  ### recording interval ###
  dt <- DT[DT>0]
  # fallback values
  if(length(dt)==0) { return(c(1,1)) }
  if(length(dt)==1) { return(c(dt,min(1,dt/2))) }
  
  # assume decimal truncation
  M <- -log10(min(dt))
  M <- ceiling(M)
  M <- max(0,M)
  M <- 10^M
  
  dt <- round(M*dt) # shift decimal place to integers
  dt <- gcd.vec(dt)/M # shift back
  
  ### fix time ###
  # longest repeating 1
  DT <- DT==0
  if(any(DT))
  {
    for(i in 2:length(DT)) { if(DT[i]) { DT[i] <- DT[i] + DT[i-1] } }
    DT <- max(DT)
    DT <- dt/(1+DT)
  }
  else
  { DT <- 0 }
  
  return(c(dt,DT))
}


# greatest common divisor of an array
gcd.vec <- function(vec)
{
  vec <- sort(vec,method='quick')
  
  GCD <- gcd(vec[1],vec[2])
  for(i in vec[-(1:2)]) { GCD <- gcd(i,GCD) }
  # i > GCD because we sorted vec first
  
  return(GCD)
}

# greatest common divisor of 2 numbers
# fastest if x > y, I think...
gcd <- function(x,y)
{
  r <- x%%y;
  return(ifelse(r, gcd(y, r), y))
}

# global variables for dop/uere/error functions (special axes)
DOP.LIST <- list(
  unknown=list(axes=NA,geo=NA,DOP=NA,VAR=NA,COV=NA,COV.geo=NA,units=NA) ,
  horizontal=list(axes=c("x","y"),geo=c("longitude","latitude"),DOP="HDOP",VAR="VAR.xy",COV=c("COV.x.x","COV.x.y","COV.y.y"),COV.geo=c("COV.major","COV.minor","COV.angle"),units="distance"),
  vertical=list(axes="z",geo="z",DOP="VDOP",VAR="VAR.z",COV=NA,COV.geo=NA,units="distance"),
  speed=list(axes=c("vx","vy"),geo=c("speed","heading"),DOP="SDOP",VAR="VAR.v",COV=c("COV.vx.vx","COV.vx.vy","COV.vy.vy"),COV.geo=NA,units="speed"),
  frequency=list(axes='f',geo=NA,DOP=NA,VAR=NA,COV=NA,COV.geo=NA,units='frequency'),
  mass=list(axes='m',geo=NA,DOP=NA,VAR=NA,COV=NA,COV.geo=NA,units='mass')
)


# are the data calibrated
is.calibrated <- function(data,type="horizontal")
{
  if(class(data)[1]=="list") { return( mean( sapply(data,is.calibrated) ) ) }
  
  DOF <- attr(data,"UERE")$DOF # RMS UERE matrix
  
  # classes in data
  CLASS <- levels(data$class)
  if(!length(CLASS))
  {
    CLASS <- "all"
    if(is.null(dimnames(DOF))) { dimnames(DOF) <- list(CLASS,type) }
  }
  
  # UERE of classes in data only
  DOF <- DOF[CLASS,type]
  DOF <- !is.na(DOF) & DOF>0
  DOF <- mean(DOF) # fraction calibrated
  
  return(DOF)
}

# match DOP type for DOP.LIST by axes argument
DOP.match <- function(axes)
{
  DOP.LIST <- DOP.LIST[-1] # skip unknown case
  NAMES <- names(DOP.LIST)
  for(i in 1:length(DOP.LIST)) { if(all(axes==DOP.LIST[[i]]$axes)) { return(NAMES[i]) } }
  # match was not found
  # warning("axes=",paste(axes,collapse=",")," not of known DOP type.")
  return("unknown")
}


# what DOP types are in the data
get.dop.types <- function(data)
{
  # types of data present
  TYPES <- names(DOP.LIST[-1])
  
  data <- listify(data)
  # all kinds of data
  NAMES <- sapply(data,names)
  NAMES <- unique(unlist(NAMES))
  
  IN <- sapply(TYPES,function(type){all(DOP.LIST[[type]]$axes %in% NAMES) || all(DOP.LIST[[type]]$geo %in% NAMES)})
  TYPES <- TYPES[IN]
  
  return(TYPES)
}


# return the RMS UERE object from set data
uere <- function(data)
{
  if(class(data)[1]=="list" && length(data)==1)
  { data <- data[[1]] }
  
  if(class(data)[1] %in% c("telemetry","variogram"))
  {
    UERE <- attr(data,"UERE")
    return(UERE)
  }
  
  UERE <- lapply(data,function(d){attr(d,"UERE")})
  SAME <- sapply(UERE[-1],function(U){identical(U,UERE[[1]])})
  if(all(SAME)) { return( UERE[[1]] ) }
  else { return(UERE) }
}


# calculate UERE value from calibration data
# look for every DOP without a VAR and calculate UERE - override option for all
uere.fit <- function(data,precision=1/2)
{
  data <- listify(data)
  
  TYPES <- get.dop.types(data) # types of DOP data present
  
  # loop over axes/types
  LIST <- lapply(TYPES,function(type){uere.type(data,precision=precision,trace=trace,type=type)})
  
  UERE <- lapply(LIST,function(L){ L$UERE })
  names(UERE) <- TYPES
  
  DOF <- lapply(LIST,function(L){ L$DOF })
  names(DOF) <- TYPES
  
  AICc <- sapply(LIST,function(L){ L$AICc })
  names(AICc) <- TYPES
  
  Zsq <- sapply(LIST,function(L){ L$Zsq })
  names(Zsq) <- TYPES
  
  VAR.Zsq <- sapply(LIST,function(L){ L$VAR.Zsq })
  names(VAR.Zsq) <- TYPES
  
  N <- sapply(LIST,function(L){ L$N })
  names(N) <- TYPES
  
  CLASS <- names(UERE[[1]])
  if(!length(CLASS)) { CLASS <- "all" }
  
  # RRRRRRRRRR WHY DOES R DROP DIMENSIONS & DIMENSION NAMES ......
  UERE <- matrix( simplify2array(UERE) , length(CLASS) , length(TYPES) )
  dimnames(UERE) <- list(CLASS,TYPES)
  DOF <- matrix( simplify2array(DOF) , length(CLASS) , length(TYPES) )
  dimnames(DOF) <- list(CLASS,TYPES)
  
  UERE <- list(UERE=UERE,DOF=DOF,AICc=AICc,Zsq=Zsq,VAR.Zsq=VAR.Zsq,N=N)
  UERE <- new.UERE(UERE)
  
  return(UERE)
}


# uere values for one set of axes
uere.type <- function(data,trace=FALSE,type='horizontal',precision=1/2,...)
{
  TOL <- .Machine$double.eps^precision
  
  axes <- DOP.LIST[[type]]$axes # axes of data type
  DOP <- DOP.LIST[[type]]$DOP # DOP type
  
  for(i in 1:length(data))
  {
    # toss out Inf DOP values (used for missing data types)
    if(DOP %in% names(data[[i]]))
    {
      IN <- (data[[i]][[DOP]] < Inf)
      data[[i]] <- data[[i]][IN,]
    }
    else # make sure data has some DOP value
    { data[[i]][[DOP]] <- 1 }
    
    # null class cannot be NULL
    if("class" %nin% names(data[[i]]))
    {
      data[[i]]$class <- "all"
      data[[i]]$class <- as.factor(data[[i]]$class)
    }
  }
  
  # all location classes
  CLASS <- lapply(data,function(D){levels(D$class)})
  CLASS <- unique(unlist(CLASS))
  
  # null UERE structure given location classes
  UERE <- rep(NA_real_,length(CLASS))
  names(UERE) <- CLASS
  # don't overwrite ARGOS when calibrating GPS
  ARGOS <- grepl('argos',tolower(CLASS))
  if(any(ARGOS) && type=='horizontal') { UERE[ARGOS] <- 1 }
  
  z <- lapply(1:length(data),function(i){get.telemetry(data[[i]],axes)})
  
  # make sure axes present in dataset (e.g., only 2D fixes)
  if(type=='speed') { IN <- sapply(z,length)>0 } else { IN <- sapply(z,length)>length(axes) }
  data <- data[IN]
  z <- z[IN]
  
  # don't have enough data to estimate any UERE
  if(!length(data))
  {
    UERE <- list(UERE=UERE)
    UERE$DOF <- numeric(length(UERE)) # sampling distribution
    UERE$AICc <- NA_real_
    UERE$Zsq <- NA_real_
    UERE$VAR.Zsq <- NA_real_
    UERE$N <- NA_real_
    
    return(UERE)
  }
  
  EST <- is.na(UERE)
  names(EST) <- CLASS
  
  # DOP values # DOP values ensured above
  DOP <- lapply(data,function(D){ D[[DOP]] })
  
  # weights
  w <- lapply(DOP,function(D){length(axes)/D^2})
  
  # class indicators
  Ci <- lapply(data,function(D){get.class(D,CLASS)}) # (animal;time)
  Cim <- lapply(data,function(D){get.class.mat(D,CLASS)}) # (animal;time,class)
  
  # ML DOF
  DOF.ML <- vapply(Cim,colSums,UERE) # (class,animals)
  dim(DOF.ML) <- c(length(UERE),length(data))
  DOF.ML <- rowSums(DOF.ML) # (class)
  names(DOF.ML) <- CLASS
  
  # initial guess of UEREs
  UERE[EST] <- 10
  
  # check for missing levels/classes
  BAD <- (DOF.ML<=0)
  if(any(BAD))
  {
    # missing UEREs should not impact calculation
    UERE[BAD] <- Inf
    EST[BAD] <- FALSE
  }
  
  # iterative fit
  DEBUG <- FALSE
  BREAK <- FALSE
  repeat
  {
    if(type=="speed") # known mean
    {
      # ML/REML degrees of freedom lost
      Kc <- numeric(length(CLASS))
      
      # known means
      mu <- array(0,c(length(data),length(axes)))
    }
    else # unknown mean
    {
      # precision weights
      Pck <- vapply(1:length(data),function(i){c(w[[i]] %*% Cim[[i]]) * (1/UERE^2)},UERE) # (class,animals)
      dim(Pck) <- c(length(UERE),length(data))
      Pc <- rowSums(Pck) # (class)
      Pk <- colSums(Pck) # (animals)
      
      # REML degrees of freedom lost
      Kc <- colSums(t(Pck)/Pk) # (class)
      
      # means
      CdUERE2 <- lapply(1:length(data),function(i){ 1/UERE[Ci[[i]]]^2 }) # (animals;time,class)
      mu <- vapply(1:length(data),function(i){t(w[[i]] * z[[i]]) %*% CdUERE2[[i]]},z[[1]][1,]) # (axes,animals)
      dim(mu) <- c(length(axes),length(data))
      mu <- t(mu)/Pk # (animals,axes)
      
      if(DEBUG)
      {
        loglike <- -(1/2)*vapply(1:length(data),function(i){
          LL <- length(axes) * sum( log(UERE[Ci[[i]]]^2/w[[i]]) )
          LL <- LL + sum( (t(z[[i]])-mu[i,])^2 %*% (w[[i]] * CdUERE2[[i]]) )
          LL <- LL + length(axes)*log( w[[i]] %*% CdUERE2[[i]] ) # REML
          return(LL) },numeric(1))
        loglike <- sum(loglike)
        message(type," log(Like) = ",loglike)
      }
    }
    DOF <- DOF.ML - Kc # (class)
    
    if(BREAK) { break }
    
    # updated UERE esitmates
    UERE2 <- vapply(1:length(data),function(i){colSums((t(z[[i]])-mu[i,])^2) %*% (w[[i]] * Cim[[i]])},UERE) # (class,animals)
    dim(UERE2) <- c(length(UERE),length(data))
    UERE2 <- rowSums(UERE2) / (length(axes)*DOF)
    UERE2 <- sqrt(UERE2)
    UERE2 <- UERE2[EST]
    
    ERROR <- abs(UERE2-UERE[EST])/max(UERE[EST],UERE2)
    ERROR <- nant(ERROR,0)
    ERROR <- max(ERROR)
    
    UERE[EST] <- UERE2
    
    if(ERROR<TOL) { BREAK <- TRUE }
  }
  
  ### AICc values ###
  if(any(BAD)) { UERE[BAD] <- NA } # missing UEREs should not impact calculation
  AICc <- vapply(1:length(data),function(i){ sum( log( (2*pi)*(UERE[Ci[[i]]]^2 / w[[i]]) ), na.rm=TRUE) },numeric(1)) # (animals)
  AICc <- length(axes) * sum(AICc)
  if(type=="speed")
  { dof <- DOF.ML }
  else
  { dof <- DOF.ML - colSums((t(Pck)/Pk)^2) } # (class)
  # missing UEREs should not impact calculation
  if(any(BAD))
  {
    Kc[BAD] <- 0
    dof[BAD] <- 0
  }
  AICc <- AICc + length(axes)^2*sum( ((DOF.ML+Kc)*dof/(length(axes)*dof-2))[EST] )
  if(any(length(axes)*dof[EST]<=2)) # does not have a mean value
  {
    AICc <- Inf
    warning("Sample size too small for AICc to exist.")
  }
  
  ### Z_reduced^2 ###
  if(any(BAD)) { UERE[BAD] <- Inf } # missing UEREs should not impact calculation
  Pi <- lapply(1:length(data),function(k){c( w[[k]] / UERE[Ci[[k]]]^2 )}) # (animal;time)
  u2 <- lapply(1:length(data),function(k){c( colMeans((t(z[[k]])-mu[k,])^2) * Pi[[k]] )}) # (animal;time)
  if(type=="speed") # known mean
  {
    alpha <- lapply(1:length(data),function(k){ array(1,length(Pi[[k]])) }) # (animal;time)
    beta <- lapply(1:length(data),function(k){ ( DOF.ML[Ci[[k]]] - 1 ) / DOF[Ci[[k]]] }) # (animal;time)
    dofi <- lapply(1:length(data),function(k){ DOF.ML[Ci[[k]]] - 1 }) # (animal;time)
  }
  else
  {
    alpha <- lapply(1:length(data),function(k){ R <- Pck[Ci[[k]],k] ; ( R - Pi[[k]] )/R }) # (animal;time)
    gamma <- lapply(1:length(data),function(k){ t(Pck[Ci[[k]],,drop=FALSE] - Pi[[k]])/outer(Pk,Pi[[k]],'-') }) # (animal;animal',time)
    beta <- lapply(1:length(data),function(k){ ( DOF.ML[Ci[[k]]] - 1 - colSums(gamma[[k]]) ) / DOF[Ci[[k]]] }) # (animal;time)
    dofi <- lapply(1:length(data),function(k){ DOF.ML[Ci[[k]]] - 1 - colSums(gamma[[k]]^2) }) # (animal;time)
  }
  t2 <- lapply(1:length(data),function(k){ clamp( beta[[k]] * u2[[k]] / ( alpha[[k]]^2 - alpha[[k]]*u2[[k]]/DOF[Ci[[k]]] ),0,Inf) }) # (animal,time)
  Z <- log(unlist(t2))/2
  # M <- -clamp(dofi-1,0,Inf)/dofi/(2*length(axes)) # asymptotic only
  # VAR <- (dofi+1)/dofi/(2*length(axes)) # asymptotic only
  d1 <- length(axes)
  d2 <- unlist(dofi)
  if(length(CLASS)==1) { d2 <- d2[1] } # minimize computation if possible
  d2 <- length(axes) * clamp(d2,0,Inf) # finish and fix for tiny samples in a class
  M2 <- function(d2)
  {
    if(d2==0) { return(Inf) }
    
    d1 <- d1/2
    d2 <- d2/2
    
    L1 <- -(log(d1)-digamma(d1))/2
    L2 <- -(log(d2)-digamma(d2))/2
    
    Q1 <- L1^2 + trigamma(d1)/4
    Q2 <- L2^2 + trigamma(d2)/4
    
    R <- Q1 + Q2 - 2*L1*L2
    
    return(R)
  }
  V2 <- sapply(d2,M2)
  Z2 <- (Z^2/V2)
  
  # fix missing UEREs
  if(any(BAD)) { UERE[BAD] <- NA }
  if(any(ARGOS) && type=='horizontal') { dof[ARGOS] <- Inf }
  
  UERE <- list(UERE=UERE)
  UERE$DOF <- dof # sampling distribution
  UERE$AICc <- AICc
  UERE$Zsq <- mtmean(Z2) # minimally trimmed mean in case of weird speed=0 estimates
  UERE$VAR.Zsq <- mtmean(Z2^2) - UERE$Zsq^2
  UERE$N <- sum(Z2<Inf)
  
  return(UERE)
}


# summarize uere object
summary.UERE <- function(object,level=0.95,...)
{
  if(class(object)[1]=="list") { return(summary.UERE.list(object,level=level,...)) }
  
  TYPE <- colnames(object$UERE)
  N <- sapply(TYPE,function(type){length(DOP.LIST[[type]]$axes)}) # (type)
  DOF <- object$DOF # (class,type)
  DOF <- t(t(DOF)*N)
  
  UERE <- sapply(1:length(object$UERE),function(i){chisq.ci(object$UERE[i]^2,DOF=DOF[i],level=level)}) #(3CI,class*type)
  UERE <- sqrt(UERE)
  dim(UERE) <- c(3,dim(object$UERE)) # (3CI,class,type)
  dimnames(UERE)[[1]] <- NAMES.CI
  dimnames(UERE)[2:3] <- dimnames(object$UERE)
  UERE <- aperm(UERE,c(2,1,3)) # (class,3CI,type)
  return(UERE)
}


# summarize list of uere objects
summary.UERE.list <- function(object,level=0.95,drop=TRUE,CI=FALSE,...)
{
  NAMES <- names(object)
  if(is.null(NAMES)) { NAMES <- 1:length(object) }
  
  # determine all types of data
  # better be consistent across joint models
  # might not be consistent within joint models' individual models
  TYPES <- NULL
  for(i in 1:length(object))
  {
    if(class(object[[i]])[1]=="UERE") # these should be all the same
    {
      TYPES <- c(TYPES,names(object[[i]]$AICc))
      TYPES <- unique(TYPES)
    }
    else if(class(object[[i]])[1]=="list") # these can differ within list
    {
      for(j in 1:length(object[[i]]))
      {
        TYPES <- c(TYPES,names(object[[i]][[j]]$AICc))
        TYPES <- unique(TYPES)
      }
    }
  }
  
  # aggregate lists of individual models to compare with joint models
  for(i in 1:length(object))
  {
    if(class(object[[i]])[1]=="list")
    {
      AICc <- rep(0,length(TYPES))
      names(AICc) <- TYPES
      N <- AICc
      Zsq <- AICc
      VAR.Zsq <- AICc
      
      for(j in 1:length(object[[i]]))
      {
        U <- object[[i]][[j]]
        IN <- TYPES %in% names(U$AICc)
        if(length(IN))
        {
          AICc[IN] <- AICc[IN] + U$AICc[IN]
          N[IN] <- N[IN] + U$N[IN]
          Zsq[IN] <- Zsq[IN] + U$N[IN] * U$Zsq[IN]
          VAR.Zsq[IN] <- VAR.Zsq[IN] + U$N[IN] * (U$VAR.Zsq[IN] + U$Zsq[IN]^2)
        }
      }
      
      Zsq <- Zsq/N
      VAR.Zsq <- VAR.Zsq/N - Zsq^2
      
      object[[i]]$AICc <- AICc
      object[[i]]$N <- N
      object[[i]]$Zsq <- Zsq
      object[[i]]$VAR.Zsq <- VAR.Zsq
    }
  }
  
  AIC <- sapply(object,function(U) { U$AICc } )
  N <- sapply(object,function(U) { U$N } )
  Zsq <- sapply(object,function(U) { U$Zsq } )
  VAR.Zsq <- sapply(object,function(U) { U$VAR.Zsq } )
  
  DIM <- dim(AIC)
  if(!length(DIM))
  {
    DIM <- c(1,length(AIC))
    dim(AIC) <- DIM
    dim(Zsq) <- DIM
    dim(VAR.Zsq) <- DIM
    dim(N) <- DIM
  }
  
  TAB <- list()
  for(i in 1:DIM[1])
  {
    AIC[i,] <- AIC[i,] - min(AIC[i,])
    
    if(CI)
    {
      # is N correct here for 1D and 2D ???
      TAB[[i]] <- sapply(1:DIM[2],function(j){ chisq.ci(Zsq[i,j],VAR=VAR.Zsq[i,j]/N[i,j],level=level) }) # (3CIS,models)
      TAB[[i]] <- cbind(AIC[i,],t(TAB[[i]]))
      colnames(TAB[[i]]) <- c("\u0394AICc","(       ","Z[red]\u00B2","       )")
    }
    else
    {
      TAB[[i]] <- cbind(AIC[i,],Zsq[i,])
      colnames(TAB[[i]]) <- c("\u0394AICc","Z[red]\u00B2")
    }
    # Encoding(colnames(TAB)) <- "UTF-8"
    
    
    rownames(TAB[[i]]) <- NAMES
    
    IND <- order(AIC[i,])
    TAB[[i]] <- TAB[[i]][IND,]
  }
  names(TAB) <- rownames(AIC)
  
  if(length(TAB)==1 && drop) { return(TAB[[1]]) }
  else { return(TAB) }
}


# calculate residuals of calibration data
# add axes argument for other uses?
residuals.calibration <- function(data,TYPES=get.dop.types(data),...)
{
  # enforce list structure
  data <- listify(data)
  
  for(TYPE in TYPES)
  {
    axes <- DOP.LIST[[TYPE]]$axes
    
    # calculate mean and residuals
    for(i in 1:length(data))
    {
      n <- nrow(data[[i]])
      
      if(is.calibrated(data[[i]])<1) { uere(data[[i]]) <- uere(data[[i]]) } # force calibration for plotting
      UERE <- uere(data[[i]])
      ERROR <- UERE$UERE[,'horizontal']
      names(ERROR) <- rownames(UERE$UERE) # R drops dimnames
      ERROR <- ctmm(error=ERROR,axes=axes)
      ERROR <- get.error(data[[i]],ERROR,calibrate=TRUE)
      ELLIPSE <- attr(ERROR,"ellipse")
      
      # locations
      z <- get.telemetry(data[[i]],axes)
      
      if(!ELLIPSE)
      {
        # now these are the weights
        w <- 1/ERROR
        W <- sum(w)
        
        # stationary mean
        mu <- c(w %*% z)/W
        
        # detrend the mean for error/residuals
        z <- t(t(z) - mu)
        # x <- rbind(x,dz)
        
        # UERE-1 standardize residuals
        z <- sqrt(w) * z
      }
      else
      {
        w <- vapply(1:n,function(j){PDsolve(ERROR[j,,])},diag(2)) # [x,x,t]
        dim(w) <- c(2,2,n)
        w <- aperm(w,c(3,1,2)) # [t,x,x]
        W <- apply(w,2:3,sum) # [x,y]
        
        # stationary mean
        mu <- vapply(1:n,function(j){w[j,,] %*% z[j,]},1:2) # [x,t]
        mu <- apply(mu,1,sum)
        mu <- c(PDsolve(W) %*% mu)
        
        # detrend the mean for error/residuals
        z <- t(t(z) - mu)
        
        # standardize/studentize
        w <- vapply(1:n,function(j){sqrtm(ERROR[j,,])},diag(2)) # [x,x,t]
        dim(w) <- c(2,2,n)
        w <- aperm(w,c(3,1,2)) # [t,x,x]
        
        z <- vapply(1:n,function(j){w[j,,] %*% z[j,]},1:2) # [x,t]
        z <- t(z) # [t,x]
      }
      
      # store back in the object
      data[[i]][,axes] <- z
      
      uere(data) <- NULL
    } # end data loop
  } # end type loop
  
  return(data)
}


## prepare error array, also return a ellipse necessity #
# circle : force variance   scalar output
# DIM : force covariance matrix output with dim [DIM,DIM]
get.error <- function(data,CTMM,circle=FALSE,DIM=FALSE,calibrate=TRUE)
{
  n <- nrow(data)
  axes <- CTMM$axes
  COLS <- names(data)
  ELLIPSE <- FALSE
  error <- rep(0,n)
  
  if(any(CTMM$error>0))
  {
    TYPE <- DOP.match(axes)
    UERE.DOF <- attr(data,"UERE")$DOF[,TYPE]
    names(UERE.DOF) <- rownames(attr(data,"UERE")$DOF)
    
    # expand to what classes are in the UERE object
    ERROR <- rep(FALSE,length(UERE.DOF))
    names(ERROR) <- names(UERE.DOF)
    ERROR[names(CTMM$error)] <- CTMM$error
    
    UERE.FIT <- ERROR & !is.na(UERE.DOF) & UERE.DOF<Inf # will we be fitting any error parameters?
    UERE.FIX <- ERROR & (is.na(UERE.DOF) | UERE.DOF==Inf)
    UERE.PAR <- names(UERE.FIT)[UERE.FIT>0] # names of fitted UERE parameters
    
    # make sure that non-fitted parameters are logicals
    if(any(!UERE.FIT)) { ERROR[!UERE.FIT] <- as.logical(ERROR[!UERE.FIT]) }
    
    # reduce to what classes are actually in the data - also fixes bad sorting
    if("class" %in% names(data))
    {
      LEVELS <- levels(data$class)
      ERROR <- ERROR[LEVELS]
      UERE.DOF <- UERE.DOF[LEVELS]
      UERE.FIT <- UERE.FIT[LEVELS]
      UERE.FIX <- UERE.FIX[LEVELS]
      UERE.PAR <- UERE.PAR[UERE.PAR %in% LEVELS]
    }
    
    CLASS <- get.class.mat(data)
    FIT <- as.logical(c(CLASS %*% UERE.FIT)) # times where errors will be fitted (possibly with prior)
    # FIX <- as.logical(c(CLASS %*% UERE.FIX)) # times where errors are fixed
    
    # DOP.LIST is global variable from uere.R
    TYPE <- DOP.LIST[[TYPE]]
    AXES <- TYPE$axes
    COV <- TYPE$COV
    VAR <- TYPE$VAR
    DOP <- TYPE$DOP
    
    ELLIPSE <- all(COV %in% COLS)
    CIRCLE <- VAR %in% COLS
    
    # location classes with fixed (known) parameters
    if(ELLIPSE) # at least some (fixed) error ellipses - ARGOS / VHF - this should also include other fixed errors (promoted to ellipses)
    {
      ellipse <- get.telemetry(data,COV[c(1,2,2,3)]) # pull matrix elements
      dim(ellipse) <- c(nrow(ellipse),2,2) # array of matrices
      
      # reduce to VAR
      if(circle)
      {
        error <- (ellipse[,1,1]+ellipse[,2,2])/2
        ELLIPSE <- FALSE
      }
    }
    else if(CIRCLE) # some fixed error circles (no ellipses)
    {
      error <- data[[VAR]] # FITs will overwrite where appropriate
    }
    
    # location classes with fitted error parameters
    if(any(UERE.FIT))
    {
      if(calibrate) # apply RMS UERE to DOP values
      { CLASS <- c( CLASS %*% ERROR ) }
      else # just calculate relative variance
      { CLASS <- c( CLASS %*% as.logical(ERROR) ) }
      
      if(DOP %in% COLS) # apply DOP values
      { CLASS <- CLASS * data[[DOP]] }
      
      CLASS <- CLASS^2/length(axes) # VAR=(RMS[UERE]*HDOP)^2/2 in 2D
      error[FIT] <- CLASS[FIT]
      # FIX was already taken care of in if(ELLIPSE) & if(CIRCLE)
    }
    
    # promote circular errors to elliptical errors
    if(ELLIPSE)
    {
      # copy over error circles
      if(any(UERE.FIT))
      {
        ellipse[FIT,1,1] <- ellipse[FIT,2,2] <- error[FIT]
        ellipse[FIT,1,2] <- ellipse[FIT,2,1] <- 0
      }
      
      error <- ellipse
    }
  } # END errors
  
  # upgrade variance scalar to covariance matrix
  if(!ELLIPSE && !circle && DIM) { error <- outer(error,diag(DIM)) } # [n,d,d] }
  
  attr(error,"ellipse") <- ELLIPSE
  return(error)
}


# try to assign one UERE to one data
# some NA value means Inf VAR; all NA values means delete VAR column
try.assign.uere <- function(data,UERE,TYPE="horizontal")
{
  # global variable DOP.LIST
  LIST <- DOP.LIST[[TYPE]]
  axes <- LIST$axes
  DOP <- LIST$DOP
  VAR <- LIST$VAR
  COV <- LIST$COV
  
  # class indices
  Ci <- get.class(data,names(UERE))
  # partition in to NA-UERE and UERE>0
  NAS <- is.na(UERE)
  POS <- !NAS
  
  # don't overwrite ARGOS when calibrating GPS
  ARGOS <- grepl('argos',tolower(levels(data$class)))
  if(TYPE=='horizontal' && any(ARGOS))
  {
    ARGOS <- levels(data$class)[ARGOS]
    SAFE <- (data$class != ARGOS)
  }
  else
  { SAFE <- rep(TRUE,length(data$t)) }
  
  if(all(NAS))
  {
    # delete the VAR/COV columns
    data[[VAR]] <- NULL
    if(any(COV %in% names(data))) { data[COV] <- NULL }
  }
  else if(any(SAFE))
  {
    # UERE[time] > 0
    if(any(NAS)) { UERE[NAS] <- 0 }
    U.POS <- UERE[Ci]
    
    if(DOP %in% names(data))
    { data[[VAR]][SAFE] <- (U.POS*data[[DOP]])[SAFE]^2/length(axes) }
    else if(!(VAR %in% names(data)))
    {
      data[[VAR]][SAFE] <- (U.POS)[SAFE]^2/length(axes)
      message("DOP values missing. Assuming DOP=1.")
    }
    
    # apply NA-UEREs
    if(any(NAS))
    {
      # indication of NA time
      U.NAS <- NAS[Ci]
      data[[VAR]][U.NAS & SAFE] <- Inf
    }
    
    # covariance matrix (dual tagged)
    if(all(COV %in% names(data)))
    {
      data[[COV[1]]][SAFE] <- data[[VAR]][SAFE]
      data[[COV[2]]][SAFE] <- 0
      data[[COV[3]]][SAFE] <- data[[VAR]][SAFE]
    }
    
    if(TYPE=="horizontal" && all(LIST$COV.geo %in% names(data)))
    {
      data$COV.major[SAFE] <- data[[VAR]][SAFE]
      data$COV.minor[SAFE] <- data[[VAR]][SAFE]
      data$COV.angle[SAFE] <- 0
    }
  }
  
  return(data)
}


# create null UERE object for @UERE slot of telemetry data
uere.null <- function(data)
{
  if("class" %in% names(data))
  { CLASS <- levels(data$class) }
  else
  { CLASS <- "all" }
  
  TYPES <- get.dop.types(data)
  
  M <- matrix(1,length(CLASS),length(TYPES))
  rownames(M) <- CLASS
  colnames(M) <- TYPES
  
  V <- M[1,]
  names(V) <- colnames(V)
  
  UERE <- list(UERE=1*M,DOF=0*M,AICc=Inf*V,Zsq=Inf*V,VAR.Zsq=Inf*V,N=0*V)
  UERE <- new.UERE(UERE)
  
  return(UERE)
}


# store the UERE
# axes determined by name(value) consistent with DOP.LIST global variable above
# NULL is for universal UERE, else "horizontal", "vertical", "speed"
"uere<-" <- function(data,value)
{
  UERE <- value
  
  # default ambiguous assignment - overrides everything
  if(class(UERE)[1]=="numeric" || class(UERE)[1]=="integer")
  {
    # in case of different location classes
    if(class(data)[1]=="list")
    {
      data <- lapply(data,function(d){"uere<-"(d,value)})
      return(data)
    }
    
    # at some point we might want to check the dimensions of value for different types
    
    uere(data) <- NULL
    UERE <- uere(data)
    UERE$UERE[,] <- value
    UERE$DOF[] <- Inf
    UERE$AICc[] <- Inf
    UERE$Zsq[] <- Inf
    UERE$VAR.Zsq[] <- Inf
    UERE$N[] <- Inf
    uere(data) <- UERE
    return(data)
  }
  else if(class(UERE)[1]=="ctmm")
  {
    # in case of different location classes
    if(class(data)[1]=="list")
    {
      data <- lapply(data,function(d){"uere<-"(d,value)})
      return(data)
    }
    
    # only calibrate axis of value
    axes <- value$axes
    AXES <- length(axes)
    TYPE <- DOP.match(axes)
    ERROR <- value$error
    CLASS <- names(ERROR)
    ECLASS <- paste("error",CLASS)
    
    if(all(ECLASS %in% rownames(value$COV)))
    {
      # sqrt(chi^2) relation
      N <- 2/AXES * ERROR^2 / (2^2*diag(value$COV[ECLASS,ECLASS,drop=FALSE]))
      names(N) <- CLASS
    }
    else
    {
      N <- UERE
      N[] <- Inf
    }
    
    UERE <- uere(data)
    CLASS <- CLASS[CLASS %in% rownames(UERE$UERE)]
    ERROR <- ERROR[CLASS]
    N <- N[CLASS]
    
    UERE$UERE[CLASS,TYPE] <- ERROR
    UERE$DOF[CLASS,TYPE] <- N
    UERE$AICc[TYPE] <- Inf
    UERE$Zsq[TYPE] <- Inf
    UERE$VAR.Zsq[TYPE] <- Inf
    UERE$N[TYPE] <- sum(N)
    uere(data) <- UERE
    return(data)
  }
  
  DOF <- UERE$DOF
  
  # promote to list and revert back if DROP
  DROP <- FALSE
  if(class(data)[1]=="telemetry" || class(data)[1]=="data.frame")
  {
    data <- list(data)
    DROP <- TRUE
  }
  
  for(i in 1:length(data))
  {
    # make sure of UERE slot compatibility
    # data[[i]] <- droplevels(data[[i]])
    # why was I doing this?
    
    # all class/type from data columns
    UERE.NULL <- uere.null(data[[i]])
    CLASS <- rownames(UERE.NULL$UERE)
    TYPES <- colnames(UERE.NULL$UERE)
    
    if(length(UERE$UERE))
    {
      # all class/type in @UERE slot
      CLASS <- c(CLASS,rownames(attr(data[[i]],"UERE")$UERE))
      TYPES <- c(TYPES,colnames(attr(data[[i]],"UERE")$UERE))
      
      # all class/type in UERE-value
      CLASS <- c(CLASS,rownames(UERE$UERE))
      TYPES <- c(TYPES,colnames(UERE$UERE))
      
      CLASS <- unique(CLASS)
      TYPES <- unique(TYPES)
    }
    
    # setup null UERE
    UERE.NULL <- matrix(1,length(CLASS),length(TYPES))
    rownames(UERE.NULL) <- CLASS
    colnames(UERE.NULL) <- TYPES
    DOF.NULL <- 0*UERE.NULL
    V <- UERE.NULL[1,]
    names(V) <- TYPES # R drops dimnames...
    AICc.NULL <- Inf*V
    Zsq.NULL <- AICc.NULL
    VAR.Zsq.NULL <- AICc.NULL
    N.NULL <- 0*V
    
    # copy over @UERE slot
    if(length(attr(data[[i]],"UERE")) && length(UERE$UERE))
    {
      CLASS <- rownames(attr(data[[i]],"UERE")$UERE)
      TYPES <- colnames(attr(data[[i]],"UERE")$UERE)
      UERE.NULL[CLASS,TYPES] <- attr(data[[i]],"UERE")$UERE
      DOF.NULL[CLASS,TYPES] <- attr(data[[i]],"UERE")$DOF
      AICc.NULL[TYPES] <- attr(data[[i]],"UERE")$AICc
      Zsq.NULL[TYPES] <- attr(data[[i]],"UERE")$Zsq
      VAR.Zsq.NULL[TYPES] <- attr(data[[i]],"UERE")$VAR.Zsq
      N.NULL[TYPES] <- attr(data[[i]],"UERE")$N
    }
    
    # copy over UERE-value (overwriting conflicts)
    if(length(UERE))
    {
      CLASS <- rownames(UERE$UERE)
      TYPES <- colnames(UERE$UERE)
      UERE.NULL[CLASS,TYPES] <- UERE$UERE
      DOF.NULL[CLASS,TYPES] <- UERE$DOF
      AICc.NULL[TYPES] <- UERE$AICc
      Zsq.NULL[TYPES] <- UERE$Zsq
      VAR.Zsq.NULL[TYPES] <- UERE$VAR.Zsq
      N.NULL[TYPES] <- UERE$N
    }
    
    # calibrate data
    TYPES <- get.dop.types(data[[i]])
    for(TYPE in TYPES)
    {
      # R does not preserve dimension names ???
      UERE.ROW <- UERE.NULL[,TYPE]
      names(UERE.ROW) <- rownames(UERE.NULL)
      # don't calibrate data if NULL assignment
      if(!is.null(value)) { data[[i]] <- try.assign.uere(data[[i]],UERE.ROW,TYPE=TYPE) }
      else { data[[i]] <- try.assign.uere(data[[i]],NA*UERE.ROW,TYPE=TYPE) }
    }
    
    # store UERE that was applied
    UERE.NULL <- list(UERE=UERE.NULL,DOF=DOF.NULL,AICc=AICc.NULL,Zsq=Zsq.NULL,VAR.Zsq=VAR.Zsq.NULL,N=N.NULL)
    UERE.NULL <- new.UERE(UERE.NULL)
    attr(data[[i]],"UERE") <- UERE.NULL
  }
  
  if(DROP) { data <- data[[1]] }
  
  # demote to data.frame
  return(data)
}


#####
# class index at times
get.class <- function(data,LEVELS=levels(data$class))
{
  if("class" %in% names(data))
  {
    Ci <- rep(NA_integer_,nrow(data))
    
    for(i in 1:length(LEVELS))
    {
      SUB <- which(data$class==LEVELS[i])
      if(length(SUB)) { Ci[SUB] <- i }
    }
  }
  else
  { Ci <- rep(1,length(data$t)) }
  
  return(Ci)
}


#####
# class indicator matrix - [time,class]
get.class.mat <- function(data,LEVELS=levels(data$class))
{
  if("class" %in% names(data))
  {
    C <- sapply(LEVELS,function(lc){data$class==lc}) # (time,class)
    dim(C) <- c(nrow(data),length(LEVELS))
    colnames(C) <- LEVELS
  }
  else
  { C <- cbind(rep(1,length(data$t))) }
  
  return(C)
}


########
get.UERE.DOF <- function(x)
{
  N <- x$DOF
  if(is.null(N) || is.na(N)) { N <- 0 }
  return(N)
}


############
# get/set dimnames of UERE objects
classnames <- function(object)
{ rownames(object$UERE) }

"classnames<-" <- function(object,value)
{
  rownames(object$UERE) <- value
  rownames(object$DOF) <- value
  return(object)
}

"typenames<-" <- function(object,value)
{
  colnames(object$UERE) <- value
  colnames(object$DOF) <- value
  names(object$AICc) <- value
  names(object$Zsq) <- value
  names(object$VAR.Zsq) <- value
  names(object$N) <- value
  return(object)
}

getMethod <- function(fn,signature,...)
{
  meth <- NULL
  # S3 and S4 method dispatching is incompatible for no reason
  meth <- methods::getMethod(fn,signature=signature,optional=TRUE,...) # try S4 first
  if(is.null(meth)) { meth <- utils::getS3method(fn,signature[1],optional=TRUE,...) } # then try S3
  # due to new CRAN policy, we can no longer have internal S3 methods
  # work around here with '_' dispatch method names
  if(is.null(meth)) { try( meth <- get(paste0(fn,"_",signature)) , silent=TRUE) }
  if(is.null(meth)) { stop('Cannot find method ',fn,' for class ',signature) }
  return(meth)
}


# base match.arg cannot match NA ???
match.arg <- function(arg,choices,...)
{
  if(is.na(arg) && any(is.na(choices))) { return(NA) }
  else { return(base::match.arg(arg,choices,...)) }
}

# sort arguments by class
# sort.arg <- function(arg,sig)
# {
#
# }


# does this thing exist and, if so, is it true
is.good <- function(x) { !is.null(x) & !is.na(x) & x }
is.bad <- function(x) { is.null(x) | is.na(x) | !x }

# not in #
"%nin%" <- function(x, table) { match(x, table, nomatch = 0) == 0 }


# positive only sequence for for() loops so that for(i in 1:0) does nothing
"%:%" <- function(x,y)
{
  if(x<=y) { x <- x:y } else { x <- NULL }
  return(x)
}


composite <- function(n) { 2^ceiling(log(n,2)) }


# sinc functions
sinc <- Vectorize( function(x,SIN=sin(x))
{
  if(x==0)
  { return(1) }
  else
  { return(SIN/x) }
} )

sinch <- Vectorize( function(x,SINH=sinh(x))
{
  if(x==0)
  { return(1) }
  else
  { return(SINH/x) }
} )


# multivariate polygamma function
mpsigamma <- function(x,deriv=0,dim=1)
{
  PSI <- 1 - 1:dim
  PSI <- x + PSI/2
  if(deriv>=0) { PSI <- sapply(PSI,function(p) psigamma(p,deriv=deriv)) }
  else if(deriv==-1) { PSI <- sapply(PSI,function(p) lgamma(p)) }
  else { stop("Derivative ",deriv+1," of log(Gamma(x)) not supported.") }
  PSI <- sum(PSI)
  return(PSI)
}

### S4 ###


### S3 ###

# forwarding function for list of a particular datatype
log.list <- function(x,...)
{
  CLASS <- class(x[[1]])[1]
  getMethod("log",CLASS)(x,...)
}
# this doesn't work outside of ctmm


# forwarding function for list of a particular datatype
mean.list <- function(x,...)
{
  CLASS <- class(x[[1]])[1]
  getMethod("mean",CLASS)(x,...)
}
#methods::setMethod("mean",signature(x="list"), function(x,...) mean.list(x,...))

# forwarding function for list of a particular datatype
median.list <- function(x,na.rm=FALSE,...)
{
  CLASS <- class(x[[1]])[1]
  getMethod("median",CLASS)(x,na.rm=na.rm,...)
}

# forwarding function for list of a particular datatype
plot.list <- function(x,...)
{
  CLASS <- class(x[[1]])[1]
  getMethod("plot",CLASS)(x,...)
}
#methods::setMethod("plot",signature(x="list"), function(x,...) plot.list(x,...))

# forwarding function for list of a particular datatype
summary.list <- function(object,...)
{
  # recurse if necessary
  CLASS <- "list"
  DATA <- object
  while(CLASS=="list")
  {
    DATA <- DATA[[1]]
    CLASS <- class(DATA)
  }
  
  getMethod("summary",CLASS)(object,...)
}

# forwarding function for list of a particular datatype
writeVector.list <- function(x,filename,...)
{
  CLASS <- class(x[[1]])[1]
  if(missing(filename))
  { getMethod("writeVector",methods::signature(x=CLASS,filename="missing"))(x,filename=filename,...) }
  else
  { getMethod("writeVector",methods::signature(x=CLASS,filename="character"))(x,filename=filename,...) }
}
methods::setMethod("writeVector",methods::signature(x="list",filename="character"), function(x,filename,...) writeVector.list(x,filename,...) )
methods::setMethod("writeVector",methods::signature(x="list",filename="missing"), function(x,filename,...) writeVector.list(x,filename,...) )


# replace NA elements
na.replace <- function(x,rep)
{
  REP <- is.na(x)
  x[REP] <- rep[REP]
  return(x)
}

########################
# 0/0 -> NaN -> to
# fixes a priori known limits
nant <- function(x,to)
{
  NAN <- is.na(x) # TRUE for NaN and NA
  if(any(NAN))
  {
    to <- array(to,length(x))
    x[NAN] <- to[NAN]
  }
  return(x)
}

# fix for infite PD matrix
# useful after nant(x,Inf)
inft <- function(x,to=0)
{
  INF <- diag(x)==Inf
  if(any(INF))
  {
    # force positive definite
    x[INF,] <- x[,INF] <- 0
    # restore infinite variance
    diag(x)[INF] <- Inf
  }
  return(x)
}


# parity tests
is.even <- Vectorize(function(x) {x %% 2 == 0})

is.odd <- Vectorize(function(x) {x %% 2 != 0})


# last element of array
last <- function(vec) { vec[length(vec)] }
first <- function(vec) { vec[1] }
# assign to last element... doesn't work
# "last<-" <- function(vec,ass)
# {
#   vec[length(vec)] <- ass
#   return(vec)
# }


# prepend to a vector
prepend <- function(x,values,before=1)
{ append(x,values,after=before-1) }


# CLAMP A NUMBER
clamp <- function(num,min=0,max=1)
{ ifelse(num<min,min,ifelse(num>max,max,num)) }


# PAD VECTOR
pad <- function(vec,size,padding=0,side="right")
{
  # this is now the pad length instead of total length
  size <- size - length(vec)
  padding <- array(padding,size)
  
  if(side=="right"||side=="r")
  { vec <- c(vec,padding) }
  else if(side=="left"||side=="l")
  { vec <- c(padding,vec) }
  
  return(vec)
}

# row pad for data frames / matrices
rpad <- function(mat,size,padding=0,side="right")
{
  mat <- cbind(mat)
  size <- size - nrow(mat)
  COL <- ncol(mat)
  padding <- array(padding,c(size,COL))
  colnames(padding) <- colnames(mat)
  
  if(side=="right"||side=="r")
  { mat <- rbind(mat,padding) }
  else if(side=="left" || side=="l")
  { mat <- rbind(padding,mat) }
  
  return(mat)
}

#remove rows and columns by name
rm.name <- function(object,name)
{
  if(length(dim(object))==2)
  { object <- object[! rownames(object) %in% name,! colnames(object) %in% name,drop=FALSE] }
  else
  { object <- object[! names(object) %in% name] }
  return(object)
}


# put in a list if not in a list already
listify <- function(x)
{
  if(is.null(x)) { return(x) }
  
  if(class(x)[1] != "list")
  {
    x <- list(x)
    names(x) <- attr(x[[1]],'info')$identity
  }
  return(x)
}


# rename elements of an object
rename <- function(object,name1,name2)
{
  IND <- which(names(object)==name1)
  names(object)[IND] <- name2
  return(object)
}

rename.matrix <- function(object,name1,name2)
{
  NAMES <- dimnames(object)[[1]]
  IND <- which(NAMES==name1)
  NAMES[IND] <- name2
  dimnames(object) <- list(NAMES,NAMES)
  return(object)
}


# glue strings together if they different
glue <- function(...)
{
  x <- c(...) # removes NULLs
  x <- unique(x) # removes matchs
  x <- paste(x) # paste different strings together
  return(x)
}


# from
capitalize <- function(s)
{
  substr(s,1,1) <- toupper(substr(s,1,1))
  s
}


simplify.formula <- function(x)
{
  #x <- as.character(x)[2]
  x <- stats::terms(x,simplify=TRUE)
  x <- as.character(x)[2]
  x <- paste("~",x)
  x <- eval(parse(text=x))
  return(x)
}


copy <- function(from,to)
{
  NAMES <- names(from)
  for(n in NAMES){ to[[n]] <- from[[n]] }
  return(to)
}

mid <- function(x)
{
  n <- length(x)
  (x[-1]+x[-n])/2
}


name.list <- function(x)
{
  if(class(x)[1]=="list" && !length(names(x)))
  {
    NAMES <- sapply(x,function(y){attr(y,"info")$identity})
    if(class(NAMES)[1]=="character") { names(x) <- NAMES }
  }
  return(x)
}

# global variable
DATA.EARTH <- list(R.EQ=6378137,R.PL=6356752.3142) # equatorial & polar radii

#setGeneric("projection", function(x,...) { standardGeneric("projection") }, signature='x')
#setGeneric("projection<-", function(x,value,...) { standardGeneric("projection<-") }, signature='x')

# setMethod('projection',signature(x='Raster'),raster::projection)
# setMethod('projection',signature(x='RasterLayer'),raster::projection)
# setMethod('projection',signature(x='RasterStack'),raster::projection)
# setMethod('projection',signature(x='RasterBrick'),raster::projection)
# setMethod('projection<-',signature(x='Raster'),raster::projection)
# setMethod('projection<-',signature(x='RasterLayer'),raster::projection)
# setMethod('projection<-',signature(x='RasterStack'),raster::projection)
# setMethod('projection<-',signature(x='RasterBrick'),raster::projection)

# range of telemetry data
projection.telemetry <- function(x,asText=TRUE)
{
  proj <- attr(x,"info")$projection
  if(!asText) { proj <- sp::CRS(proj,doCheckCRSArgs=FALSE) }
  return(proj)
}

project <- function(x,from=DATUM,to=DATUM)
{
  if(to==from) { return(x) }
  
  # x <- sp::SpatialPoints(x,proj4string=sp::CRS(from))
  # x <- sp::spTransform(x,sp::CRS(to))
  # x <- sp::coordinates(x)
  
  x <- data.frame(x) # super annoying
  from <- sf::st_crs(from)
  to <- sf::st_crs(to)
  
  x <- sf::st_as_sf(x,coords=1:2,crs=from)
  x <- sf::st_transform(x,crs=to)
  x <- sf::st_coordinates(x)
  
  colnames(x) <- c('x','y')
  
  return(x)
}


# change the projection of one telemetry object

# unit vector pointing north at projected locations R
# x = partially projected telemetry data dim(n,2)
# return dim(n,2)
northing <- function(x,proj,angle=FALSE)
{
  # WGS-84 ellipsoid
  R.EQ <- DATA.EARTH$R.EQ
  R.PL <- DATA.EARTH$R.PL
  # approximate 1-meter-North latitude displacements
  d.lambda <- 1/sqrt((R.EQ*sin(x$latitude))^2+(R.PL*cos(x$latitude))^2)
  # could use grad() but would be slowwwww....
  d.lambda <- d.lambda*(360/2/pi) # arg, degrees!
  u <- cbind(x$longitude,x$latitude + d.lambda)
  u <- project(u,to=proj)
  colnames(u) <- c("x","y")
  # difference vectors pointing North ~1 meters
  u <- u - get.telemetry(x) # [n,2]
  # difference vectors pointing North 1 meters exact
  u <- u / sqrt(rowSums(u^2)) # [n,2]
  
  if(angle) { u <- atan2(u[,'y'],u[,'x']) * (360/(2*pi)) } # R plotting functions require degrees
  
  return(u)
}


# rotate northing to heading
# u = dim(2,n)
# return = dim(n,2)
rotate.north <- function(u,heading)
{
  heading <- heading * (2*pi/360) # ack, degrees!
  
  u <- u[,1] + 1i*u[,2] # velocity vector
  R <- exp(-1i*heading) # rotation matrix
  u <- R*u
  u <- cbind(Re(u),Im(u))
  
  return(u)
}


# put projection into character format
format_projection <- function(proj,datum="WGS84")
{
  if(class(proj)[1]=="CRS")
  { proj <- as.character(proj) }
  else if(class(proj)[1] != "character")
  {
    # pull out geodesic coordinates and format into matrix
    proj <- as.matrix(rbind(proj)[,c("longitude","latitude")])
    
    if(nrow(proj)==1)
    { proj <- paste0("+proj=aeqd  +lon_0=",proj[1,1]," +lat_0=",proj[1,2]," +datum=",datum) }
    else if(nrow(proj)==2)
    { proj <- paste0("+proj=tpeqd +lon_1=",proj[1,1]," +lat_1=",proj[1,2]," +lon_2=",proj[2,1]," +lat_2=",proj[2,2]," +datum=",datum) }
    else if(nrow(proj)==3)
    { proj <- paste0("+proj=chamb +lon_1=",proj[1,1]," +lat_1=",proj[1,2]," +lon_2=",proj[2,1]," +lat_2=",proj[2,2]," +lon_3=",proj[3,1]," +lat_3=",proj[3,2]," +datum=",datum) }
    else
    { stop("PROJ4 does not support ",nrow(proj)," foci projections.") }
  }
  
  # put in canonical format
  proj <- sp::CRS(proj)
  proj <- as.character(proj)
  
  validate.projection(proj)
  return(proj)
}


# only allow compatible projections
validate.projection <- function(projection)
{
  if(class(projection)[1]=="character") { projection <- sp::CRS(projection) } # this adds missing longlat specification
  if(class(projection)[1]=="CRS") { projection <- as.character(projection) }
  
  if(grepl("longlat",projection,fixed=TRUE) || grepl("latlong",projection,fixed=TRUE))
  { stop("A projected coordinate system must be specified.") }
  
  if(grepl("units=",projection,fixed=TRUE) && !grepl("units=m",projection,fixed=TRUE))
  { stop("Units of distance other than meters not supported.") }
}


# projection check data against grid
validate.grid <- function(data,grid)
{
  if(class(grid)[1] %in% c("UD","RasterLayer") && !is.null(projection(data)) && !is.null(projection(grid)))
  {
    proj1 <- projection(data)
    proj2 <- projection(grid)
    
    # put into canonical format
    proj1 <- sp::CRS(proj1)
    proj2 <- sp::CRS(proj2)
    
    proj1 <- as.character(proj1)
    proj2 <- as.character(proj2)
    
    if(proj1 != proj2) { stop("Grid projection does not match data projection.") }
  }
}


# check for consistent projections
check.projections <- function(object)
{
  PROJ <- projection(object)
  if(length(PROJ)==0) { stop("Missing projection.") }
  if(length(PROJ)>1) { stop("Inconsistent projections.") }
  return(PROJ)
}

# median of dataset
median.telemetry <- function(x,na.rm=FALSE,...)
{
  if(class(x)[1]=="list")
  {
    if(length(x)>1)
    {
      id <- mean_info(x)$identity
      # is there a common projection to preserve
      proj <- sapply(x,projection.telemetry)
      proj <- unlist(proj) # need for multiple NULLs
      
      if(length(proj)<length(x))
      { proj <- NULL }
      else
      {
        proj <- unique(proj)
        if(length(proj)>1) { proj <- NULL }
      }
      
      x <- lapply(x,function(d){ median.telemetry(d,...) }) # this will only ever return long-lat & x-y columns
      x <- do.call(rbind,x)
      x <- as.data.frame(x)
      x <- ctmm::new.telemetry(x,info=list(identity=id,projection=proj),UERE=new.UERE())
    }
    else
    { x <- x[[1]] }
    # now flow through to per-individual analysis below
  }
  
  id <- paste0("median of ",attr(x,'info')$identity)
  # tz <- attr(x,'info')$timezone
  
  # t <- stats::median(x$t)
  # timestamp <- as.character(as.POSIXct(t,tz=tz,origin="1970/01/01"))
  
  if(all(c("longitude","latitude") %in% names(x)))
  {
    proj <- projection(x)
    
    mu <- median_longlat(x)
    x <- data.frame(longitude=mu[,"longitude"],latitude=mu[,"latitude"])
    x <- new.telemetry(x,info=list(identity=id,projection=proj),UERE=new.UERE())
    
    ctmm::projection(x) <- proj # NULL is handled
  }
  else # fallback for turtle data (k>1 not supported)
  {
    x <- cbind(x$x,x$y)
    mu <- apply(x,2,stats::median)
    mu <- c(Gmedian::Gmedian(x,init=mu,...))
    
    x <- data.frame(x=mu[1],y=mu[2])
    x <- new.telemetry(x,info=list(identity=id),UERE=new.UERE())
  }
  
  return(x)
}

median_longlat <- function(data,k=1,...)
{
  data <- ellipsoid2cartesian(data)
  
  # what is the centroid of the data in 3D
  mu <- apply(data,2,stats::median)
  STAT <- Gmedian::GmedianCov(data,init=mu,scores=0,nstart=10,...)
  mu <- c(STAT$median)
  COV <- STAT$covmedian
  
  if(k==1 || nrow(data)==1)
  {
    mu <- cartesian2ellipsoid(mu)
    return(mu)
  }
  
  # k==2 below
  if(nrow(data)==1)
  { mu <- rbind(data,data) }
  else if(nrow(data)==2)
  { mu <- data }
  else
  {
    # find the longest axis of variation
    COV <- eigen(COV)$vectors[,1]
    # distances along long axis
    COV <- t(t(data)-mu) %*% COV
    # centroid along long axis
    mu <- stats::median(COV)
    # detrended distances along long axis
    COV <- c(COV) - mu
    
    # positive half
    SUB <- data[COV>=0,,drop=FALSE]
    if(nrow(SUB)==1) { mu1 <- SUB }
    else
    {
      mu1 <- apply(SUB,2,stats::median)
      mu1 <- c(Gmedian::Gmedian(SUB,init=mu1))
    }
    
    # negative half
    SUB <- data[COV<=0,,drop=FALSE]
    if(nrow(SUB)==1) { mu2 <- SUB }
    else
    {
      mu2 <- apply(SUB,2,stats::median)
      mu2 <- c(Gmedian::Gmedian(SUB,init=mu2))
    }
    
    # 2-mode cluster
    mu <- rbind(mu1,mu2)
    if(nrow(mu)>3) { mu <- Gmedian::kGmedian(data,ncenters=mu,nstart=10,...)$centers }
    
    # make sure separation is not reduced by clustering
    if(sum((mu1-mu2)^2)>sum((mu[1,]-mu[2,])^2))
    {
      mu1 -> mu[1,]
      mu2 -> mu[2,]
    }
  }
  mu <- cartesian2ellipsoid(mu)
  
  # order from west to east
  if(mu[1,1] > mu[2,1])
  {
    mu[1,] -> mu2
    mu[2,] -> mu1
    
    mu <- rbind(mu1,mu2)
  }
  
  colnames(mu) <- c("longitude","latitude")
  return(mu)
}

ellipsoid2cartesian <- function(data)
{
  # assume Movebank data.frame
  lon <- data$longitude * (2*pi/360)
  lat <- data$latitude * (2*pi/360)
  s <- data$z # distance from surface
  if(is.null(s)) { s <- 0 }
  
  # convert to x,y,z coordinates (3D)
  z <- (DATA.EARTH$R.PL + s) * sin(lat)
  r <- (DATA.EARTH$R.EQ + s) * cos(lat)
  
  x <- r * cos(lon)
  y <- r * sin(lon)
  
  data <- cbind(x,y,z)
  return(data)
}
cartesian2ellipsoid <- function(mu)
{
  mu <- rbind(mu)
  
  x <- mu[,1]
  y <- mu[,2]
  z <- mu[,3]
  
  lon <- atan2(y,x)
  
  r <- sqrt(x^2 + y^2)
  
  r <- r/DATA.EARTH$R.EQ
  z <- z/DATA.EARTH$R.PL
  
  lat <- atan2(z,r)
  
  mu <- cbind(longitude=lon,latitude=lat)  * (360/(2*pi))
  
  return(mu)
}

get.telemetry <- function(data,axes=c("x","y"))
{
  # z <- "[.data.frame"(data,axes)
  # z <- as.matrix(z)
  # colnames(z) <- axes
  # return(z)
  if(all(axes %in% names(data)))
  { data <- as.matrix(data.frame(data)[, axes], dimnames = axes) }
  else
  { data <- numeric(0) }
  
  return(data)
}

unit.telemetry <- function(data,length=1,time=1,axes=c('x','y'))
{
  if(class(data)[1]=="phylometry")
  {
    lag <- attr(data,"lag")/time
    data[,axes] <- data[,axes]/length
    attr(data,"lag") <- lag
    # class(data) <- "phylometry"
    return(data)
  }
  # TELEMETRY CLASS BELOW
  
  convert <- function(NAMES,scale) { for(NAME in NAMES) { if(NAME %in% names(data)) { data[[NAME]] <<- data[[NAME]]/scale } } }
  
  convert('t',time)
  convert('light.time',time)
  convert('dark.time',time)
  convert('sundial.rate',1/time)
  
  if(any(axes %in% c('x','y','z')))
  {
    convert(DOP.LIST$horizontal$axes,length)
    convert(DOP.LIST$vertical$axes,length)
    convert(DOP.LIST$speed$axes,length/time)
    
    convert(DOP.LIST$horizontal$VAR,length^2)
    convert(DOP.LIST$horizontal$COV,length^2)
    
    convert(DOP.LIST$vertical$VAR,length^2)
    
    convert(DOP.LIST$speed$VAR,(length/time)^2)
    convert(DOP.LIST$speed$COV,(length/time)^2)
  }
  else
  { for(axis in axes) { convert(axis,length) } }
  
  # calibration constants
  attr(data,"UERE")$UERE <- attr(data,"UERE")$UERE/length # don't logicals get divided this way?
  
  return(data)
}

plot.env <- new.env()
new.plot <- function(data=NULL,CTMM=NULL,UD=NULL,R=NULL,col.bg="white",col.R="green",legend=FALSE,level.UD=0.95,level=0.95,units=TRUE,fraction=1,add=FALSE,xlim=NULL,ylim=NULL,ext=NULL,cex=NULL,...)
{
  RESIDUALS <- !is.null(data) && !is.null(attr(data[[1]],"info")$residual)
  
  if(is.null(cex)) { cex <- graphics::par('cex') }
  if(class(cex)[1]=='list') { cex <- sapply(cex,stats::median) }
  cex <- stats::median(cex)
  
  dist <- list()
  dist$name <- "meters"
  dist$scale <- 1
  axes <- c("x","y")
  
  if(!is.null(ext))
  {
    xlim <- ext$x
    ylim <- ext$y
  }
  
  if(!add)
  {
    ext <- NULL
    
    # bounding locations from data
    if(!is.null(data))
    { ext <- rbind(ext,extent(data)[,axes]) }
    
    # bounding locations from UDs
    if(!is.null(UD))
    {
      if(class(UD[[1]])[1]=='RS') # backup extent for RS objects
      {
        ext <- data.frame(x=1:2,y=1:2)
        rownames(ext) <- c('min','max')
        ext['min','x'] <- min(sapply(UD,function(U){U$r$x[1]}))
        ext['max','x'] <- min(sapply(UD,function(U){last(U$r$x)}))
        ext['min','y'] <- min(sapply(UD,function(U){U$r$y[1]}))
        ext['max','y'] <- min(sapply(UD,function(U){last(U$r$y)}))
      }
      else
      { ext <- rbind(ext,extent(UD,level=level,level.UD=level.UD)[,axes]) }
    }
    
    # bounding locations from Gaussian CTMM
    if(!is.null(CTMM) & !any(is.na(level.UD)))
    { ext <- rbind(ext,extent(CTMM,level=level,level.UD=level.UD)[,axes]) }
    
    # # bounding locations from standard normal quantiles
    # if(RESIDUALS && !any(is.na(level.UD)))
    # {
    #   Z <- c(1,1) * sqrt(-2*log(1-max(level.UD)))
    #   ext <- rbind(ext,-Z)
    #   ext <- rbind(ext,+Z)
    # }
    
    # combine ranges
    ext <- extent.telemetry(ext)
    
    # bounding box
    mu <- c(mean(ext$x),mean(ext$y))
    buff <- c(diff(ext$x),diff(ext$y))/2
    
    # now zoom in/out to some fraction of the grid
    buff <- fraction*buff
    
    ext$x <- mu[1] + buff[1]*c(-1,1)
    ext$y <- mu[2] + buff[2]*c(-1,1)
    
    # try to obey xlim/ylim if provided
    if(!is.null(xlim) || !is.null(ylim))
    {
      max.diff <- max(diff(xlim),diff(ylim))*c(-1,1)/2
      
      if(is.null(ylim))
      { ylim <- mu[2] + max.diff }
      else if(is.null(xlim))
      { xlim <- mu[1] + max.diff }
      
      ext$x <- xlim
      ext$y <- ylim
    }
    
    # Get best unit scale
    dist <- unit(unlist(ext),"length",SI=!units)
    
    xlab <- paste("x ", "(", dist$name, ")", sep="")
    ylab <- paste("y ", "(", dist$name, ")", sep="")
    
    # residuals have no units
    if(RESIDUALS)
    {
      xlab <- "x"
      ylab <- "y"
    }
    
    ext <- ext/dist$scale
    
    # empty base layer plot
    plot(ext, xlab=xlab, ylab=ylab, col=grDevices::rgb(1,1,1,0), asp=1, cex=cex)
    
    # plot background color
    lim <- graphics::par('usr')
    xlim <- lim[1:2]
    dx <- diff(xlim)
    ylim <- lim[3:4]
    dy <- diff(ylim)
    # dx <- dy <- 0
    graphics::rect(xlim[1]-dx/2,ylim[1]-dy/2,xlim[2]+dx/2,ylim[2]+dy/2,border=col.bg,col=col.bg)
    # this can cover the plot box
    graphics::box()
    
    # plot information for further layering
    projection <- unique(c(projection(data),projection(CTMM),projection(UD))) # some objects could be NULL
    if(length(projection)>1 && !RESIDUALS) { stop("Multiple projections not yet supported.") }
    assign("projection",projection,envir=plot.env)
    # dimensional type
    assign("x.dim","length",pos=plot.env)
    assign("y.dim","length",pos=plot.env)
    # unit name
    assign("x.units",dist$name,pos=plot.env)
    assign("y.units",dist$name,pos=plot.env)
    # unit conversion
    assign("x.scale",dist$scale,pos=plot.env)
    assign("y.scale",dist$scale,pos=plot.env)
    
    scale <- dist$scale
  } # end !add
  else # get distance information from environment
  {
    name <- unique( c( get0('x.units',plot.env), get0('y.units',plot.env) ) )
    scale <- unique( c( get0('x.scale',plot.env), get0('y.scale',plot.env) ) )
    if(length(name)==1 && length(scale)==1)
    {
      dist$name <- name
      dist$scale <- scale
    }
    projection <- get('projection',plot.env)
  }
  
  # PLOT RASTER / SUITABILITY
  if(!is.null(R))
  { plot_R(R,col=col.R,legend=legend,projection=projection,scale=scale) }
  
  return(dist)
}

extent.telemetry <- function(x,level=1,...)
{
  level <- max(level)
  
  alpha <- (1-level)/2
  probs <- c(alpha,1-alpha)
  RANGE <- data.frame(row.names=c('min','max'))
  
  for(COL in colnames(x))
  {
    if(class(x[[COL]])[1]=="numeric")
    {
      if(COL=="longitude") # use circular statistics
      { RANGE[[COL]] <- quantile.longitude(x[[COL]],probs=probs,na.rm=TRUE) }
      else
      { RANGE[[COL]] <- stats::quantile(x[[COL]],probs=probs,na.rm=TRUE) }
    }
  }
  
  return(RANGE)
}
setMethod('extent', signature(x='telemetry'), extent.telemetry)
setMethod('extent', signature(x='data.frame'), extent.telemetry)


# extent of matrix/extent
extent.matrix <- function(x,level=1,...)
{
  level <- max(level)
  
  NAMES <- colnames(x)
  x <- data.frame(x)
  colnames(x) <- NAMES # preserve NULL
  extent(x,level=level,...)
}

# some default and non-default units
UNITS <- list()
UNITS$mean <- list()
UNITS$mean$DAY <- 86400.002 # mean solar day
UNITS$mean$YEAR <- 365.24217 * UNITS$mean$DAY # mean tropical year
UNITS$mean$MONTH <- 2.8 + 60*( 44 + 60*( 12 + 24*29 ) ) # mean synodic month
UNITS$calendar <- list()
UNITS$calendar$DAY <- 86400
UNITS$calendar$YEAR <- 365 * UNITS$calendar$DAY
UNITS$calendar$MONTH <- 30 * UNITS$calendar$DAY

# create units dictionary
generate.units <- function()
{
  alias <- list()
  scale <- c()
  
  add <- function(a,s)
  {
    n <- length(alias)
    alias[[n+1]] <<- canonical.name(a)
    scale[n+1] <<- s
  }
  
  OP <- getOption("time.units")
  
  # TIME
  add(c("\u03BCs","microsecond","microseconds"),1E-6)
  add(c("ms","milisecond","miliseconds"),1/1000)
  add(c("s","sec","sec.","second","seconds"),1)
  add(c("min","minute","minutes"),60)
  add(c("h","hr","hour","hours"),60^2)
  add(c("day","days"),UNITS[[OP]]$DAY) #day
  add(c("wk","week","weeks"),7*UNITS[[OP]]$DAY) # week
  add(c("mon","month","months"),UNITS[[OP]]$MONTH) # month
  add(c("yr","year","years"),UNITS[[OP]]$YEAR) # year
  add(c("ka","ky","millennium","millenniums","millennia","kiloannum","kiloannums","kiloyear","kiloyears"),1000*UNITS[[OP]]$YEAR)
  add(c("ma","my","megaannum","megaannums","megaanna","megayear","megayears","millionennium","millionenniums","millionennia"),1000^2*UNITS[[OP]]$YEAR)
  add(c("ae","ga","gy","gyr","aeon","aeons","eon","eons","gigayear","gigayears","gigaannum","gigaannums","giggaanna"),1000^3*UNITS[[OP]]$YEAR)
  
  # Distance conversions
  add(c("\u03BCm","micron","microns","micrometer","micrometers"),1E-6)
  add(c("mm","milimeter","milimeters"),1/1000)
  add(c("cm","centimeter","centimeters"),1/100)
  add(c("m","meter","meters"),1)
  add(c("km","kilometer","kilometers"),1000)
  add(c("in","inch","inches"),0.3048/12)
  add(c("ft","foot","feet"),0.3048)
  add(c("yd","yard","yards"),0.3048*3)
  add(c("mi","mile","miles"),0.3048*5280)
  
  # Area conversions
  add(c("\u03BCm\u00B2","micron\u00B2","microns\u00B2","micrometer\u00B2","micrometers\u00B2","\u03BCm^2","\u03BCm.^2","micron^2","microns^2","micrometer^2","micrometers^2","square micron","square microns","square micrometer","square micrometers","micron squared","microns squared","micrometer squared","micrometers squared"),1E-12)
  add(c("mm\u00B2","milimeter\u00B2","milimeters\u00B2","mm^2","mm.^2","milimeter^2","milimeters^2","square milimeter","square milimeters","milimeter squared","milimeters squared"),1/1000^2)
  add(c("cm\u00B2","centimeter\u00B2","centimeters\u00B2","cm^2","cm.^2","centimeter^2","centimeters^2","square centimeter","square centimeters","centimeter squared","centimeters squared"),1/100^2)
  add(c("m\u00B2","meter\u00B2","meters\u00B2","m^2","m.^2","meter^2","meters^2","square meter","square meters","meter squared","meters squared"),1)
  add(c("ha","hectare","hectares","hm\u00B2","hectometer\u00B2","hectometre\u00B2","hectometers\u00B2","hectometres\u00B2","hm^2","hectometer^2","hectometre^2","hectometers^2","hectometres^2","square hm","square hectometer","square hectometre","square hectometers","square hectometres"),100^2)
  add(c("km\u00B2","kilometer\u00B2","kilometers\u00B2","km^2","km.^2","kilometer^2","kilometers^2","square kilometer","square kilometers","kilometer squared","kilometers squared"),1000^2)
  add(c("in\u00B2","inch\u00B2","inches\u00B2","in^2","in.^2","inch^2","inches^2"),(0.3048/12)^2)
  add(c("ft\u00B2","foot\u00B2","feet\u00B2","ft^2","ft.^2","foot^2","feet^2","square foot","square feet","foot squared","feet squared"),0.3048^2)
  add(c("yd\u00B2","yard\u00B2","yards\u00B2","yd^2","yd.^2","yard^2","yards^2","square yard","square yards","yard squared","yards squared"),(0.3048*3)^2)
  add(c("mi\u00B2","mile\u00B2","miles\u00B2","mi^2","mi.^2","mile^2","miles^2","square mile","square miles","mile squared","miles squared"),(0.3048*5280)^2)
  
  # speed
  add(c("mps","m/s","m/sec","meter/sec","meter/second","meters/second"),1)
  add(c("kmph","kph","km/h","km/h","km/hr","kilometer/hour","kilometers/hour"),0.277777777777777777777)
  add(c("mph","mi/h","mi/hr","mile/h","mile/hr","mile/hour","miles/hour"),0.44704)
  add(c("fps","ft/s","ft/sec","feet/second"),0.3048)
  add(c('kt','kn','knot','knots'),1.852 * 0.277777777777777777777)
  add(c("km/s","kmps","km/sec"),1000)
  add(c("cm/s","cmps","cm/sec"),1/100)
  add(c("mm/s","mmps","mm/sec"),1/1000)
  add(c("\u03BCm/s","\u03BCmps","\u03BCm/sec"),1E-6)
  
  # frequency
  add(c("Hz","hertz"),1)
  add(c("kHz","kilohertz"),1000)
  add(c("MHz","megahertz"),1000^2)
  add(c("GHz","gigahertz"),1000^3)
  add(c("THz","terahertz"),1000^4)
  add(c("per min","1/min","min^-1","min\u207B\u00B9","per minute","1/minute","minute^-1","minute\u207B\u00B9"),1/60)
  add(c("per hr","1/hr","hr^-1","hr\u207B\u00B9","per hour","1/hour","hour^-1","hour\u207B\u00B9"),1/60^2)
  add(c("per day","1/day","day^-1","day\u207B\u00B9"),1/UNITS[[OP]]$DAY)
  add(c("per mon","1/mon","mon^-1","mon\u207B\u00B9","per month","1/month","month^-1","month\u207B\u00B9"),1/UNITS[[OP]]$MONTH)
  add(c("per yr","1/yr","yr^-1","yr\u207B\u00B9","per year","1/year","year^-1","year\u207B\u00B9"),1/UNITS[[OP]]$YEAR)
  add(c("per ky","1/ky","ky^-1","ky\u207B\u00B9","per ka","1/ka","ka^-1","ka\u207B\u00B9","per millennium","1/millennium","millennium^-1","millennium\u207B\u00B9","per kiloannum","1/kiloannum","kiloannum^-1","kiloannum\u207B\u00B9","per kiloyear","1/kiloyear","kiloyear^-1","kiloyear\u207B\u00B9"),1/UNITS[[OP]]$YEAR/1000)
  add(c("per my","1/my","my^-1","my\u207B\u00B9","per ma","1/ma","ma^-1","ma\u207B\u00B9","per megaannum","1/megaannum","megaannum^-1","megaannum\u207B\u00B9","per megayear","1/megayear","megayear^-1","megayear\u207B\u00B9","per millionennium","1/millionennium","millionennium^-1","millionennium\u207B\u00B9"),1/UNITS[[OP]]$YEAR/1000^2)
  add(c("per ae","1/ae","ae^-1","ae\u207B\u00B9","per ga","1/ga","ga^-1","ga\u207B\u00B9","per gy","1/gy","gy^-1","gy\u207B\u00B9","per gyr","1/gyr","gyr^-1","gyr\u207B\u00B9","per aeon","1/aeon","aeon^-1","aeon\u207B\u00B9","per eon","1/eon","eon^-1","eon\u207B\u00B9","per gigayear","1/gigayear","gigayear^-1","gigayear\u207B\u00B9","per gigaannum","1/gigaannum","gigaannum^-1","gigaannum\u207B\u00B9"),1/UNITS[[OP]]$YEAR/1000^3)
  
  # mass
  add(c("g","gm","gram","grams"),1/1000) # kg is SI
  add(c("kg","kilogram","kilograms"),1) # kg is SI
  add(c("Mg","t","tonne","tonnes","mt","m ton","m tons","metric ton","metric tons"),1000)
  add(c("mg","milligram","milligrams"),1/1000)
  add(c("\u03BCg","microgram","micrograms"),1/1000^2)
  add(c("ng","nanogram","nanograms"),1/1000^3)
  add(c("lb","lbs","pound","pounds"),0.45359237)
  add(c("oz","ounce","ounces"),0.45359237/16)
  add(c("slug","slugs"),14.59390)
  add(c("st","stone","stones"),0.45359237*14)
  add(c("ton","tons"),0.45359237*2000) # NA ton (not UK)
  
  # memory
  add(c("byte","bytes","B"),1)
  add(c("Kb","KiB"),1024)
  add(c("Mb","MiB"),1024^2)
  add(c("Gb","GiB"),1024^3)
  add(c("Tb","TiB"),1024^4)
  add(c("Pb","PiB"),1024^5)
  add(c("Eb","EiB"),1024^6)
  add(c("Zb","ZiB"),1024^7)
  add(c("Yb","YiB"),1024^8)
  
  return(list(alias=alias,scale=scale))
}
UNIT <- list() # generated onLoad
#UNIT <- generate.units()


dimfig <- function(data,dimension,thresh=1,...)
{
  UNITS <- unit(data,dimension=dimension,thresh=thresh,...)
  data <- data/UNITS$scale
  R <- list(data=data,units=c(UNITS$name,UNITS$abrv))
  return(R)
}

# CHOOSE BEST UNITS FOR A LIST OF DATA
# thresh is threshold for switching to coarser unit
# concise gives abbreviated names
unit <- function(data,dimension,thresh=1,concise=FALSE,SI=FALSE,...)
{
  if(SI) { data <- 1.01 ; thresh <- 1 } # will always choose base units
  OP <- getOption("time.units")
  
  if(dimension %in% c("length",'distance'))
  {
    name.list <- c("microns","milimeters","centimeters","meters","kilometers")
    abrv.list <- c("\u03BCm","mm","cm","m","km")
    scale.list <- c(1E-6,1/1000,1/100,1,1000)
  }
  else if(dimension=="area")
  {
    name.list <- c("square microns","square milimeters","square centimeters","square meters","hectares","square kilometers")
    abrv.list <- c("\u03BCm\u00B2","mm\u00B2","cm\u00B2","m\u00B2","hm\u00B2","km\u00B2")
    scale.list <- c(1E-12,1/1000^2,1/100^2,1,100^2,1000^2)
  }
  else if(dimension=="time")
  {
    name.list <- c("microseconds","miliseconds","seconds","minutes","hours","days","months","years","millennia","mega-anna","aeons")
    abrv.list <- c("\u03BCs","ms","sec","min","hr","day","mon","yr","ka","Ma","AE")
    scale.list <- c(1E-6,1/1000,1,60*c(1,60)) # through minutes
    scale.list[6] <- UNITS[[OP]]$DAY
    scale.list[7] <- UNITS[[OP]]$MONTH
    scale.list[8] <- UNITS[[OP]]$YEAR
    scale.list[9] <- 1000 * UNITS[[OP]]$YEAR
    scale.list[10] <- 1000^2 * UNITS[[OP]]$YEAR
    scale.list[11] <- 1000^3 * UNITS[[OP]]$YEAR
  }
  else if(dimension %in% c("speed",'velocity'))
  {
    name.list <- c("microns/day","milimeters/day","centimeters/day","meters/day","kilometers/day")
    abrv.list <- c("\u03BCm/day","mm/day","cm/day","m/day","km/day")
    scale.list <- c(1E-6,1/1000,1/100,1,1000)/UNITS[[OP]]$DAY
    
    # SI units fix
    if(SI)
    {
      name.list <- "meters/second"
      abrv.list <- "m/s"
      scale.list <- 1
    }
  }
  else if(dimension=="diffusion")
  {
    name.list <- c("square microns/day","square milimeters/day","square centimeters/day","square meters/day","hectares/day","square kilometers/day")
    abrv.list <- c("\u03BCm\u00B2/day","mm\u00B2/day","cm\u00B2/day","m\u00B2/day","hm\u00B2/day","km\u00B2/day")
    scale.list <- c(1E-12,1/1000^2,1/100^2,1,100^2,1000^2)/UNITS[[OP]]$DAY
    
    # SI units fix
    if(SI)
    {
      name.list <- "square meters/second"
      abrv.list <- "m\u00B2/s"
      scale.list <- 1
    }
  }
  else if(dimension=="frequency")
  {
    name.list <- c("per aeon","per mega-anna","per millennia","per year","per month","per day","per hour","per minute","hertz","kilohertz","megahertz","gigahertz","terahertz")
    abrv.list <- c("AE\u207B\u00B9","Ma\u207B\u00B9","ka\u207B\u00B9","yr\u207B\u00B9","mon\u207B\u00B9","day\u207B\u00B9","hr\u207B\u00B9","min\u207B\u00B9","Hz","kHz","MHz","GHz","THz")
    scale.list <- c(1/UNITS[[OP]]$YEAR/1000^3,1/UNITS[[OP]]$YEAR/1000^2,1/UNITS[[OP]]$YEAR/1000,1/UNITS[[OP]]$YEAR,1/UNITS[[OP]]$MONTH,1/UNITS[[OP]]$DAY,1/60^2,1/60,1,1000,1000^3,1000^6,1000^9,1000^12)
  }
  else if(dimension=="mass")
  {
    name.list <- c("nanograms","micrograms","milligrams","grams","kilograms","tonnes")
    abrv.list <- c("ng","\u03BCg","mg","gm","kg","Mg")
    scale.list <- c(1/1000^4,1/1000^3,1/1000^2,1/1000,1,1000)
  }
  else # units not recognized
  {
    R <- list(scale=1,name=NULL)
    return(R)
  }
  
  data <- data[!is.na(data)]
  data <- data[abs(data)<Inf]
  if(length(data)) { max.data <- max(abs(data)) } else { max.data <- 1 }
  
  if(concise) { name.list <- abrv.list }
  
  # choose most parsimonious units
  I <- (max.data >= thresh * scale.list)
  if(any(I))
  {
    I <- (1:length(I))[I]
    I <- last(I)
  }
  else { I <- 1 }
  
  name <- name.list[I]
  abrv <- abrv.list[I]
  scale <- scale.list[I]
  
  return(list(scale=scale,name=name,abrv=abrv))
}


### determine parsimonious units for parameter CIs
# preference point estimate, but fall back on high-CI if point estimate is zero
unit.par <- function(par,...)
{
  PAR <- par[2:3]
  PAR <- PAR[PAR>.Machine$double.eps]
  if(length(PAR)) { PAR <- min(PAR) } else { PAR <- 0 }
  
  return( unit(PAR,...) )
}


## rescale the units of telemetry object
unit.telemetry <- function(data,length=1,time=1,axes=c('x','y'))
{
  if(class(data)[1]=="phylometry")
  {
    lag <- attr(data,"lag")/time
    data[,axes] <- data[,axes]/length
    attr(data,"lag") <- lag
    # class(data) <- "phylometry"
    return(data)
  }
  # TELEMETRY CLASS BELOW
  
  convert <- function(NAMES,scale) { for(NAME in NAMES) { if(NAME %in% names(data)) { data[[NAME]] <<- data[[NAME]]/scale } } }
  
  convert('t',time)
  convert('light.time',time)
  convert('dark.time',time)
  convert('sundial.rate',1/time)
  
  if(any(axes %in% c('x','y','z')))
  {
    convert(DOP.LIST$horizontal$axes,length)
    convert(DOP.LIST$vertical$axes,length)
    convert(DOP.LIST$speed$axes,length/time)
    
    convert(DOP.LIST$horizontal$VAR,length^2)
    convert(DOP.LIST$horizontal$COV,length^2)
    
    convert(DOP.LIST$vertical$VAR,length^2)
    
    convert(DOP.LIST$speed$VAR,(length/time)^2)
    convert(DOP.LIST$speed$COV,(length/time)^2)
  }
  else
  { for(axis in axes) { convert(axis,length) } }
  
  # calibration constants
  attr(data,"UERE")$UERE <- attr(data,"UERE")$UERE/length # don't logicals get divided this way?
  
  return(data)
}


## rescale the units of dimensionful parameters
unit.ctmm <- function(CTMM,length=1,time=1)
{
  if(length(CTMM$tau)){ CTMM$tau <- CTMM$tau/time }
  CTMM$omega <- CTMM$omega * time
  CTMM$circle <- CTMM$circle * time
  
  # all means scale with length the same way... but not time
  if("mu" %in% names(CTMM))
  {
    CTMM$mu <- CTMM$mu/length
    drift <- get(CTMM$mean)
    CTMM <- drift@scale(CTMM,time)
  }
  
  if(!is.null(CTMM$timelink.cycle))
  { CTMM$timelink.cycle <- CTMM$timelink.cycle/time }
  
  # if(class(CTMM$error)[1]=='numeric')
  { CTMM$error <- CTMM$error/length } # don't divide logicals
  
  if("sigma" %in% names(CTMM))
  {
    CTMM$sigma <- scale.covm(CTMM$sigma,1/length^2)
    
    # variance -> diffusion adjustment
    if(!CTMM$range)
    { CTMM$sigma <- scale.covm(CTMM$sigma,time) }
  }
  
  if("COV.mu" %in% names(CTMM)) { CTMM$COV.mu <- CTMM$COV.mu/length^2 }
  
  if("COV" %in% names(CTMM))
  {
    NAMES <- dimnames(CTMM$COV)[[1]]
    
    PAR <- c('major','minor')
    PAR <- PAR[PAR %in% NAMES]
    if(length(PAR))
    {
      CTMM$COV[PAR,] <- CTMM$COV[PAR,]/length^2
      CTMM$COV[,PAR] <- CTMM$COV[,PAR]/length^2
      
      if(!CTMM$range)
      {
        CTMM$COV[PAR,] <- CTMM$COV[PAR,]*time
        CTMM$COV[,PAR] <- CTMM$COV[,PAR]*time
      }
    }
    
    tau <- CTMM$tau
    tau <- tau[tau<Inf]
    if(length(tau))
    {
      tau <- NAMES[grepl("tau",NAMES)]
      if(length(tau))
      {
        CTMM$COV[tau,] <- CTMM$COV[tau,]/time
        CTMM$COV[,tau] <- CTMM$COV[,tau]/time
      }
      
      if("omega" %in% NAMES)
      {
        CTMM$COV["omega",] <- CTMM$COV["omega",] * time
        CTMM$COV[,"omega"] <- CTMM$COV[,"omega"] * time
      }
    }
    
    ERROR <- NAMES[grepl("error",NAMES)] # error estimate covariance
    if(length(ERROR))
    {
      CTMM$COV[ERROR,] <- CTMM$COV[ERROR,]/length
      CTMM$COV[,ERROR] <- CTMM$COV[,ERROR]/length
    }
    
    if("circle" %in% NAMES)
    {
      CTMM$COV["circle",] <- CTMM$COV["circle",] * time
      CTMM$COV[,"circle"] <- CTMM$COV[,"circle"] * time
    }
  }
  
  return(CTMM)
}


######################
unit.UD <- function(UD,length=1)
{
  UD$r <- lapply(UD$r,function(x){ x/length })
  UD$PDF <- UD$PDF * length^2
  UD$dr <- UD$dr / length
  UD$H <- UD$H / length^2
  
  return(UD)
}


##################
unit.variogram <- function(SVF,time=1,area=1)
{
  SVF$lag <- SVF$lag / time
  SVF$SVF <- SVF$SVF / area
  if("MSE" %in% names(SVF)) { SVF$MSE <- SVF$MSE / area }
  
  return(SVF)
}


# convert units
`%#%` <- function(x,y)
{
  # convert to si units
  if(is.numeric(x))
  {
    num <- x
    name <- y
    pow <- +1
  }
  else # convert from si units
  {
    if(!is.numeric(y)) { return(x %#% 1 %#% y) }
    
    num <- y
    name <- x
    pow <- -1
  }
  
  name <- canonical.name(name)
  if(name=="") { return(num) }
  
  name <- strsplit(name,'*',fixed=TRUE)[[1]]
  if(length(name)>1)
  {
    if(pow==1)
    { for(i in 1:length(name)) { num <- num %#% name[i] } }
    else if(pow==-1)
    { for(i in 1:length(name)) { num <- name[i] %#% num } }
    return(num)
  }
  
  name <- strsplit(name,"/",fixed=TRUE)[[1]]
  if(length(name)>1)
  {
    if(pow==1)
    {
      num <- num %#% name[1]
      for(i in 2:length(name)) { num <- name[i] %#% num }
    }
    else if(pow==-1)
    {
      num <- name[1] %#% num
      for(i in 2:length(name)) { num <- num %#% name[i] }
    }
    return(num)
  }
  
  alias <- UNIT$alias
  scale <- UNIT$scale
  for(i in 1:length(alias))
  {
    if(name %in% alias[[i]]) { return(num*scale[i]^pow) }
  }
  stop(paste("Unit",name,"unknown."))
}

# interpret a string as number with units
ustring <- function(x)
{
  x <- canonical.name(x)
  
  y <- strsplit(x,"[a-z,A-Z]")[[1]][1]
  n <- nchar(y)+1
  y <- as.numeric(y)
  x <- substr(x,n,nchar(x))
  
  x <- y %#% x
  return(x)
}

