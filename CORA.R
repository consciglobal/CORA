library(ggplot2) # for plotting
library(dplyr)
library(zoo)
library(amt) # for movement and overlap
library(sf)
library(furrr)
library(purrr)
library(ctmm)
cora <- read.csv("~/Documents/CSG/CORA/data/CORA_Raw.csv")

#no HDOP > 10 in data set
cora <- cora[cora$Fix > 2,]
cora <- cora[cora$Altitude < 4000,]
cora <- cora[cora$Longitude < -80,]

#look at a single bird to try to see seasonal clusters
CA     <- cora[cora$Animal_ID==103,]
CA$datetime <- as.POSIXct(CA$Date_Time, format="%m/%d/%Y %H:%M:%S", tz="UTC")
CA$dt <- as.integer(CA$datetime) 
df     <- data.frame(long=CA$Longitude, lat=CA$Latitude, city=CA$datetime, month=CA$Month_LocalTime, dt=CA$dt, day=CA$JulianDate_LocalTime)
geo.dist = function(df) {
  require(geosphere)
  d <- function(i,z){         # z[1:2] contain long, lat
    dist <- rep(0,nrow(z))
    dist[i:nrow(z)] <- distHaversine(z[i:nrow(z),1:2],z[i,1:2])
    return(dist)
  }
  dm <- do.call(cbind,lapply(1:nrow(df),d,df))
  return(as.dist(dm))
}
d      <- geo.dist(df)   # distance matrix
hc     <- hclust(d)      # hierarchical clustering
plot(hc)
df$clust <- cutree(hc,k=2)

df1 <- df[df$month == 12,]
ggplot(df1)+
  geom_point(aes(x=long, y=lat, color=factor(clust),shape=as.factor(month), alpha=day/max(day)), size=4) +
  scale_color_discrete("Cluster")+
  coord_fixed()

#subset the data to play with/look at speed and direction for identifying bad points
#cora2 <- cora %>% mutate(i = rollapply(KPH > 120, 2, all, align = "left", fill = FALSE))
#speed <- sort(c(sapply(1:5,function(x) c(which(cora2$i)-x, (which(cora2$i)+x))), which(cora2$i))) #look at points around high speed
#should it be the point that has a high speed, and then a high speed following it, where the subsequent point also has a direction change?
#possibly also that point has a direction change?
#cora2$lat <- c(0,diff(cora2$Latitude))
#cora2$lon <- c(0,diff(cora2$Longitude))
#which row has a change in sign in the difference from the row before it (i.e., where does lat/long change from increasing/decreasing?)
#gives you the row where the direction change occurs (e.g. lat or long was decreasing, until this row which indicates an increase)
#dirchange <- unique(c(which(c(0,diff(sign(cora2$lat)))!=0), which(c(0,diff(sign(cora2$lon)))!=0))) 
#cora2 <- cora2[-speed,]
#could possibly integrate heading and the calculated distance from last point here too
#####

cora$animalmo <- paste(cora$Animal_ID, cora$Month_LocalTime, cora$Year_LocalTime, sep="_")

#check to see number of records per day per bird
cora$animalday <- paste(cora$Animal_ID, cora$Date)
table(cora$animalday)
group_len <- aggregate(x= cora$Latitude,
                       # Specify group indicator
                       by = list(cora$Date, cora$Animal_ID),      
                       # Specify function (i.e. mean)
                       FUN = length)
names(group_len) <- c("Date", "Bird", "Records")

turtles_track <- cora %>%
  filter(!is.na(Latitude), !is.na(Longitude)) %>%
  mutate(ts = as.POSIXct(Date_Time, format="%m/%d/%Y %H:%M:%S", tz="GMT")) %>%
  make_track(Longitude, Latitude, ts, id = animalmo, HDOP=HDOP, VDOP=VDOP, crs = 4326) %>%
  # Use nest() to allow us to deal with multiple animals (5 in sample set)
  nest(data = -"id") %>%
  arrange(id)

posslm1 = possibly(.f = hr_akde, otherwise = "Error")
#posslm2 = possibly(.f = hr_isopleths, otherwise = "Error")

plan(multisession, workers=27)
startTime <- Sys.time()
hr <- turtles_track %>% 
  mutate(hrmcp = map(data, function(x) {tryCatch({hr_mcp(x, levels = c(1.0))}, error= function(cond){st_sf(st_sfc())})})) %>%
         mutate(hr_akde = future_map(data, ~posslm1(.x, fit_ctmm(.x, "auto")))) #I think it has to be in its own mutate() call to run parallel
         
#hr <- hr%>%
   # mutate(isopleth = map(hr_akde, ~posslm2(.x)))

#TOO RESOURCE INTENSIVE! What's going on here?
#hr <- hr %>% mutate(
#  hr_akde_error = future_map(data, ~posslm1(.x, fit_ctmm(.x, "auto", uere = 1.67))))

endTime <- Sys.time() 
print(endTime - startTime)
#plot(hr[10,]$hrmcp[[1]])
#test <- plot(hr[1,]$isopleth[[1]][3])
#save(hr, file="hr3.RData")

