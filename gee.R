source("~/Documents/work/code/CSG/CORA/CORA_output.R")
library(rgee)
library(data.table)
library(ggplot2)
library(ggpattern)
library(units)
library(dplyr)

ee$Initialize(project='ee-gorzo')
srtm <- ee$Image("USGS/NLCD/NLCD2016")
landcover <- ee$Image(srtm$select("landcover"))
landcoverVis <- list(
  min = 0.0,
  max = 95.0,
  palette = c(
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "466b9f",
    "d1def8",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "dec5c5",
    "d99282",
    "eb0000",
    "ab0000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "b3ac9f",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "68ab5f",
    "1c5f2c",
    "b5c58f",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "af963c",
    "ccb879",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "dfdfc2",
    "d1d182",
    "a3cc51",
    "82ba9e",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "dcd939",
    "ab6c28",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "000000",
    "b8d9eb",
    "000000",
    "000000",
    "000000",
    "000000",
    "6c9fb8"
  )
)
Map$setCenter(-95, 38, 5)
Map$addLayer(landcover, landcoverVis, "Landcover")

#this comes from CORA_output.R
polygons <- sf_as_ee(fc)

geometry_rect <- ee$Geometry$Rectangle(
  coords = c(-122, 32, -113, 40),
  proj = "EPSG:4326",
  geodesic = FALSE
)

#Map$addLayer(geometry_rect)
#nlcd <- landcover$reproject(sa_extent$projection())

stats <- landcover$reduceRegions(
  collection = polygons,
  reducer = ee$Reducer$frequencyHistogram()$combine(
    reducer2 = ee$Reducer$count(),
    sharedInputs = TRUE
  )
  #geometry = geometry_rect,
  #scale = 10,
  #bestEffort = TRUE # Use maxPixels if you care about scale.
)

myshp <- ee_as_sf(stats)
dictcount <- lapply(myshp$histogram,jsonlite::fromJSON)
#myshp$vals=dictcount

#task <- ee$batch$Export$table$toDrive(
#  collection = ee$FeatureCollection(ee$Feature(NULL, stats)),
#  description = "exported_stats_demo",
#  fileFormat = "CSV"
#)
#task$start()
#ee_monitoring(task)
#ee_extract()
#exported_stats <- ee_drive_to_local(task = task,dsn = "exported_stats.csv")
#read.csv(exported_stats)
#pop <- ee_as_rast(image=landcover, region=geometry, via="drive",  maxPixels = 1150000000)
#load your workspace that has hr in it

#write code to do raster extract to get vector of cell values per shapefile
lc <- function(x) {
  #sum(unname(unlist(x[names(x) %in% c("21", "22")])))
  ag <- sum(unname(unlist(x[names(x) %in% c("81","82")])))
  develop <- sum(unname(unlist(x[names(x) %in% c("21","22", "23", "24")])))
  veg <- sum(unname(unlist(x[names(x) %in% c("41","42", "43", "51", "52", "71", "72", "73", "74", "90", "95")])))
return(data.frame(agriculture=ag,developed=develop,native_veg=veg))}

lulc <- rbindlist(lapply(dictcount, lc))
myshp2 <- cbind(myshp,lulc)
myshp2$ag_percent <- myshp2$agriculture / myshp2$count 
myshp2$develop_percent <- myshp2$developed / myshp2$count 
myshp2$veg_percent <- myshp2$native_veg / myshp2$count 
raven_out <- data.frame(ID=myshp2$id, area_m2 = myshp2$area, ag_percent = myshp2$ag_percent, develop_percent = myshp2$develop_percent, veg_percent = myshp2$veg_percent)
raven_out$bird <- sapply(strsplit(raven_out$ID, "_"), "[[", 1)
raven_out$month <- as.integer(sapply(strsplit(raven_out$ID, "_"), "[[", 2))
raven_out$year <- as.integer(sapply(strsplit(raven_out$ID, "_"), "[[", 3))
raven_out$area_km2 <- raven_out$area_m2 / 1000000
group_mean<- aggregate(x= raven_out$area_km2,
                       # Specify group indicator
                       by = list(raven_out$month, raven_out$year),      
                       # Specify function (i.e. mean)
                       FUN = function(x) { c(MEAN=mean(x), min=min(x), max=max(x))})
group_sd<- aggregate(x= raven_out$area_km2,
                     # Specify group indicator
                     by = list(raven_out$month, raven_out$year),      
                     # Specify function (i.e. mean)
                     FUN = sd)
ag_mean<- aggregate(x= list(raven_out$ag_percent, raven_out$develop_percent, raven_out$veg_percent),
                     # Specify group indicator
                     by = list(raven_out$month, raven_out$year),      
                     # Specify function (i.e. mean)
                     FUN = mean)
group_mean$sd <- group_sd$x
summary_stats_mcp <- data.frame(Month=as.integer(group_mean[,1]), Year=as.integer(group_mean[,2]), Area_Mean_km2 = group_mean$x[,1], Area_Min = group_mean$x[,2], Area_Max = group_mean$x[,3], Area_SD=group_mean[,4], Ag_Mean = ag_mean[,3], Develop_Mean = ag_mean[,4], Veg_Mean = ag_mean[,5])
#summary_stats$YearMo <- paste(summary_stats$Year, summary_stats$Month, sep="-")
summary_stats_mcp <- summary_stats_mcp[order(summary_stats_mcp$Year,summary_stats_mcp$Month),]
summary_stats_land <- summary_stats_mcp[,c("Month", "Year",  "Ag_Mean", "Develop_Mean", "Veg_Mean")]
summary_stats_land$YearMo <- paste(summary_stats$Year, summary_stats$Month, sep="-")
summary_stats_land$Month <- NULL
summary_stats_land$Year <- NULL
long <- reshape(summary_stats_land, varying = c("Ag_Mean", "Develop_Mean", "Veg_Mean"), v.names = "val", timevar = "LULC", times = c("Ag_Mean", "Develop_Mean", "Veg_Mean"), direction="long")
ggplot(long,aes(x=YearMo, y=as.numeric(val), group=LULC)) +
  geom_line(aes(color=LULC)) +
  ylab("Proportion of Home Range") + xlab("Month")#+ scale_y_continuous(limits = c(-0.5, NA)) + 
#  scale_y_continuous(trans='log10') +
#  theme(axis.text = element_text(size = 12, color="black"), strip.text = element_text(size = 12), axis.title = element_text(size = 12), legend.position = 'none', axis.line=element_line(colour="black"), axis.title.x=element_blank())

#write.csv(raven_out, "raven_lulc.csv")

#res = monthglm(area_km2~ag_percent + develop_percent + veg_percent, monthvar="month", data=raven_out)
library(season)
res = cosinor(area_km2~ag_percent + develop_percent + veg_percent, date="month", type="monthly", data=raven_out)
#summary(res)
library(lme4)
library(lmerTest)
forestmodel <- lmer(area_km2 ~ ag_percent + develop_percent + veg_percent + as.character(month) + (1|bird) , data=raven_out) #+ (1 + ag_percent| month) 
require(lattice)
dotplot(ranef(forestmodel))

#write.csv(summary_stats_mcp, "~/Documents/CORA/summary_stats.csv")


