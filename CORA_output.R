library(purrr)
library(amt)
load("~/Documents/work/data/CSG/CORA/hr.RData") #UNCOMMENT TO GET THIS TO WORK
#allpts <- rbindlist(h1$data)
#bound <- c(min(allpts$x_), min(allpts$y_), max(allpts$x_), max(allpts$y_))
posslm2 = possibly(.f = hr_isopleths, otherwise = "Error")
hr <- hr%>%
  mutate(isopleth = map(hr_akde, ~posslm2(.x)))

mcps <- purrr::map(hr$hrmcp, 1)

test <- sapply(mcps, is.data.frame)
mcps_good <- mcps[test]
ids <- hr$id[test]
ided <- Map(cbind, mcps_good, id = ids)
fc <- bind_rows(ided)

good_iso <- sapply(hr$isopleth, is.data.frame)
akde_iso <- hr$isopleth[good_iso]
ids <- hr$id[good_iso]
ided <- Map(cbind, akde_iso, id = ids)
akde_cora <- bind_rows(ided)
akde_cora_95 <- akde_cora[akde_cora$what=="estimate",]
raven <- data.frame(ID=akde_cora_95$id, area = akde_cora_95$area)
units(raven$area) = 'km2'
raven$log_area <- log(raven$area)

raven$bird <- sapply(strsplit(raven$ID, "_"), "[[", 1)
raven$mo <- as.factor(as.integer(sapply(strsplit(raven$ID, "_"), "[[", 2)))
raven$year <- sapply(strsplit(raven$ID, "_"), "[[", 3)
raven$moyear <- paste(raven$mo, raven$year, sep="_")
#raven_ord <- raven[order(raven$area, decreasing=TRUE),]

#ggplot(raven,aes(x=mo, y=as.numeric(area))) +
#  geom_boxplot_pattern() + facet_grid(cols = vars(year)) +
#  ylab("Area") + #+ scale_y_continuous(limits = c(-0.5, NA)) + 
#  scale_y_continuous(trans='log10') +
#  theme(axis.text = element_text(size = 12, color="black"), strip.text = element_text(size = 12), axis.title = element_text(size = 12), legend.position = 'none', axis.line=element_line(colour="black"), axis.title.x=element_blank())

group_mean<- aggregate(x= raven$area,
                       # Specify group indicator
                       by = list(raven$mo, raven$year),      
                       # Specify function (i.e. mean)
                       FUN = function(x) { c(MEAN=mean(x), min=min(x), max=max(x))})
group_sd<- aggregate(x= raven$area,
                     # Specify group indicator
                     by = list(raven$mo, raven$year),      
                     # Specify function (i.e. mean)
                     FUN = sd)
group_mean$sd <- group_sd$x
summary_stats <- data.frame(Month=group_mean[,1], Year=group_mean[,2], Area_Mean = group_mean$x[,1], Area_Min = group_mean$x[,2], Area_Max = group_mean$x[,3], Area_SD=group_mean[,4])
#names(group_mean) <- c("Month", "Year", "Area_Mean", "Area_Min", "Area_Max", "Area_SD")
# Print the resulting summary data frame
#print(group_mean)
group_mean %>%
  ggplot(aes(x=Month, y=Area, group=Year, color=Year)) +
  geom_line()