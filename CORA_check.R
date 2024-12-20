library(purrr)
library(amt)
cora <- read.csv("~/Documents/work/data/CSG/CORA/CORA_Raw.csv")

#no HDOP > 10 in data set
cora <- cora[cora$Fix > 2,]
cora <- cora[cora$Altitude < 4000,]
cora <- cora[cora$Longitude < -80,]
cora$animalmo <- paste(cora$Animal_ID, cora$Month_LocalTime, cora$Year_LocalTime, sep="_")
load("~/Documents/work/data/CSG/CORA/hr.RData")

posslm2 = possibly(.f = hr_isopleths, otherwise = "Error")
hr <- hr%>%
  mutate(isopleth = map(hr_akde, ~posslm2(.x)))

mcps <- purrr::map(hr$hrmcp, 1)

test <- sapply(mcps, is.data.frame)
mcps_bad <- hr[which(!test),]
