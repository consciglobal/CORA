library(dplyr)
library(zoo)
library(ggplot2) # for plotting
library(amt) # for movement and overlap
library(sf)
#library(furrr)
library(purrr)
load("~/Documents/CSG/CORA/data/hr.RData")

startTime <- Sys.time()
posslm1 = possibly(.f = hr_akde, otherwise = "Error")
#plan(multisession, workers=31)
hr <- hr %>% mutate(
  hr_akde_error = map(data, ~posslm1(.x, fit_ctmm(.x, "auto", uere = 1.67))))
  #isopleth = map(hr_akde, ~posslm2(.x)))
endTime <- Sys.time() 
print(endTime - startTime)

save(hr, file="hr3.RData")

#plot(hr[1,]$isopleth[[1]][3])
