rm(list = ls())
library(openCR)

source(file.path("./fitOpenCR_FUN.R"))

load(file.path("./ModelData/capt_pop.rda"))
load(file.path("./ModelData/LTpoints.rda"))


#####Fit the models. Don't need to re-run if already done.######
LTest <- fitOpenCR(data.dir = file.path("./ModelData"),
          results.dir = file.path("./ModelResults/"),
          searchPoints=LTpoints,
          capt_pop=capt_pop,
          nocc=3,
          buffer=1500,
          site="LakeTyers",
          nProcessors=16)


#source process model results script, this will load model fits that are already created in the results dir
results.dir = file.path("./ModelResults", "LakeTyers")
source(file.path("./process.modelRes.R"))

#Examine AIC table and decide which models are identifiable and have substantial support.


#Estimate derived parameters and tabulate model results
res1 <- process.modelRes(MODEL=JS_CL_HHR, site="LakeTyers", 
                         results.dir = results.dir) 

res2 <- process.modelRes(MODEL=JS_CL_BVN_HHR, site="LakeTyers", 
                        results.dir = results.dir) 

res3 <- process.modelRes(MODEL=JS_CL_HHR_t, site="LakeTyers", 
                        results.dir = results.dir) 

res4 <- process.modelRes(MODEL=JS_CL, site="LakeTyers", 
                         results.dir = results.dir) 



##Polygon area size over which density was estimated.
load(file = file.path(results.dir, paste0("capt_hist_", "LakeTyers",".rda")))
load(file = file.path(results.dir, paste0("mask_", "LakeTyers", ".rda")))
sf_mask <- st_as_sf(mask, coords=c("x", "y"), crs=32755)
mask_as_grid <- st_make_grid(sf_mask, cellsize = 100, what="polygons", square=FALSE, offset = c(0,0))
sum(st_area(mask_as_grid))/1e6


