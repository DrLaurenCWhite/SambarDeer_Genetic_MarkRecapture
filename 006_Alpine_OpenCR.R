rm(list = ls())
library(openCR)

source(file.path("./fitOpenCR_FUN.R"))

load(file.path("./ModelData/capt_pop.rda"))
load(file.path("./ModelData/ANP2points.rda"))

#####Fit the models. Don't need to re-run if already done.######
ANPest <- fitOpenCR(data.dir = file.path("./ModelData"),
                    results.dir = file.path("./ModelResults/"),
                    searchPoints=ANP2points,
                    capt_pop=capt_pop,
                    nocc=2,
                    buffer=1500,
                    site="Alpine",
                    nProcessors=16)


#source process model results script, this will load model fits that are already created in the results dir
results.dir = file.path("./ModelResults", "Alpine")
source(file.path("./process.modelRes.R"))

#Examine AIC table and decide which models are identifiable and have substantial support.


#Estimate derived parameters and tabulate model results
res <- process.modelRes(MODEL=JS_CL_HHR, site="Alpine", 
                         results.dir =results.dir, nProcessors=16)

res2 <- process.modelRes(MODEL=JS_CL_HHR_t, site="Alpine", 
                         results.dir =results.dir, nProcessors=16)


##Polygon area size over which density was estimated.
load(file = file.path(results.dir, paste0("capt_hist_", "Alpine",".rda")))
load(file = file.path(results.dir, paste0("mask_", "Alpine", ".rda")))
sf_mask <- st_as_sf(mask, coords=c("x", "y"), crs=32755)
mask_as_grid <- st_make_grid(sf_mask, cellsize = 100, what="polygons", square=FALSE, offset = c(0,0))
sum(st_area(mask_as_grid))/1e6
