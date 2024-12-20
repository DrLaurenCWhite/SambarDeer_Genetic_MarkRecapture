rm(list = ls())
library(openCR)

source(file.path("./fitOpenCR_FUN.R"))

load(file.path("./ModelData/capt_pop.rda"))
load(file.path("./ModelData/SR2points.rda"))


#####Fit the models. Don't need to re-run if already done.######
SRest <- fitOpenCR(data.dir = file.path("./ModelData"),
          results.dir = file.path("./ModelResults/"),
          searchPoints=SR2points,
          capt_pop=capt_pop,
          nocc=3,
          buffer=1500,
          site="SnowyRiver",
          nProcessors=16)

#source process model results script, this will load model fits that are already created in the results dir
results.dir = file.path("../ModelResults", "SnowyRiver")
source(file.path("./process.modelRes.R"))

#Examine AIC table and decide which models are identifiable and have substantial support.

#Estimate derived parameters and tabulate model results
res <- process.modelRes(MODEL=JS_CL_BVN_HHR_t, site="SnowyRiver", 
                        results.dir = results.dir, nProcessors = 7)

##Polygon area size over which density was estimated.
load(file = file.path(results.dir, paste0("capt_hist_", "SnowyRiver",".rda")))
load(file = file.path(results.dir, paste0("mask_", "SnowyRiver", ".rda")))
sf_mask <- st_as_sf(mask, coords=c("x", "y"), crs=32755)
mask_as_grid <- st_make_grid(sf_mask, cellsize = 100, what="polygons", square=FALSE, offset = c(0,0))
sum(st_area(mask_as_grid))/1e6






















