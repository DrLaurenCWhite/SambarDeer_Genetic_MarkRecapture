library(tidyverse)
library(readxl)
library(secr)
library(sf)
library(terra)
library(openCR)
library(data.table)
#library(nimble)

source(file.path("./Helper.FUNs.R"))

fitOpenCR <- function(data.dir = file.path("../ModelData"),
                      results.dir = file.path("../ModelResults"),
                      searchPoints=SR2points,
                      capt_pop=capt_pop,
                      nocc=3,
                      buffer=1500,
                      site="SnowyRiver",
                      nProcessors=6) {
   results.dir <- file.path(results.dir, site)
  dir.create(results.dir)

  GBForpoints<- as.data.frame(searchPoints@coords)
  GBpoints<- st_as_sf(GBForpoints, coords=c("Longitude", "Latitude"), crs=32755)

  sspace<- st_convex_hull(st_union(GBpoints))
  sspace<- st_buffer(sspace, buffer)
  GBgrid<- make_grid(sspace, cell_diameter = 100, clip=TRUE)

  # Grab capture data
  capGB <- capt_pop[[site]]
  capGB <- do.call(rbind, capGB)
  if(nocc==2) {
    capGB[, `5`:=NULL]
    capGB[, `6`:=NULL]
  }
  
  capGB[, Treatment := factor(Treatment, levels = c("Pre_Control","Post_Control"))]
  capGB<- st_as_sf(capGB, coords=c("Long","Lat"), crs=4326)

  capGB<- st_transform(capGB, crs=32755)

  # Construct detection matrix
  cap_tr<- st_intersection(capGB, GBgrid)
  cap_tr<- cap_tr %>% arrange(IndID)
  eff_tr<- st_intersection(GBpoints, GBgrid)
  eff_tr<- eff_tr %>% group_by(HexID) %>% summarise(n=n()) %>% st_centroid()

  if(nocc == 3) {
  caphist<- cap_tr %>% select(Treatment=Treatment, ID=IndID, X1,X2,X3,X4,X5,X6, Detector=HexID) %>% st_drop_geometry()
  } else {
    caphist<- cap_tr %>% select(Treatment=Treatment, ID=IndID, X1,X2,X3,X4, Detector=HexID) %>% st_drop_geometry()
  }
  caphist<- caphist %>% pivot_longer(cols = starts_with("X"), names_to = "Occ", values_to = "cap")
  if(nocc == 3) {
    caphist<- caphist %>% mutate(Occasion = case_when(Occ == "X1" ~ 1,
                                                      Occ == "X2" ~ 2,
                                                      Occ == "X3" ~ 3,
                                                      Occ == "X4" ~ 1,
                                                      Occ == "X5" ~ 2,
                                                      Occ == "X6" ~ 3))
  } else {
    caphist<- caphist %>% mutate(Occasion = case_when(Occ == "X1" ~ 1,
                                                      Occ == "X2" ~ 2,
                                                      Occ == "X3" ~ 1,
                                                      Occ == "X4" ~ 2))
  }

  caphist<- caphist %>% filter(cap > 0) %>% select(-Occ, -cap)
  #caphist$Session <- droplevels(caphist$Session, exclude = "Autumn_2017")
  hex_samp<- sort(unique(c(eff_tr$HexID, cap_tr$HexID)))
  traps<- GBgrid %>% filter(HexID %in% hex_samp)
  coords<- st_centroid(traps) %>% st_coordinates() %>% as.data.frame()
  traps<- traps %>% select(Detector = HexID) %>% mutate(x = coords$X, y = coords$Y)
  traps<- traps %>% st_drop_geometry()

  caphist<- caphist %>% mutate(trap = match(Detector, traps$Detector)) %>% select(-Detector)
  searched_cells<- read.traps(data=traps, detector="count")

  mask<- make.mask(searched_cells, type="polygon", poly=sspace, spacing = 100)

  CH<- make.capthist(as.data.frame(caphist), searched_cells, fmt="trapID", noccasions=nocc)

  CH <- join(CH)
  summary(CH)

  # Check intervals are correctly set
  attr(CH, "intervals")

  #win.graph(10,10)
  png(filename = file.path(results.dir,paste0(site, ".png")),
      width = 18, height = 15, units = "cm", res = 300)
  plot(CH, track=TRUE)
  dev.off()
  message(paste("Saved plot of trap history in", results.dir))

  #### Fit models ####
  #### Conditional likelihood ####
  # CJS_CL<- openCR.fit(CH, type = "CJSsecr", mask=mask, model=list(D ~ 1), ncores = 4)

  # save each model as it finishes so that if it crashes we can resume from where it left off.
  # This could be done more elegantly with tryCatch but I can't be bother to set it up (and I have often failed in the past)

  if(!file.exists(file.path(results.dir, "JS_CL.rda"))) {
    message("Started fitting the model JS_CL")
    JS_CL<- openCR.fit(CH, type = "JSSAsecrfCL", mask=mask, ncores = nProcessors)
    save(JS_CL, file=file.path(results.dir, "JS_CL.rda"))
  }

  if(!file.exists(file.path(results.dir, "JS_CL_Ind.rda"))) {
    message("Started fitting the model JS_CL_Ind")
  JS_CL_Ind <- openCR.fit(CH, type = "JSSAsecrfCL", mask=mask,
                                       movementmodel = "IND",
                                       # detectfn = "HHR",
                                       ncores = nProcessors)
  save(JS_CL_Ind, file=file.path(results.dir, "JS_CL_Ind.rda"))
  }

  if(!file.exists(file.path(results.dir, "JS_CL_BVN.rda"))) {
    message("Started fitting the model JS_CL_BVN")
  JS_CL_BVN <- openCR.fit(CH, type = "JSSAsecrfCL", mask=mask,
                                       movementmodel = "BVN",
                                       #detectfn = "HHR",
                                       ncores = nProcessors)
  save(JS_CL_BVN, file=file.path(results.dir, "JS_CL_BVN.rda"))
  }
  
  if(!file.exists(file.path(results.dir, "JS_CL_HHR.rda"))) {
    message("Started fitting the model JS_CL_HHR")
    JS_CL_HHR <- openCR.fit(CH, type = "JSSAsecrfCL", mask=mask,
                                #movementmodel = "IND",
                                detectfn = "HHR",
                                ncores = nProcessors)
    save(JS_CL_HHR, file=file.path(results.dir, "JS_CL_HHR.rda"))
  }

  if(!file.exists(file.path(results.dir, "JS_CL_Ind_HHR.rda"))) {
    message("Started fitting the model JS_CL_Ind_HHR")
  JS_CL_Ind_HHR <- openCR.fit(CH, type = "JSSAsecrfCL", mask=mask,
                                           movementmodel = "IND",
                                           detectfn = "HHR",
                                           ncores = nProcessors)
  save(JS_CL_Ind_HHR, file=file.path(results.dir, "JS_CL_Ind_HHR.rda"))
  }

  if(!file.exists(file.path(results.dir, "JS_CL_BVN_HHR.rda"))) {
       message("Started fitting the model JS_CL_BVN_HHR")
  JS_CL_BVN_HHR <- openCR.fit(CH, type = "JSSAsecrfCL", mask=mask,
                                          movementmodel = "BVN",
                                           detectfn = "HHR",
                                           ncores = nProcessors)
  save(JS_CL_BVN_HHR, file=file.path(results.dir, "JS_CL_BVN_HHR.rda"))
  }

  if(!file.exists(file.path(results.dir, "JS_CL_t.rda"))) {
       message("Started fitting the model JS_CL_t")
    JS_CL_t <- openCR.fit(CH, type = "JSSAsecrfCL", mask=mask,
                                               model=list(lambda0 ~ t),
                                               #movementmodel = "BVN",
                                               #detectfn = "HHR",
                                               ncores = nProcessors)
  save(JS_CL_t, file=file.path(results.dir, "JS_CL_t.rda"))
  }

  if(!file.exists(file.path(results.dir, "JS_CL_Ind_t.rda"))) {
       message("Started fitting the model JS_CL_Ind_t")
    JS_CL_Ind_t <- openCR.fit(CH, type = "JSSAsecrfCL", mask=mask,
                                             model=list(lambda0 ~ t),
                                             movementmodel = "IND",
                                             #detectfn = "HHR",
                                             ncores = nProcessors)
  save(JS_CL_Ind_t, file=file.path(results.dir, "JS_CL_Ind_t.rda"))
  }
  
  if(!file.exists(file.path(results.dir, "JS_CL_BVN_t.rda"))) {
    message("Started fitting the model JS_CL_BVN_t")
    JS_CL_BVN_t <- openCR.fit(CH, type = "JSSAsecrfCL", mask=mask,
                              model=list(lambda0 ~ t),
                              movementmodel = "BVN",
                              #detectfn = "HHR",
                              ncores = nProcessors)
    save(JS_CL_BVN_t, file=file.path(results.dir, "JS_CL_BVN_t.rda"))
  }
  
  if(!file.exists(file.path(results.dir, "JS_CL_HHR_t.rda"))) {
    message("Started fitting the model JS_CL_HHR_t")
    JS_CL_HHR_t <- openCR.fit(CH, type = "JSSAsecrfCL", mask=mask,
                                  model=list(lambda0 ~ t),
                                  #movementmodel = "IND",
                                  detectfn = "HHR",
                                  ncores = nProcessors)
    save(JS_CL_HHR_t, file=file.path(results.dir, "JS_CL_HHR_t.rda"))
  }

  if(!file.exists(file.path(results.dir, "JS_CL_Ind_HHR_t.rda"))) {
    message("Started fitting the model JS_CL_Ind_HHR_t")
    JS_CL_Ind_HHR_t <- openCR.fit(CH, type = "JSSAsecrfCL", mask=mask,
                                  model=list(lambda0 ~ t),
                                  movementmodel = "IND",
                                  detectfn = "HHR",
                                  ncores = nProcessors)
    save(JS_CL_Ind_HHR_t, file=file.path(results.dir, "JS_CL_Ind_HHR_t.rda"))
  }
  
  if(!file.exists(file.path(results.dir, "JS_CL_BVN_HHR_t.rda"))) {
       message("Started fitting the model JS_CL_BVN_HHR_t")
    JS_CL_BVN_HHR_t <- openCR.fit(CH, type = "JSSAsecrfCL", mask=mask,
                                                 model=list(lambda0 ~ t),
                                                 movementmodel = "BVN",
                                                 detectfn = "HHR",
                                                 ncores = nProcessors)
  save(JS_CL_BVN_HHR_t, file=file.path(results.dir, "JS_CL_BVN_HHR_t.rda"))
  }

  modelList <- c("JS_CL", "JS_CL_Ind", "JS_CL_BVN",
    "JS_CL_HHR", "JS_CL_Ind_HHR", "JS_CL_BVN_HHR",
    "JS_CL_t", "JS_CL_Ind_t", "JS_CL_BVN_t", 
    "JS_CL_HHR_t", "JS_CL_Ind_HHR_t", "JS_CL_BVN_HHR_t")


  modelInMemory <- sapply(modelList, exists)
  if(sum(!modelInMemory)>0) {
    for(model in paste0(modelList, ".rda")[!modelInMemory]) {
      load(file = file.path(results.dir, model))
    }
  }

  aic <- AIC(JS_CL, JS_CL_Ind, JS_CL_BVN,
             JS_CL_HHR, JS_CL_Ind_HHR, JS_CL_BVN_HHR,
             JS_CL_t, JS_CL_Ind_t, JS_CL_BVN_t, 
             JS_CL_HHR_t, JS_CL_Ind_HHR_t, JS_CL_BVN_HHR_t)
  write.csv(aic, file.path(results.dir, paste0("AIC_table_", site, ".csv")))
  save(CH, file = file.path(results.dir, paste0("capt_hist_", site, ".rda")))
  save(mask, file = file.path(results.dir, paste0("mask_", site, ".rda")))

  return(list(captureHistory=CH, mask=mask))
}
