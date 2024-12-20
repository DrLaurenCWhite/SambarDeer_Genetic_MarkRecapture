modelList <- c("JS_CL", "JS_CL_Ind", "JS_CL_BVN",
               "JS_CL_HHR", "JS_CL_Ind_HHR", "JS_CL_BVN_HHR",
               "JS_CL_t", "JS_CL_Ind_t", "JS_CL_BVN_t", 
               "JS_CL_HHR_t", "JS_CL_Ind_HHR_t", "JS_CL_BVN_HHR_t")



# Load models in memory
#modelInMemory <- sapply(modelList, exists)
#if(sum(!modelInMemory)>0) {
  for(model in paste0(modelList, ".rda")) {
    load(file = file.path(results.dir, model))
  }
#}


process.modelRes <- function(MODEL, site, results.dir=results.dir, nProcessors=6) {
  # checked whether there are problem of identifiability (no variance estimates or
  # rank of Hessian matrix < # param)
  
  isModelOK <- length(MODEL$betanames) == length(MODEL$betanames) -
    length(which(MODEL$fit$hessian[1, ] == 0)) 
  if(isModelOK) message("No problem detected with the model") else
    warning("Check the model there is a problem!")
  
  # Set movement model
  mvm <- MODEL$movementmodel
  
  # Set the detection fn
  detf <- openCR:::.openCRstuff$DFN[MODEL$detectfn - 13]
  
  # Get the model formula. There must be a more elegant way, but this is what I got at this stage
  ac<-AIC(MODEL)
  mod <- lapply(strsplit(ac$model, split = " ")[[1]], formula)
  
  # This gets derived param (Density) but without ci
  der <- derived(MODEL, all.levels = T)
  
  print(der)
  
  # CI for primary params
  pred <- predict(MODEL, type="response")
  
  
  # confidence interval of derived params
  load(file = file.path(results.dir, paste0("capt_hist_", site,".rda")))
  load(file = file.path(results.dir, paste0("mask_", site, ".rda")))
  betas <- MODEL$fit$estimate
  names(betas) <- MODEL$betanames
  
  
  # Density goes first
  logD <- log(der$estimates$D)
  betaD <- logD
  for(i in 2:length(logD)) {
    betaD[i] <- logD[i]- logD[1]
  }
  names(betaD) <- c("D", paste0("D.t", 2:length(sessionlabels(CH))))
  
  Dvar <- openCR.fit(CH, type = "JSSAsecrD", mask=mask, method = "none",
                     model=c(mod, formula(D ~ t)),
                     movementmodel = mvm,
                     detectfn = detf, start= c(betas, betaD),
                     ncores = nProcessors)
  
  
    
  # get params on the right scale
  predD <- predict(Dvar, type="response")
  predD
  pred[["D"]] <- predD[["D"]]
  
  # lambda
  betaLam <- log(der$estimates$lambda)
  names(betaLam) <- "lambda"
  Lamvar <- openCR.fit(CH, type = "JSSAsecrlCL", mask=mask, method = "none",
                     model=c(mod, formula(lambda ~ 1)),
                     movementmodel = mvm,
                     detectfn = detf, start= c(betas, betaLam[1]),
                     ncores = nProcessors)
  
  predLam <- predict(Lamvar, type="response")
  pred[["lambda"]] <- predLam[["lambda"]]
  
  # Put results together
  res <- data.table::rbindlist(pred, idcol = "Param")
  # Density in sqkm
  res[Param == "D", c("estimate", "SE.estimate", "lcl", "ucl") := lapply(.SD, "*", 100), 
      .SDcols = c("estimate", "SE.estimate", "lcl", "ucl")]
  # Work out what rows to keep
  rkeep <- apply(res[, c("estimate", "SE.estimate", "lcl", "ucl"), with=F], 1, 
                 FUN = function(x) sum(is.na(x)))
  res <- res[rkeep != 4, ]
  
  # Round to second decimal dig
  res[, c("estimate", "SE.estimate", "lcl", "ucl") := lapply(.SD, round, 2),
      .SDcols = c("estimate", "SE.estimate", "lcl", "ucl")]
  write.csv(res, file = file.path(results.dir, paste(
    paste("Par_est_CI", site, mvm, detf, 
          if(as.character(MODEL$model$lambda0)[2] == "t") "t", sep="_"),
    "csv", sep=".")), row.names = F)
  print(res)
  return(res)
}