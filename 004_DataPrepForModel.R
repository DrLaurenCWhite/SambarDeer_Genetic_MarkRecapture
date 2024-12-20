rm(list = ls())

library(rgeos)
library(sp)
library(data.table)
library(gpx)

data.dir <- "./Data"
output.dir <- "../ModelData"

#### Generate capture histories ####
load(file.path(data.dir, "matches.rda"))
names(matches)
class(matches)

lapply(matches[, c("Pop", "Treatment", "Trip")], unique)
matches[, Pop := gsub(" ", "", Pop)]
matches[, Treatment := gsub("-", "_", Treatment)]

matches[, length(unique(IndID))]
matches[, length(unique(SampleID))]
matches[, length(unique(LibraryID))]
matches[, length(unique(ExtractID))]

capt_wide <- dcast(matches, formula = IndID + SampleID +  Lat + Long + Pop + Treatment~Trip,
                   fun.aggregate = function(x) as.integer(length(unique(x))>0))

save(capt_wide, file = file.path(output.dir, "capt_wide.rda"))

capt_pop <- vector("list", length = 3)
names(capt_pop) <- capt_wide[, unique(Pop)]

for(pop in capt_wide[, unique(Pop)]) {
  capt_session <- vector("list", length = capt_wide[, length(unique(Treatment))])
  names(capt_session) <- capt_wide[, unique(Treatment)] # This ensures all session are considered but some will be NULL for some combinations
  for(ses in capt_wide[, unique(Treatment)]) {
    capt_session[[ses]] <- capt_wide[Pop == pop & Treatment == ses,]
  }
  capt_pop[[pop]] <- capt_session
}

save(capt_pop, file = file.path(output.dir, "capt_pop.rda"))

# To quickly check samples locations
matchesSR <- matches[Pop=="SnowyRiver",]
matchesLT <- matches[Pop=="LakeTyers",]
matchesANP <- matches[Pop=="Alpine",]

write.csv(matches, file.path(output.dir, "matches.csv"), row.names = F)

#------------------------------------------------------------------------------#

#### Process search tracks ####
process.search <- function(gps, nPoints=ceiling(nrow(gps)/2), lat="lat", lon="lon") {
  coords <- as.matrix(gps[, c(lon, lat)])
  coords <- spTransform(SpatialPoints(coords,
                                      proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")),
                        CRS("+proj=utm +zone=55 +south ellps=WGS84"))
  searchPoints<- spsample(Line(coords), nPoints, type="regular")
  return(searchPoints)
}

##### Snowy River #####
SR2 <- read_gpx("./GPXData/SnowyRiver_Tracklogs.gpx")

SR2pointsl <- vector("list", length = length(SR2[[2]]))
for(i in seq_along(SR2[[2]])) {
  SR2pointsl[[i]] <- process.search(SR2[[2]][[i]], nPoints=340, 
                                    lat = "Latitude", lon = "Longitude")
}
SR2points <- do.call(rbind, SR2pointsl)
plot(SR2points)

save(SR2points, file = file.path(output.dir, "SR2points.rda"))

##### ANP #####
ANP2 <- read_gpx("./GPXData/Alpine_Tracklogs.gpx")

ANP2pointsl <- vector("list", length = length(ANP2[[2]]))
for(i in seq_along(ANP2[[2]])) {
  ANP2pointsl[[i]] <- process.search(ANP2[[2]][[i]], nPoints=50, 
                                     lat = "Latitude", lon = "Longitude")
}

ANP2points <-  do.call(rbind, ANP2pointsl)
plot(ANP2points)

save(ANP2points, file = file.path(output.dir, "ANP2points.rda"))


##### LT #####
LT1 <- read_gpx("./GPXData/LakeTyers_Tracklogs_1.gpx")
LT1pointsl <- vector("list", length = length(LT1[[2]]))
for(i in seq_along(LT1[[2]])) {
  LT1pointsl[[i]] <- process.search(LT1[[2]][[i]], nPoints=50, 
                                    lat = "Latitude", lon = "Longitude")
}

LT2 <- read_gpx("./GPXData/LakeTyers_Tracklogs_2.gpx")
LT2pointsl <- vector("list", length = length(LT2[[2]]))
for(i in seq_along(LT2[[2]])) {
  LT2pointsl[[i]] <- process.search(LT2[[2]][[i]], nPoints=50, 
                                    lat = "Latitude", lon = "Longitude")
}

LT1points <- do.call(rbind, LT1pointsl)
LT2points <- do.call(rbind, LT2pointsl)
LTpoints <- rbind(LT1points, LT2points)

plot(LTpoints)

save(LTpoints, file = file.path(output.dir, "LTpoints.rda"))

