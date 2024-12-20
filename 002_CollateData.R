rm(list = ls())

library(data.table)
library(ggplot2)
library(dartR)
library(PopGenReport)
library(poppr)
library(stringr)
source("UTM2deg.R")


#### Import haplo-type data ####
dat=fread("./Data/Haplotype_table.csv", header=TRUE)
dat$Sample_ID <- sub("_[^_]+$", "", dat$Sample_ID)
dat.long=melt(dat,id.vars = "Sample_ID",value.name="Genotype",variable.name = "Locus")
dat$Plate=str_extract(dat$Sample_ID, "[^-]+")

#### Import library meta data ####
meta = fread('./Data/Library_MetaData.csv', na.strings = c("NA", ""))

#### Import sample meta data ####
samp = fread('./Data/Sample_MetaData.csv', na.strings = c("NA", ""))

#### Convert UTM to lat long ####
samp[is.na(Lat)][,c("Lat", "Long")] = samp[is.na(Lat),UTM2deg(as.numeric(Easting), as.numeric(Northing), zone = 55)]

#### count samples ####
samp[CollectionBuffer=="Longmires", .N] #Samples collected
samp[CollectionBuffer=="Longmires" & !is.na(ExtractID), .N] #Samples extracted
length(unique(meta[SampleType=="Swab"]$ExtractID)) #number of samples sequenced

meta[SampleType=="Swab" & ReplicateLibraries=="Y", .N] #number of swab library replicates
meta[SampleType=="Swab" & is.na(ReplicateLibraries), .N] #number of swab libraries without replicates
meta[SampleType=="Swab", .N] #number of libraries total
meta[, .N, by=SampleType] #libraries by sample type 


#merge data types and count again to check
dat.long=merge(meta, dat.long, by.x="LibraryID", by.y="Sample_ID", all = TRUE)
length(unique(dat.long$LibraryID))
length(unique(dat.long[SampleType=="Swab"]$SampleID))
length(unique(dat.long[SampleType=="Swab" & ReplicateLibraries=="Y"]$SampleID))
length(unique(dat.long[SampleType=="Swab" & is.na(ReplicateLibraries)]$SampleID))

dat.long=merge(dat.long, samp,  by=c("SampleID", "ExtractID"), all.x=TRUE)
length(unique(dat.long$LibraryID))
length(unique(dat.long[SampleType=="Swab"]$SampleID))
length(unique(dat.long[SampleType=="Swab" & ReplicateLibraries=="Y"]$SampleID))
length(unique(dat.long[SampleType=="Swab" & is.na(ReplicateLibraries)]$SampleID))
length(unique(dat.long[SampleType=="Concentration-blank" | SampleType=="Extraction-blank" |
                         SampleType=="Index-blank" | SampleType=="Amplicon-blank"]$SampleID))


#### Calculate and filter on missing data/sample type ####
#Calculate the % missing per Loci
#34 loci genotyped
dat.long[is.na(Genotype), .N/length(unique(dat.long$LibraryID)), by=.(Locus)]
ggplot(dat.long[is.na(Genotype), .N/length(unique(dat.long$LibraryID)), by=.(Locus)], aes(x=V1)) +
  geom_histogram()

#Exclude loci with >25% missing data
excludeLocus=dat.long[is.na(Genotype), .N/length(unique(dat.long$LibraryID)), by=.(Locus)][V1>0.25]$Locus
droplevels(excludeLocus)
#Excluded 12 loci with >25% missing data (these same 12 needed to be excluded even if individuals were filtered first)
dat.long=dat.long[!Locus %in% excludeLocus]


#Calculate % missing data per library
#22 loci genotyped
dat.long[is.na(Genotype), .N, by=.(LibraryID, SampleID, SampleType)][, mean(N/22), by=SampleType]


#Exclude everything except for swab samples (collected in longmires buffer, not ethanol)
dat.swab=dat.long[SampleType=="Swab"]
length(unique(dat.swab$LibraryID)) #904 swab libraries genotyped total
length(unique(dat.swab[ReplicateLibraries=="Y"]$LibraryID)) #143 of which are replicates
length(unique(dat.swab[is.na(ReplicateLibraries)]$LibraryID)) #761 unique samples genotyped


#Exclude libraries with >25% missing data
ggplot(dat.swab[is.na(Genotype), .N, by=.(LibraryID, SampleID, SampleType)], aes(x=N/22)) + 
  geom_histogram()
excludeLib=dat.swab[is.na(Genotype), .N/22, by=.(LibraryID, SampleID, SampleType)][V1>0.25]$LibraryID
length(excludeLib)
#Excluded 140 Libraries with >25% missing data

dat.swab=dat.swab[!LibraryID %in% excludeLib]
length(unique(dat.swab$LibraryID)) #764 swab libraries kept after filtering (including some replicates)


#average and sd percent missing
dat.swab[is.na(Genotype), .N, by=.(LibraryID, SampleID, SampleType)][, mean(N/22), by=SampleType]
dat.swab[is.na(Genotype), .N, by=.(LibraryID, SampleID, SampleType)][, sd(N/22), by=SampleType]




#### Format and export ####
dat.swab$QualScore = nchar(dat.swab$Quality)
data.wide <- dcast(dat.swab, 
                   SampleID + ExtractID + LibraryID + SampleType + Nanodrop_Conc + ExtractqPCRConc +
                     ReplicateLibraries + Location + Treatment + Trip +PelletID + Lat +Long + 
                     Date + QualScore  + Plate~ Locus, value.var="Genotype")

data.wide=data.wide[,c(3,8, 12:13, 1:2, 4:6, 9:11, 14:39)]
colnames(data.wide)[c(1,2)]=c("ind", "pop")

data.wide$check=duplicated(data.wide$SampleID)

data.wide[check==FALSE, .N]
data.wide[check==TRUE, .N]
data.wide[,17:38]
data.wide[check==FALSE, .N, by=.(pop, Treatment)]
data.wide[check==FALSE, .N, by=.(pop)]

ggplot(data.wide[check==FALSE, .N, by=.(pop, Treatment, Trip)], aes(x=Trip, y=N)) +
  geom_bar(stat="identity") + facet_grid(.~pop)


#export as csv#
#MetaData
write.csv(data.wide[,1:15], './Data/MetaData_All.csv', row.names = FALSE)

#All swab libraries and replicates
write.csv(data.wide[,c(1:2,17:38)], './Data/Genotypes_ALLSwab.csv', row.names = FALSE)

#Swabs without replicates
write.csv(data.wide[check==FALSE][,c(1:2,17:38)], 
          './Data/Genotypes_SwabNoReps.csv', row.names = FALSE)


