rm(list = ls())

library(data.table)
library(ggplot2)
library(dartR)
library(PopGenReport)
library(poppr)
library(stringr)
library(tidyr)
source("UTM2deg.R")
source("pid.R")

data.wide=fread('./Data/MetaData_All.csv')

#### Replicate error rate ####
#Read in with replicates
deerpelletsALLgi <- read.genetable('./Data/Genotypes_ALLSwab.csv', ind = 1, pop = 2, lat = NULL, long = NULL, other.min = NULL, other.max = NULL, oneColPerAll = FALSE, sep = "/", NA.char=NA)

#Calculate distance score as count between all pairs of libraries
dist=diss.dist(deerpelletsALLgi, percent = FALSE, mat = FALSE)

dist.dt <- melt(as.matrix(dist), na.rm=TRUE) #ignore warning, its fine
colnames(dist.dt) = c("ID1", "ID2", "dist")
dist.dt = as.data.table(
  dist.dt%>%
    rowwise() %>%
    mutate(key = paste(sort(c(ID1, ID2)), collapse = ".")))
dist.dt <- copy(dist.dt[!ID1==ID2])
dist.dt <- copy(unique(setDT(dist.dt)[order(key, -ID1)], by = "key"))

dist.dt=merge(dist.dt, data.wide[,1:7], by.x="ID1", by.y="ind")
colnames(dist.dt)[5:10]=c("pop1", "Lat1", "Long1", "SampleID1", "ExtractID1", "SampleType1")
dist.dt=merge(dist.dt, data.wide[,1:7], by.x="ID2", by.y="ind")
colnames(dist.dt)[11:16]=c("pop2", "Lat2", "Long2", "SampleID2", "ExtractID2", "SampleType2")

#number of loci=22. number of alleles=44
nrow(dist.dt[SampleID1==SampleID2]) #113 replicates in final set
mean(dist.dt[SampleID1==SampleID2]$dist/44)
median(dist.dt[SampleID1==SampleID2]$dist/44)
sd(dist.dt[SampleID1==SampleID2]$dist/44)

ggplot(dist.dt[SampleID1==SampleID2], aes(x=dist/44)) + geom_histogram(bins=10) + xlab("Proportion of Allelic Mismatches") +
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=15))


#distribution of allelic mismatches across non-replicate libraries. 
hist(dist.dt[SampleID1!=SampleID2]$dist)



#### Recapture History ####
deerpelletsALL <- read.genetable('./Data/Genotypes_SwabNoReps.csv', ind = 1, pop = 2, lat = NULL, long = NULL, other.min = NULL, other.max = NULL, oneColPerAll = FALSE, sep = "/", NA.char="NA")

#Check for monomorphic loci
informloci(deerpelletsALL, cutoff = 0, MAF = 0.05, quiet = FALSE)  

# 22 loci = 44 allele. 2 point difference is 2/44 ~ 5% threshold
mlg <- mlg.filter(deerpelletsALL, threshold = 2, algorithm = "nearest_neighbour")
length(unique(mlg)) #411 unique IDs

mlg.crosspop(deerpelletsALL) #No recaptures across populations

matches <- data.table(IndID=mlg,
                      SampleLabel=row.names(deerpelletsALL@tab),
                      Pop=pop(deerpelletsALL))
matches[, length(unique(IndID)), by=Pop] #unique individuals per population

#Summarise recaptures
setkey(matches, IndID)
chk1<-matches[, .N, by=IndID]
chk1[, mean(N)] #mean number of recaptures
chk1[, sd(N)] #SD

chk2<-matches[, .N, by=c("IndID", "Pop")]
chk2[,.N, by=IndID][N>1] #No individuals captured across populations
chk2[, .(Min=min(N), Median=median(N), Mean=mean(N), Max=max(N)), by=c("Pop")]

matches <- merge(matches, data.wide[,c(1:15)], by.x="SampleLabel", by.y="ind", all.x = TRUE)

chk4<-matches[, .N, by=c("IndID", "Pop", "Treatment")]
chk4[,.N, by=IndID][N>1] #28 individuals re-captured pre and post control

chk4[, .(Min=min(N), Median=median(N), Mean=mean(N), Max=max(N)), by=c("Pop", "Treatment")]

colnames(matches)
colnames(matches)[1]="LibraryID"
colnames(matches)[15]="SamplingDate"
colnames(matches[,c(1:3, 5:9, 12:17)])
matches=matches[,c(1:3, 5:9, 12:17)]
save(matches, file = "./Data/matches.rda")




#### PID ####

#remove recaptures for PID calculation
keep=matches[!duplicated(IndID)]$LibraryID
deerpelletsUni=deerpelletsALL[i=keep]

#Calculate pid
pid=pid_calc(deerpelletsUni)
pid$pid_comb
pid$pidsibs_comb


#Permutation of pid 
# perm=pid_permute(deerpelletsUni)
# perm
# 
# print(perm$median_values, n=22)
# 
# permres=as.data.table(perm$results)
# permres[loci==21 & statistic=="PID"]$value
# 
# summary((permres[loci==21 & statistic=="PID"]$value))
# summary((permres[loci==21 & statistic=="PIDsibs"]$value))
