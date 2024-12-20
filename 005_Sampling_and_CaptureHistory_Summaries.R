rm(list = ls())

library(data.table)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(RColorBrewer)
library(sf)

#### Load ####
#Load matches object which contains sample metadata and ID as determined by their multi-locus genotype.
load(file.path("./Data/matches.rda"))
lapply(matches[, c("Pop", "Treatment", "Trip")], unique)
matches[, Pop := gsub(" ", "", Pop)]
matches[, Treatment := gsub("-", "_", Treatment)]
matches$Treatment <- factor(matches$Treatment, levels = c("Pre_Control", "Post_Control"))

#### Count samples and individuals ####
#Count samples again
matches[, length(unique(IndID))]
matches[, length(unique(SampleID))]
matches[, length(unique(ExtractID))]

#### Plot samples collected ####
#labels for plots
site.labs <- c("Bogong High Plains", "Lake Tyers", "Snowy River") #labels
names(site.labs) <- c("Alpine", "LakeTyers", "SnowyRiver") #values in dt

barp <- ggplot(matches, aes(as.factor(Trip), fill=Treatment)) + geom_histogram(stat = "count") +
  facet_grid(Pop~., scales = "free_x", labeller = labeller(Pop = site.labs)) + 
  ylab("Number of Samples Collected") + xlab("Trip Number") +
  scale_fill_startrek(labels = c("Pre-control", "Post-control"))
barp
# ggsave("../Figures/barplot_samplesByTreatment_trip.png", plot = barp, 
#        width = 15, height = 12, units = "cm", dpi = 300)


#### Summarise number of recaptures ####
nCaptures <- matches[,.N, by=IndID]
nCaptures[, .(Min=min(N), Median=median(N), Mean=mean(N), Max=max(N))]
nCapturesPop<-matches[,.N, by=c("IndID", "Pop")]
nCapturesPop[, .(Min=min(N), Median=median(N), Mean=mean(N), Max=max(N)), by=c("Pop")]

nCapturesBrkDown <- matches[,.N, by=c("IndID", "Pop", "Treatment", "Trip")]
nCapturesBrkDown$Treatment <- factor(nCapturesBrkDown$Treatment, levels = c("Pre_Control", "Post_Control"))
nCapturesBrkDown[, .N, by=.(Pop, Treatment)]
nCapturesBrkDown

#### Plot number of unique individuals ####
barp2 <- ggplot(nCapturesBrkDown, aes(as.factor(Trip),fill=Treatment)) + geom_histogram(stat = "count") +
  facet_grid(Pop~., scales = "free_x", labeller = labeller(Pop = site.labs)) + 
  ylab("Number of Individuals Dectected") + xlab("Trip Number") +
  scale_fill_startrek(labels = c("Pre-control", "Post-control")) + ylim(c(0,124))
barp2
# ggsave("../Figures/barplot_IndByTreatment_trip.png", plot = barp2,
#        width = 15, height = 12, units = "cm", dpi = 300)


ggsave("../Figures/barplots_combined.png", plot = ggarrange(barp, barp2, common.legend = TRUE, legend="bottom", labels="AUTO"),
       width=25, height=12, units='cm', dpi=300)



#### Number of individuals sampled across populations ####
acrossPops <- dcast(matches, formula = IndID~Pop,
                    fun.aggregate = function(x) length(unique(x))>0)
acrossPops[, acrosspops := apply(acrossPops[, 2:4, with=F], 1, sum)]
multipops <- matches[IndID %in% acrossPops[acrosspops>1, IndID], ] # Empty

multipops[, length(unique(IndID))] # 0

#### Number of individuals sampled across sessions (eg. pre and post control) ####
acrossSession <- dcast(matches, formula = IndID + Pop~Treatment,
                       fun.aggregate = function(x) length(unique(x))>0)
acrossSession[, acrosssession := apply(acrossSession[, -(1:2), with=F], 1, sum)]
multisessions <- matches[IndID %in% acrossSession[acrosssession>1, IndID], ]
setkey(multisessions, IndID)

multisessions[, length(unique(IndID)), by=Pop] #per pop

#### Number of individuals sampled across visits within a session ####
acrossVisits <- dcast(matches, formula = IndID + Pop + Treatment~Trip,
                      fun.aggregate = function(x) length(unique(x))>0)
acrossVisits[, acrossvisit := apply(acrossVisits[, -(1:3), with=F], 1, sum)]
multivisits <- matches[IndID %in% acrossVisits[acrossvisit>1, IndID], ]
setkey(multivisits, IndID)
multivisits

multivisits[, length(unique(IndID)), by=c("Pop", "Treatment")]


