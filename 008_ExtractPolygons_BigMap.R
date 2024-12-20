rm(list = ls())
library(sf)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(ggpubr)
library(ggspatial)
library(grid)


#### Extract polygon areas over which density was estimated ####
# #SnowyRiver
results.dir = file.path("./ModelResults", "SnowyRiver")
load(file = file.path(results.dir, paste0("capt_hist_", "SnowyRiver",".rda")))
load(file = file.path(results.dir, paste0("mask_", "SnowyRiver", ".rda")))
sf_mask <- st_as_sf(mask, coords=c("x", "y"), crs=32755)
mask_as_grid <- st_make_grid(sf_mask, cellsize = 100, what="polygons", square=FALSE, offset = c(0,0))
sum(st_area(mask_as_grid))/1e6
st_write(mask_as_grid, "./Polygons/SnowyRiver_Polygon.shp")

# #Alpine
results.dir = file.path("./ModelResults", "Alpine")
##Polygon area size over which density was estimated.
load(file = file.path(results.dir, paste0("capt_hist_", "Alpine",".rda")))
load(file = file.path(results.dir, paste0("mask_", "Alpine", ".rda")))
sf_mask <- st_as_sf(mask, coords=c("x", "y"), crs=32755)
mask_as_grid <- st_make_grid(sf_mask, cellsize = 100, what="polygons", square=FALSE, offset = c(0,0))
sum(st_area(mask_as_grid))/1e6
st_write(mask_as_grid, "./Polygons/Alpine_Polygon.shp")
 
 
# #LakeTyers
results.dir = file.path("./ModelResults", "LakeTyers")
##Polygon area size over which density was estimated.
load(file = file.path(results.dir, paste0("capt_hist_", "LakeTyers",".rda")))
load(file = file.path(results.dir, paste0("mask_", "LakeTyers", ".rda")))
sf_mask <- st_as_sf(mask, coords=c("x", "y"), crs=32755)
mask_as_grid <- st_make_grid(sf_mask, cellsize = 100, what="polygons", square=FALSE, offset = c(0,0))
sum(st_area(mask_as_grid))/1e6
st_write(mask_as_grid, "./Polygons/LakeTyers_Polygon.shp")


#### Read in data for mapping ####
m <- read_sf("./Maps/Australian States Shapefile/States Map.shp")
parks <- read_sf("./Maps/PLM25/PLM25.shp")

SR = read_sf("./Polygons/SnowyRiver_Polygon.shp")
Al = read_sf("./Polygons/Alpine_Polygon.shp")
LT = read_sf("./Polygons/LakeTyers_Polygon.shp")
LTbox = st_as_sfc(st_bbox(LT))
Albox = st_as_sfc(st_bbox(Al))
SRbox = st_as_sfc(st_bbox(SR))

load(file.path("./Data/matches.rda"))

matches[, Pop := gsub(" ", "", Pop)]
matches[, Treatment := gsub("-", "_", Treatment)]
matches$Treatment <- factor(matches$Treatment, levels = c("Pre_Control", "Post_Control"))

mapme=matches[LibraryID!="DA2-A11", .(Lat=mean(Lat), Long=mean(Long)), by=Pop]


#### Plot ####
p1=ggplot(m) + geom_sf(fill=NA) + theme_bw() + xlim(c(140, 150)) + ylim(c(-39.5, -35)) +
  geom_sf(data=parks, aes(fill=NAME), col=NA) + 
  geom_spatial_rect(aes(xmin=146.3, xmax=148.8, ymin=-38, ymax=-36.4), fill="transparent", lty=2) +
  annotation_scale(location = "bl", width_hint = 0.4) +
  scale_fill_brewer(palette = "Dark2", labels=c("Alpine NP", "Lake Tyers SP", "Snowy River NP"), name="") +
  theme(legend.position=c(0.85, 0.2),
        axis.text = element_text(size=12),
        axis.title = element_blank(),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))

p2=ggplot() + geom_sf(data=m, fill=NA) + geom_sf(data=parks, aes(fill=NAME), col=NA, alpha=0.8) + 
  geom_sf(data=LTbox, aes(fill=NA), linewidth=1, col="black") +
  geom_sf(data=Albox, aes(fill=NA), linewidth=1, col="black") + 
  geom_sf(data=SRbox, aes(fill=NA), linewidth=1, col="black") +
  #scale_y_continuous(breaks=C(-37.8, -37.4, -37, -36.6)) +
  annotation_scale(location = "bl", width_hint = 0.4) +
  annotate("text", x=147.6, y= -37, label="Bogong High Plains", fontface =2, size=4.5) +
  annotate("text", x=148.08, y= -37.2, label="Snowy River", fontface =2, size=4.5) +
  annotate("text", x=147.83, y= -37.8, label="Lake Tyers", fontface =2, size=4.5) +
  scale_fill_brewer(palette = "Dark2", labels=c("Alpine NP", "Lake Tyers SP", "Snowy River NP"), name="") +
  theme_bw() +
  xlim(c(146.4, 148.7)) + ylim(c(-37.9, -36.5)) +
  theme(legend.position="none",
        axis.text = element_text(size=12),
        axis.title = element_blank())


fig=ggarrange(p1, p2, ncol=1, widths = c(10, 5))


mapfig=annotate_figure(fig, left = textGrob("Longitude", rot = 90, vjust = 4, gp = gpar(cex = 1.3)),
                       bottom = textGrob("Latitude", gp = gpar(cex = 1.3)))



ggsave("../Figures/mainmap.png", mapfig, width=26.5, height=26.5, units='cm', dpi=300, bg="white")



#plot of Australia with Vic outlined. Use Inkscape to add this to mainmap.
ggplot(m) + 
  geom_sf(aes(fill=NULL)) +
  geom_rect(aes(xmin=140, ymin=-39.5, 
                xmax=150, ymax=-35),
            lwd=1, col="black", fill=NA) +
  ylab("Latitude") + xlab("Longitude") +
  theme_bw() +
  xlim(c(114, 154)) + ylim(c(-44, -11)) +
  theme(#legend.position = "none",
    axis.text=element_blank(),
    axis.title = element_blank(),
    legend.spacing.y = unit(0.5, 'cm'),
    plot.background = element_rect(fill=NA, colour=NA))

