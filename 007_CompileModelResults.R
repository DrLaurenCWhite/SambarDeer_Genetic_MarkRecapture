rm(list = ls())
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

sites <- c("Alpine", "LakeTyers", "SnowyRiver")

fl <- list.files(file.path("./ModelResults", sites), pattern = "^Par_est_CI", full.names = T)

#Best models:
#Alpine: Alpine_static_HHR_
#SnowyRiver: Snowy_River_BVN_HHR_t
#LakeTyers: LakeTyers_static_HHR_

fl <- fl[grep(pattern = "SnowyRiver_BVN|Alpine_static_HHR_\\.csv|LakeTyers_static_HHR_\\.csv", x = fl)]

res <- lapply(fl, read.csv)
  
names(res) <- sites
res <- rbindlist(res, use.names = T, fill = T, idcol="Population")
  
est_vec <- res[, paste(estimate, paste0("(", lcl, " - ", ucl, ")"))]
  
res_range <- cbind(res[, 1:3], est_vec)
res_range[, session := factor(session, levels = c("Pre_Control", "Post_Control"))]
res_cast <- dcast(res_range, formula =  Population + session ~ Param)

setnames(res_cast, "session", "Treatment")

write.csv(res_cast, file = file.path("./Tables/BestModelResults.csv"), row.names = F)



#Plot best model results : density and lambda
SR = fread("./ModelResults/SnowyRiver/Par_est_CI_SnowyRiver_BVN_HHR_t.csv")
LT = fread("./ModelResults/LakeTyers/Par_est_CI_LakeTyers_static_HHR_.csv")
AL = fread("./ModelResults/Alpine/Par_est_CI_Alpine_static_HHR_.csv")

SR = SR[Param %in% c("phi", "lambda", "D")]
LT = LT[Param %in% c("phi", "lambda", "D")]
AL = AL[Param %in% c("phi", "lambda", "D")]

SR$Site="Snowy River"
LT$Site = "Lake Tyers"
AL$Site = "Alpine"

dt=rbind(SR, LT, AL)

dt = dt[, c(1:3, 5:7)]

dt.w = dcast(dt, Site + session ~ Param, value.var=c("estimate", "ucl", "lcl"))

dt.w$session <- factor(dt.w$session, levels = c("Pre_Control", "Post_Control"))
dt.w[Site=="Alpine", Site:="Bogong High Plains"]

p1=ggplot(dt.w, aes(x=session, y=estimate_D, group=1, col=Site)) + geom_point(size=4) + 
  geom_linerange(aes(ymin=lcl_D, ymax=ucl_D), linewidth=1) +
  scale_color_brewer(palette="Dark2") + stat_summary(fun.y=sum, geom="line", lty=5) +
  facet_grid(.~Site) + ylab(bquote('Deer Density '(km^2))) + xlab("") +
  scale_x_discrete(labels=c("Pre-Control", "Post-Control")) +
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        axis.text = element_text(size=15),
        axis.title = element_text(size=15),
        strip.background = element_blank(),
        strip.text.x = element_blank())

p2=ggplot(dt.w[session=="Pre_Control"], aes( x = Site, y=estimate_lambda, col=Site))+ geom_point(size=4) + 
  geom_linerange(aes(ymin=lcl_lambda, ymax=ucl_lambda), linewidth=1) + 
  geom_hline(yintercept = 1, lty=5, col="red", linewidth=1, alpha=0.5) +
  scale_color_brewer(palette="Dark2") + ylab("Population Rate of Change") +
  xlab("") + 
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        axis.text = element_text(size=15),
        axis.title = element_text(size=15))


ggarrange(p1, p2, ncol=1)

ggsave("./Figures/MainModelRes.png", plot = ggarrange(p1, p2, ncol=1),
       width=25, height=15, units='cm', dpi=300)


