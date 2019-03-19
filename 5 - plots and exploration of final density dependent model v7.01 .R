# exploration of final density dependent habitat selection model
##################################################################################################################################################################
# set wd ----
wrk.dir = "T:\\HMP\\HMP work\\Analysis"
model.outputs.path = "T:\\HMP\\HMP work\\Analysis\\objects.corrections\\"
setwd(wrk.dir)

##################################################################################################################################################################
# load libraries ----
library(viridis)
library(mgcv)
library(gtools)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(data.table)
library(MuMIn)
library(TTR)
library(smooth)
library(Mcomp)
library(GGally)
library("cowplot")

##################################################################################################################################################################
# load .csvs ----
par(mfrow = c(1,1))
def.par = par()
unsex.trap = read.csv("unsex.trap.v7.01.csv") 
#unsex.trap$sesid = as.factor(unsex.trap$sesid)
##################################################################################################################################################################
# run functions code -- needs to be done ----
source("functions.code.R")
##################################################################################################################################################################
#load models ------
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.2.0")
##################################################################################################################################################################
#
#summary(bam.st.habseason.hab.density.2.0) 
traps.per.night = table(unsex.trap$date)

traps.per.typical.night = table(unsex.trap$date[unsex.trap$sub.sesh.night.num<5])
##################################################################################################################################################################
# plot of catpture rates by trap night and sub session ----

intrest.covar = "sub.sesh.night.num"
expl.lab =   "Trap night within trapping session"
ylims = c(0 , 0.012)
interested.length = 700

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                     interested.covar = intrest.covar, interested.length = interested.length, 
                                     by.fac = F, interesting.factor.name = NULL)
my.data = my.data1[rep(1:dim(my.data1)[1], 4),] # rep that data four times
my.data$sub.ses.num = rep(1:4, each = dim(my.data1)[1])

pred.bam.st.habseason.hab.density.2.0 = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")
my.data$fit = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit))
my.data$upper.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit + 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))
my.data$lower.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit - 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))

my.data$expl = my.data[,which(names(my.data)==intrest.covar )]
my.data$expl = seq(1,max(my.data$sub.ses.num)*max(my.data$sub.sesh.night.num), length= dim(my.data)[1] )

my.data$sub.ses.num <- as.character(my.data$sub.ses.num)

ggplot(my.data, aes(expl))+
    geom_ribbon(aes(ymin = lower.95, ymax = upper.95, fill = sub.ses.num) #, alpha = 0.6
                ) +
    geom_line(aes(y = fit, x = expl, group = sub.ses.num ), size = 1, col = "black")+
    scale_fill_viridis("Sub-session\nnumber", 1:4, discrete = T)+
    labs(x = expl.lab, y = "Trapping rate") +
    #scale_y_continuous(limits=ylims) +
    scale_x_continuous(limits=c(0,max(my.data$expl))) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )

##################################################################################################################################################################
# shelter.plot -----

intrest.covar = "salat.sum.means"
expl.lab =   'Coastal shelter index'
interested.length = 7000

breaks.by = 100
break.divider = 1

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                     interested.covar = intrest.covar, interested.length = interested.length, 
                                     by.fac = F, interesting.factor.name = NULL)
my.data = my.data1
my.data[,which(names(my.data)==intrest.covar)] = seq(0, max(unsex.trap[,which(names(unsex.trap)==intrest.covar)]), length = dim(my.data)[1]) 


pred.bam.st.habseason.hab.density.2.0 = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")
my.data$fit = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit))
my.data$upper.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit + 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))
my.data$lower.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit - 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))

my.data$expl = my.data[,which(names(my.data)==intrest.covar)]
rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])
ylims = c(0 , round(max(my.data$upper.95), digits = 3)+0.001)
x.max = max(my.data$expl)


shelter.plot = ggplot()+
    geom_ribbon(aes(x = my.data$expl, ymin = my.data$lower.95, ymax = my.data$upper.95), fill = "grey80"#, alpha = 0.4
    ) +
    geom_line(aes(y = my.data$fit, x = my.data$expl), size = 1, col = "black")+
    labs(x = expl.lab, y = "Trapping rate") +
    geom_rug( aes(x = rugs),sides = "b",# position = "jitter",
              size = 0.01) +
    scale_y_continuous(limits=ylims) +
    scale_x_continuous(breaks=seq(0,x.max, by = breaks.by ), labels= seq(0,max(my.data$expl), by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )

main.plot = ggplot()+
    geom_ribbon(aes(x = my.data$expl, ymin = my.data$lower.95, ymax = my.data$upper.95), fill = "grey80"#, alpha = 0.4
    ) +
    geom_line(aes(y = my.data$fit, x = my.data$expl), size = 1, col = "black")+
    labs(x = "", y = "Trapping rate") +
    geom_rug( aes(x = rugs),sides = "b",# position = "jitter",
              size = 0.01) +
    scale_y_continuous(limits=ylims) +
    scale_x_continuous(breaks=seq(0,x.max, by = breaks.by ), labels= seq(0,max(my.data$expl), by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )
main.plot

bottom.plot = ggplot(unsex.trap[unsex.trap[,which(names(unsex.trap)==intrest.covar)] !=0,], aes(combo.heather_all.mk2))+
    geom_density(fill = "grey90") +
    labs(x = expl.lab, y = "Density") +
    theme_bw() +
    scale_x_continuous(breaks=seq(0,max(my.data$expl), by = breaks.by ), labels= seq(0,max(my.data$expl), by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
    scale_y_continuous (labels = NULL, breaks = NULL) +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )
bottom.plot

grid.arrange(main.plot, bottom.plot, ncol=1, nrow=2, heights=c(4, 1))

##################################################################################################################################################################
# heather.plot -----

intrest.covar = "combo.heather_all.mk2"
expl.lab =   'Heather area within 1.5km ('~km^2*')'
interested.length = 7000

breaks.by = 1000000
break.divider = 1000000

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                     interested.covar = intrest.covar, interested.length = interested.length, 
                                     by.fac = F, interesting.factor.name = NULL)
my.data = my.data1
my.data[,which(names(my.data)==intrest.covar)] = seq(0, max(unsex.trap[,which(names(unsex.trap)==intrest.covar)]), length = dim(my.data)[1]) 


pred.bam.st.habseason.hab.density.2.0 = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")
my.data$fit = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit))
my.data$upper.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit + 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))
my.data$lower.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit - 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))

my.data$expl = my.data[,which(names(my.data)==intrest.covar)]
rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])
ylims = c(0 , round(max(my.data$upper.95), digits = 3)+0.001)
x.max = max(my.data$expl)


heather.plot = ggplot()+
    geom_ribbon(aes(x = my.data$expl, ymin = my.data$lower.95, ymax = my.data$upper.95), fill = "grey80"#, alpha = 0.4
    ) +
    geom_line(aes(y = my.data$fit, x = my.data$expl), size = 1, col = "black")+
    labs(x = expl.lab, y = "Trapping rate") +
    geom_rug( aes(x = rugs),sides = "b",# position = "jitter",
              size = 0.01) +
    scale_y_continuous(limits=ylims) +
    scale_x_continuous(breaks=seq(0,x.max, by = breaks.by ), labels= seq(0,max(my.data$expl), by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )

main.plot = ggplot()+
    geom_ribbon(aes(x = my.data$expl, ymin = my.data$lower.95, ymax = my.data$upper.95), fill = "grey80"#, alpha = 0.4
    ) +
    geom_line(aes(y = my.data$fit, x = my.data$expl), size = 1, col = "black")+
    labs(x = "", y = "Trapping rate") +
    geom_rug( aes(x = rugs),sides = "b",# position = "jitter",
              size = 0.01) +
    scale_y_continuous(limits=ylims) +
    scale_x_continuous(breaks=seq(0,x.max, by = breaks.by ), labels= seq(0,max(my.data$expl), by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )
main.plot

bottom.plot = ggplot(unsex.trap[unsex.trap[,which(names(unsex.trap)==intrest.covar)] !=0,], aes(combo.heather_all.mk2))+
    geom_density(fill = "grey90") +
    labs(x = expl.lab, y = "Density") +
    theme_bw() +
    scale_x_continuous(breaks=seq(0,max(my.data$expl), by = breaks.by ), labels= seq(0,max(my.data$expl), by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
    scale_y_continuous (labels = NULL, breaks = NULL) +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )
bottom.plot

grid.arrange(main.plot, bottom.plot, ncol=1, nrow=2, heights=c(4, 1))



##################################################################################################################################################################
# builtup.plot -----

intrest.covar = "builtup.mk2"
expl.lab =   'Infrastructure area within 1.5km ('~km^2*')'
interested.length = 7000

breaks.by = 1000000
break.divider = 1000000

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                     interested.covar = intrest.covar, interested.length = interested.length, 
                                     by.fac = F, interesting.factor.name = NULL)
my.data = my.data1
my.data[,which(names(my.data)==intrest.covar)] = seq(0, max(unsex.trap[,which(names(unsex.trap)==intrest.covar)]), length = dim(my.data)[1]) 


pred.bam.st.habseason.hab.density.2.0 = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")
my.data$fit = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit))
my.data$upper.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit + 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))
my.data$lower.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit - 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))

my.data$expl = my.data[,which(names(my.data)==intrest.covar)]
rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])
ylims = c(0 , round(max(my.data$upper.95), digits = 3)+0.001)
x.max = max(my.data$expl)


builtup.plot = ggplot()+
    geom_ribbon(aes(x = my.data$expl, ymin = my.data$lower.95, ymax = my.data$upper.95), fill = "grey80"#, alpha = 0.4
    ) +
    geom_line(aes(y = my.data$fit, x = my.data$expl), size = 1, col = "black")+
    labs(x = expl.lab, y = "Trapping rate") +
    geom_rug( aes(x = rugs),sides = "b",# position = "jitter",
              size = 0.01) +
    scale_y_continuous(limits=ylims) +
    scale_x_continuous(breaks=seq(0,x.max, by = breaks.by ), labels= seq(0,max(my.data$expl), by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )

main.plot = ggplot()+
    geom_ribbon(aes(x = my.data$expl, ymin = my.data$lower.95, ymax = my.data$upper.95), fill = "grey80"#, alpha = 0.4
    ) +
    geom_line(aes(y = my.data$fit, x = my.data$expl), size = 1, col = "black")+
    labs(x = "", y = "Trapping rate") +
    geom_rug( aes(x = rugs),sides = "b",# position = "jitter",
              size = 0.01) +
    scale_y_continuous(limits=ylims) +
    scale_x_continuous(breaks=seq(0,x.max, by = breaks.by ), labels= seq(0,max(my.data$expl), by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )
main.plot

bottom.plot = ggplot(unsex.trap[unsex.trap[,which(names(unsex.trap)==intrest.covar)] !=0,], aes(combo.heather_all.mk2))+
    geom_density(fill = "grey90") +
    labs(x = expl.lab, y = "Density") +
    theme_bw() +
    scale_x_continuous(breaks=seq(0,max(my.data$expl), by = breaks.by ), labels= seq(0,max(my.data$expl), by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
    scale_y_continuous (labels = NULL, breaks = NULL) +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )
bottom.plot

grid.arrange(main.plot, bottom.plot, ncol=1, nrow=2, heights=c(4, 1))


##################################################################################################################################################################
# beach.plot -----

intrest.covar = "combo.beach_saltmarsh.mk2"
expl.lab =   'Beach area within 1.5km ('~km^2*')'
interested.length = 7000

breaks.by = 1000000
break.divider = 1000000

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                     interested.covar = intrest.covar, interested.length = interested.length, 
                                     by.fac = F, interesting.factor.name = NULL)
my.data = my.data1
my.data[,which(names(my.data)==intrest.covar)] = seq(0, max(unsex.trap[,which(names(unsex.trap)==intrest.covar)]), length = dim(my.data)[1]) 


pred.bam.st.habseason.hab.density.2.0 = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")
my.data$fit = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit))
my.data$upper.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit + 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))
my.data$lower.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit - 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))

my.data$expl = my.data[,which(names(my.data)==intrest.covar)]
rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])
ylims = c(0 , round(max(my.data$upper.95), digits = 3)+0.001)
x.max = max(my.data$expl)


beach.plot = ggplot()+
    geom_ribbon(aes(x = my.data$expl, ymin = my.data$lower.95, ymax = my.data$upper.95), fill = "grey80"#, alpha = 0.4
    ) +
    geom_line(aes(y = my.data$fit, x = my.data$expl), size = 1, col = "black")+
    labs(x = expl.lab, y = "Trapping rate") +
    geom_rug( aes(x = rugs),sides = "b",# position = "jitter",
              size = 0.01) +
    scale_y_continuous(limits=ylims) +
    scale_x_continuous(breaks=seq(0,x.max, by = breaks.by ), labels= seq(0,max(my.data$expl), by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )

main.plot = ggplot()+
    geom_ribbon(aes(x = my.data$expl, ymin = my.data$lower.95, ymax = my.data$upper.95), fill = "grey80"#, alpha = 0.4
    ) +
    geom_line(aes(y = my.data$fit, x = my.data$expl), size = 1, col = "black")+
    labs(x = "", y = "Trapping rate") +
    geom_rug( aes(x = rugs),sides = "b",# position = "jitter",
              size = 0.01) +
    scale_y_continuous(limits=ylims) +
    scale_x_continuous(breaks=seq(0,x.max, by = breaks.by ), labels= seq(0,max(my.data$expl), by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )
main.plot

bottom.plot = ggplot(unsex.trap[unsex.trap[,which(names(unsex.trap)==intrest.covar)] !=0,], aes(combo.heather_all.mk2))+
    geom_density(fill = "grey90") +
    labs(x = expl.lab, y = "Density") +
    theme_bw() +
    scale_x_continuous(breaks=seq(0,max(my.data$expl), by = breaks.by ), labels= seq(0,max(my.data$expl), by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
    scale_y_continuous (labels = NULL, breaks = NULL) +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )
bottom.plot

grid.arrange(main.plot, bottom.plot, ncol=1, nrow=2, heights=c(4, 1))

##################################################################################################################################################################
# bog.plot -----

intrest.covar = "bog.mk2"
expl.lab =   'Bog area within 1.5km ('~km^2*')'
interested.length = 7000

breaks.by = 1000000
break.divider = 1000000

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                     interested.covar = intrest.covar, interested.length = interested.length, 
                                     by.fac = F, interesting.factor.name = NULL)
my.data = my.data1
my.data[,which(names(my.data)==intrest.covar)] = seq(0, max(unsex.trap[,which(names(unsex.trap)==intrest.covar)]), length = dim(my.data)[1]) 


pred.bam.st.habseason.hab.density.2.0 = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")
my.data$fit = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit))
my.data$upper.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit + 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))
my.data$lower.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit - 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))

my.data$expl = my.data[,which(names(my.data)==intrest.covar)]
rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])
ylims = c(0 , round(max(my.data$upper.95), digits = 3)+0.001)
x.max = max(my.data$expl)


bog.plot = ggplot()+
    geom_ribbon(aes(x = my.data$expl, ymin = my.data$lower.95, ymax = my.data$upper.95), fill = "grey80"#, alpha = 0.4
    ) +
    geom_line(aes(y = my.data$fit, x = my.data$expl), size = 1, col = "black")+
    labs(x = expl.lab, y = "Trapping rate") +
    geom_rug( aes(x = rugs),sides = "b",# position = "jitter",
              size = 0.01) +
    scale_y_continuous(limits=ylims) +
    scale_x_continuous(breaks=seq(0,x.max, by = breaks.by ), labels= seq(0,max(my.data$expl), by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )

main.plot = ggplot()+
    geom_ribbon(aes(x = my.data$expl, ymin = my.data$lower.95, ymax = my.data$upper.95), fill = "grey80"#, alpha = 0.4
    ) +
    geom_line(aes(y = my.data$fit, x = my.data$expl), size = 1, col = "black")+
    labs(x = "", y = "Trapping rate") +
    geom_rug( aes(x = rugs),sides = "b",# position = "jitter",
              size = 0.01) +
    scale_y_continuous(limits=ylims) +
    scale_x_continuous(breaks=seq(0,x.max, by = breaks.by ), labels= seq(0,max(my.data$expl), by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )
main.plot

bottom.plot = ggplot(unsex.trap[unsex.trap[,which(names(unsex.trap)==intrest.covar)] !=0,], aes(combo.heather_all.mk2))+
    geom_density(fill = "grey90") +
    labs(x = expl.lab, y = "Density") +
    theme_bw() +
    scale_x_continuous(breaks=seq(0,max(my.data$expl), by = breaks.by ), labels= seq(0,max(my.data$expl), by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
    scale_y_continuous (labels = NULL, breaks = NULL) +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )
bottom.plot

grid.arrange(main.plot, bottom.plot, ncol=1, nrow=2, heights=c(4, 1))

##################################################################################################################################################################
# plot habitat effects without interactions ---- 
grid.arrange(shelter.plot, heather.plot, builtup.plot, beach.plot,bog.plot , ncol=2, nrow=3)






##################################################################################################################################################################
#
# MAIN DENSITY PREDICTION density -----

intrest.covar = "log.dens"
expl.lab =  'Regional abundance'
interested.length = 7000
divider = 1
dens.to.axis = c(100,500,1000, 2000)


my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                     interested.covar = intrest.covar, interested.length = interested.length, 
                                     by.fac = F, interesting.factor.name = NULL)
my.data = my.data1
#my.data[,which(names(my.data)==intrest.covar)]  = log(1+ seq(0, exp(max(unsex.trap[,which(names(unsex.trap)==intrest.covar)])), length = dim(my.data)[1] ))


pred.bam.st.habseason.hab.density.2.0 = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")
my.data$fit = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit))
my.data$upper.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit + 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))
my.data$lower.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit - 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))

#my.data$expl = exp(my.data[,which(names(my.data)==intrest.covar)]  )-1
my.data$expl = exp(my.data[,which(names(my.data)==intrest.covar)]  )-1 

rugs = exp(unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)]))-1 
ylims = c(0 , round(max(my.data$upper.95), digits = 3)+0.001)

ggplot()+
    geom_ribbon(aes(x = my.data$expl, ymin = my.data$lower.95, ymax = my.data$upper.95), fill = "grey80"#, alpha = 0.4
    ) +
    geom_line(aes(y = my.data$fit, x = my.data$expl), size = 1, col = "black")+
    labs(x = expl.lab, y = "Trapping rate") +
    geom_rug( aes(x = rugs),sides = "b",# position = "jitter",
              size = 0.01) +
    scale_y_continuous(limits=ylims) +
    scale_x_continuous(breaks=dens.to.axis, labels= dens.to.axis ) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )


##################################################################################################################################################################
#
#
# ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
# ti(julian, log.builtup.mk2, bs = c("cc", "cr"), d = c(1,1)) +
# ti(julian, coast.dist, bs = c("cc", "cr"), d = c(1,1)) +
# ti(julian, log.pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
# ti(julian, combo.bog_all.mk2, bs = c("cc", "cr"), d = c(1,1)) +
the.interested.length = 300
##########################################################################################
# main effect of seaonality ----
intrest.covar = "julian"
intrest.term = "ti(julian)"
expl.lab =   "Day of year"
interested.length = the.interested.length

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                     interested.covar = intrest.covar, interested.length = interested.length, 
                                     by.fac = F, interesting.factor.name = NULL)
my.data = my.data1
#my.data[,which(names(my.data)==intrest.covar)]  = log(1+ seq(0, exp(max(unsex.trap[,which(names(unsex.trap)==intrest.covar)])), length = dim(my.data)[1] ))
pred.bam.st.habseason.hab.density.2.0.dum = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "terms")

dum.se = as.numeric(pred.bam.st.habseason.hab.density.2.0.dum$se.fit[,which(attributes(pred.bam.st.habseason.hab.density.2.0.dum$fit)$dimnames[[2]] == intrest.term)])
pred.bam.st.habseason.hab.density.2.0 = pred.bam.st.habseason.hab.density.2.0.dum

my.data$fit = as.numeric( pred.bam.st.habseason.hab.density.2.0.dum$fit[,which(attributes(pred.bam.st.habseason.hab.density.2.0.dum$fit)$dimnames[[2]] == intrest.term)])
my.data$upper.95 = my.data$fit - (dum.se*1.96)
my.data$lower.95 = my.data$fit + (dum.se*1.96)
#my.data$expl = exp(my.data[,which(names(my.data)==intrest.covar)]  )-1
my.data$expl = my.data[,which(names(my.data)==intrest.covar)]  

rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])
ylims = c(0 , round(max(my.data$upper.95), digits = 3)+0.001)

    ggplot()+
    geom_ribbon(aes(x = my.data$expl, ymin = my.data$lower.95, ymax = my.data$upper.95), fill = "grey80"#, alpha = 0.4
                )+
    geom_line(aes(y = my.data$fit, x = my.data$expl), size = 1, col = "black")+
    labs(x = expl.lab , y = "Effect size") +
    scale_y_continuous() +
    scale_x_continuous(breaks= c(1/12, 3/12, 5/12, 7/12, 9/12,11/12)*365, labels= c("Feb", "Apr", "June", "Aug","Oct", "Dec")) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )

##########################################################################################
# main effect of DENSITY -- LOGGED -- ----
intrest.covar = "log.dens"
intrest.term = "ti(log.dens)"
expl.lab =  "Regional abundance (N)"
interested.length = the.interested.length
divider = 1
dens.to.axis = c(100,500,1000, 2000)

#min.expl = min(unsex.trap[,names(unsex.trap)==intrest.covar ])
#max.expl = max(unsex.trap[,names(unsex.trap)==intrest.covar ])
#max.expl = 250 # 2 traps over 250m, 5 sessions
min.expl = exp(min(unsex.trap[,names(unsex.trap)==intrest.covar ]))-1
max.expl = exp(max(unsex.trap[,names(unsex.trap)==intrest.covar ]))-1

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                     interested.covar = intrest.covar, interested.length = interested.length, 
                                     by.fac = F, interesting.factor.name = NULL)
my.data1$exp.expl = seq(min.expl, max.expl, length = dim(my.data1)[1])
my.data1[, which(names(my.data1)== intrest.covar)] = log(my.data1$exp.expl+1)

my.data = my.data1
#my.data = my.data1[rep(1:(dim(my.data1)[1]), times = interested.length ),]
#my.data[, which(names(my.data)== itneracting)] = rep (seq(min.inter, max.inter, length = interested.length), each = interested.length)

pred.bam.st.habseason.hab.density.2.0.dum = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "terms")

dum.se = as.numeric(pred.bam.st.habseason.hab.density.2.0.dum$se.fit[,which(attributes(pred.bam.st.habseason.hab.density.2.0.dum$fit)$dimnames[[2]] == intrest.term)])
pred.bam.st.habseason.hab.density.2.0 = pred.bam.st.habseason.hab.density.2.0.dum

my.data$fit = as.numeric( pred.bam.st.habseason.hab.density.2.0.dum$fit[,which(attributes(pred.bam.st.habseason.hab.density.2.0.dum$fit)$dimnames[[2]] == intrest.term)])
my.data$upper.95 = my.data$fit - (dum.se*1.96)
my.data$lower.95 = my.data$fit + (dum.se*1.96)
#my.data$expl = exp(my.data[,which(names(my.data)==intrest.covar)]  )-1
my.data$expl = exp(my.data[,which(names(my.data)==intrest.covar)]  )-1
#my.data$expl = my.data[,which(names(my.data)==intrest.covar)]  
rugs = exp(unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)]))-1
#rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])

ylims = c(0 , round(max(my.data$upper.95), digits = 3)+0.001)


main.plot = 
    ggplot()+
    geom_ribbon(aes(x = my.data$expl, ymin = my.data$lower.95, ymax = my.data$upper.95), fill = "grey80"#, alpha = 0.4
    )+
    geom_line(aes(y = my.data$fit, x = my.data$expl), size = 1, col = "black")+
    labs(x = expl.lab, y = "Effect size") +
    geom_rug( aes(x = rugs),sides = "b",# position = "jitter",
              size = 0.01) +
    scale_y_continuous() +
    scale_x_continuous(breaks=  dens.to.axis, labels= dens.to.axis) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )
main.plot

bottom.plot = ggplot( unsex.trap, aes(block.dens))+
    geom_density(fill = "grey90") +
    labs(x = expl.lab, y = "Density") +
    theme_bw() +
    scale_x_continuous(breaks=  dens.to.axis, labels= dens.to.axis) +
    scale_y_continuous (labels = NULL, breaks = NULL) +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )
grid.arrange(main.plot, bottom.plot, ncol=1, nrow=2, heights=c(4, 1))



##########################################################################################
# main effect of altitude ----
intrest.covar = "ALT3"
intrest.term = "ti(ALT3)"
expl.lab =   "Altitude (m)"
interested.length = the.interested.length

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                     interested.covar = intrest.covar, interested.length = interested.length, 
                                     by.fac = F, interesting.factor.name = NULL)
my.data = my.data1
#my.data[,which(names(my.data)==intrest.covar)]  = log(1+ seq(0, exp(max(unsex.trap[,which(names(unsex.trap)==intrest.covar)])), length = dim(my.data)[1] ))
pred.bam.st.habseason.hab.density.2.0.dum = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "terms")

dum.se = as.numeric(pred.bam.st.habseason.hab.density.2.0.dum$se.fit[,which(attributes(pred.bam.st.habseason.hab.density.2.0.dum$fit)$dimnames[[2]] == intrest.term)])
pred.bam.st.habseason.hab.density.2.0 = pred.bam.st.habseason.hab.density.2.0.dum

my.data$fit = as.numeric( pred.bam.st.habseason.hab.density.2.0.dum$fit[,which(attributes(pred.bam.st.habseason.hab.density.2.0.dum$fit)$dimnames[[2]] == intrest.term)])
my.data$upper.95 = my.data$fit - (dum.se*1.96)
my.data$lower.95 = my.data$fit + (dum.se*1.96)
#my.data$expl = exp(my.data[,which(names(my.data)==intrest.covar)]  )-1
my.data$expl = my.data[,which(names(my.data)==intrest.covar)]  

x.max = 300 #max(my.data$expl)
breaks.by = 100
break.divider = 1
rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])
ylims = c(0 , round(max(my.data$upper.95), digits = 3)+0.001)




#main.plot = 
    ggplot()+
    geom_ribbon(aes(x = my.data$expl, ymin = my.data$lower.95, ymax = my.data$upper.95), fill = "grey80"#, alpha = 0.4
    )+
    geom_line(aes(y = my.data$fit, x = my.data$expl), size = 1, col = "black")+
    labs(x = expl.lab, y = "Effect size") +
    geom_rug( aes(x = rugs),sides = "b",# position = "jitter",
              size = 0.01) +
    scale_y_continuous() +
        scale_x_continuous(breaks=seq(0,x.max, by = breaks.by ), labels= seq(0,max(my.data$expl), by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
        theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )

  ##########################################################################################
  # main effect of pasture.mk2----
  intrest.covar = "pasture.mk2"
  expl.lab =  'Pasture area within 1.5km ('~km^2*')'
  intrest.term = "ti(pasture.mk2)"
  
  interested.length = the.interested.length
  divider = 1000000
  
  #min.expl = min(unsex.trap[,names(unsex.trap)==intrest.covar ])
  #max.expl = max(unsex.trap[,names(unsex.trap)==intrest.covar ])
  #max.expl = 250 # 2 traps over 250m, 5 sessions
  min.expl = min(unsex.trap[,names(unsex.trap)==intrest.covar ])
  max.expl = max(unsex.trap[,names(unsex.trap)==intrest.covar ])
  
  my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                       interested.covar = intrest.covar, interested.length = interested.length, 
                                       by.fac = F, interesting.factor.name = NULL)
  my.data = my.data1
  
  pred.bam.st.habseason.hab.density.2.0.dum = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "terms")
  
  dum.se = as.numeric(pred.bam.st.habseason.hab.density.2.0.dum$se.fit[,which(attributes(pred.bam.st.habseason.hab.density.2.0.dum$fit)$dimnames[[2]] == intrest.term)])
  pred.bam.st.habseason.hab.density.2.0 = pred.bam.st.habseason.hab.density.2.0.dum
  
  my.data$fit = as.numeric( pred.bam.st.habseason.hab.density.2.0.dum$fit[,which(attributes(pred.bam.st.habseason.hab.density.2.0.dum$fit)$dimnames[[2]] == intrest.term)])
  my.data$upper.95 = my.data$fit + (dum.se*1.96)
  my.data$lower.95 = my.data$fit - (dum.se*1.96)
  #my.data$expl = exp(my.data[,which(names(my.data)==intrest.covar)]  )-1
  my.data$expl = my.data[,which(names(my.data)==intrest.covar)]  
 
  
  x.max = max(my.data$expl)
  breaks.by = 1000000
  break.divider = 1000000
  rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])
  ylims = c(0 , round(max(my.data$upper.95), digits = 3)+0.001)
  
  
  
  
  #main.plot = 
  ggplot()+
      geom_ribbon(aes(x = my.data$expl, ymin = my.data$lower.95, ymax = my.data$upper.95), fill = "grey80"#, alpha = 0.4
      )+
      geom_line(aes(y = my.data$fit, x = my.data$expl), size = 1, col = "black")+
      labs(x = expl.lab, y = "Effect size") +
      geom_rug( aes(x = rugs),sides = "b",# position = "jitter",
                size = 0.01) +
      scale_y_continuous() +
      scale_x_continuous(breaks=seq(0,x.max, by = breaks.by ), labels= seq(0,x.max, by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
      theme_bw() +
      theme(plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color = 'black'),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12)
      )
  
  ##########################################################################################
  # main effect of rough.grass.mk2----
  intrest.covar = "rough.grass.mk2"
  expl.lab =  'Rough grass area within 1.5km ('~km^2*')'
  intrest.term = "ti(rough.grass.mk2)"
  
  interested.length = the.interested.length
  divider = 1000000
  
  #min.expl = min(unsex.trap[,names(unsex.trap)==intrest.covar ])
  #max.expl = max(unsex.trap[,names(unsex.trap)==intrest.covar ])
  #max.expl = 250 # 2 traps over 250m, 5 sessions
  min.expl = min(unsex.trap[,names(unsex.trap)==intrest.covar ])
  max.expl = max(unsex.trap[,names(unsex.trap)==intrest.covar ])
  
  my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                       interested.covar = intrest.covar, interested.length = interested.length, 
                                       by.fac = F, interesting.factor.name = NULL)
  my.data = my.data1
  
  pred.bam.st.habseason.hab.density.2.0.dum = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "terms")
  
  dum.se = as.numeric(pred.bam.st.habseason.hab.density.2.0.dum$se.fit[,which(attributes(pred.bam.st.habseason.hab.density.2.0.dum$fit)$dimnames[[2]] == intrest.term)])
  pred.bam.st.habseason.hab.density.2.0 = pred.bam.st.habseason.hab.density.2.0.dum
  
  my.data$fit = as.numeric( pred.bam.st.habseason.hab.density.2.0.dum$fit[,which(attributes(pred.bam.st.habseason.hab.density.2.0.dum$fit)$dimnames[[2]] == intrest.term)])
  my.data$upper.95 = my.data$fit + (dum.se*1.96)
  my.data$lower.95 = my.data$fit - (dum.se*1.96)
  #my.data$expl = exp(my.data[,which(names(my.data)==intrest.covar)]  )-1
  my.data$expl = my.data[,which(names(my.data)==intrest.covar)]  
  
  rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])
  
  x.max = max(my.data$expl)
  breaks.by = 1000000
  break.divider = 1000000
  rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])
  ylims = c(0 , round(max(my.data$upper.95), digits = 3)+0.001)
  
  
  
  
  #main.plot = 
  ggplot()+
      geom_ribbon(aes(x = my.data$expl, ymin = my.data$lower.95, ymax = my.data$upper.95), fill = "grey80"#, alpha = 0.4
      )+
      geom_line(aes(y = my.data$fit, x = my.data$expl), size = 1, col = "black")+
      labs(x = expl.lab, y = "Effect size") +
      geom_rug( aes(x = rugs),sides = "b",# position = "jitter",
                size = 0.01) +
      scale_y_continuous() +
      scale_x_continuous(breaks=seq(0,x.max, by = breaks.by ), labels= seq(0,x.max, by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
      theme_bw() +
      theme(plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color = 'black'),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12)
      )
  
  
##########################################################################################
  # main effect of water.edge.length----
  
  
  intrest.covar = "water.edge.length"
  intrest.term = "ti(water.edge.length)"
  expl.lab =   "Water edge length (km)"
  interested.length = the.interested.length
  
  my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                       interested.covar = intrest.covar, interested.length = interested.length, 
                                       by.fac = F, interesting.factor.name = NULL)
  my.data = my.data1
  #my.data[,which(names(my.data)==intrest.covar)]  = log(1+ seq(0, exp(max(unsex.trap[,which(names(unsex.trap)==intrest.covar)])), length = dim(my.data)[1] ))
  pred.bam.st.habseason.hab.density.2.0.dum = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "terms")
  
  dum.se = as.numeric(pred.bam.st.habseason.hab.density.2.0.dum$se.fit[,which(attributes(pred.bam.st.habseason.hab.density.2.0.dum$fit)$dimnames[[2]] == intrest.term)])
  pred.bam.st.habseason.hab.density.2.0 = pred.bam.st.habseason.hab.density.2.0.dum
  
  my.data$fit = as.numeric( pred.bam.st.habseason.hab.density.2.0.dum$fit[,which(attributes(pred.bam.st.habseason.hab.density.2.0.dum$fit)$dimnames[[2]] == intrest.term)])
  my.data$upper.95 = my.data$fit - (dum.se*1.96)
  my.data$lower.95 = my.data$fit + (dum.se*1.96)
  #my.data$expl = exp(my.data[,which(names(my.data)==intrest.covar)]  )-1
  my.data$expl = my.data[,which(names(my.data)==intrest.covar)]  
  
  x.max = max(my.data$expl)
  breaks.by = 10000
  break.divider = 1000
  rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])
  ylims = c(0 , round(max(my.data$upper.95), digits = 3)+0.001)
  
  
  #main.plot = 
  ggplot()+
      geom_ribbon(aes(x = my.data$expl, ymin = my.data$lower.95, ymax = my.data$upper.95), fill = "grey80"#, alpha = 0.4
      )+
      geom_line(aes(y = my.data$fit, x = my.data$expl), size = 1, col = "black")+
      labs(x = expl.lab, y = "Effect size") +
      geom_rug( aes(x = rugs),sides = "b",# position = "jitter",
                size = 0.01) +
      scale_y_continuous() +
      scale_x_continuous(breaks=seq(0,x.max, by = breaks.by ), labels= seq(0,x.max, by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
      theme_bw() +
      theme(plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color = 'black'),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12)
      )
  
  
  
##########################################################################################
#
#
##########################################################################################
# TENSOR OF SEASONALITIES -----
##########################################################################################
  # tensor effect of altitude ----
  intrest.covar = "ALT3"
  expl.lab =   "Altitude (m)"
  intrest.term = "ti(julian,ALT3)"
  itneracting = "julian"
  breaks.by = 100
  break.divider = 1
  interacting.lab = "Day of year"
  dist.too.far = 0.05
  divider = 1
  
  interested.length = the.interested.length
  min.inter = min(unsex.trap[,names(unsex.trap)==itneracting ])
  max.inter = max(unsex.trap[,names(unsex.trap)==itneracting ])
  min.expl = min(unsex.trap[,names(unsex.trap)==intrest.covar ])
  max.expl = max(unsex.trap[,names(unsex.trap)==intrest.covar ])
  
  my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap, interested.covar = intrest.covar, 
                                       interested.length = interested.length, by.fac = F, interesting.factor.name = NULL)
  # my.data1$exp.expl = seq(min.expl, max.expl, length = dim(my.data1)[1])
  # my.data1[, which(names(my.data1)== intrest.covar)] = my.data1$exp.expl
  my.data1[, which(names(my.data1) %in% c(intrest.covar ))] = seq(min.expl, max.expl, length = dim(my.data1)[1])
  
  my.data = my.data1
  my.data = my.data1[rep(1:(dim(my.data1)[1]), times = interested.length ),]
  my.data[, which(names(my.data)== itneracting)] = rep (seq(min.inter, max.inter, length = interested.length), each = interested.length)
  
  #my.data[,which(names(my.data)==intrest.covar)]  = log(1+ seq(0, exp(max(unsex.trap[,which(names(unsex.trap)==intrest.covar)])), length = dim(my.data)[1] ))
  pred.bam.st.habseason.hab.density.2.0.dum = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "terms")
  
  pred.bam.st.habseason.hab.density.2.0 = pred.bam.st.habseason.hab.density.2.0.dum
  
  my.data$fit = as.numeric( pred.bam.st.habseason.hab.density.2.0.dum$fit[,which(attributes(pred.bam.st.habseason.hab.density.2.0.dum$fit)$dimnames[[2]] == intrest.term)])
  my.data$expl = my.data[,which(names(my.data)==intrest.covar)]  
  x.max = 300#max(my.data$expl)
  
  rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])
  ylims = c(0 , round(max(my.data$upper.95), digits = 3)+0.001)
  
  plotr = my.data[, which(names(my.data) %in% c(intrest.covar, itneracting ))]
  plotr$fit = my.data$fit
  
  plotr$expl = as.numeric(my.data[,which(names(my.data) %in% c(intrest.covar))]/divider)
  plotr$intr = as.numeric(my.data[,which(names(my.data) == (itneracting))])
  
  og.dat = unsex.trap[unsex.trap[,which(names(unsex.trap) %in% c( intrest.covar))]<max.expl,]
  og.dat$intr = og.dat[,which(names(og.dat) %in% c( itneracting))]
  og.dat$expl = og.dat[,which(names(og.dat) %in% c( intrest.covar))]
  rug = unique(og.dat[,which(names(og.dat) == intrest.covar)])/divider
  
  caps= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c(intrest.covar, itneracting))]
  caps$itneracting= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c( itneracting))]
  caps$intrest.covar= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c( intrest.covar))]/divider
  
  vars.too.far<- exclude.too.far(plotr$expl , plotr$intr,
                                 og.dat$expl/divider , og.dat$intr,
                                 dist.too.far)
  cap.point.size = 0.01
  
  # # # #
  #col.scale
  if(abs(range(plotr$fit[!vars.too.far])[2])> abs(range(plotr$fit[!vars.too.far])[1])){
      col.scale.end = 1
      col.scale.start = 0.5 - (abs(range(plotr$fit[!vars.too.far])[1])/abs(range(plotr$fit[!vars.too.far])[2]))
  } else{
      col.scale.start = 0
      col.scale.end = 0.5 + (abs(range(plotr$fit[!vars.too.far])[2])/abs(range(plotr$fit[!vars.too.far])[1]))
  }
  # # # #
  ggplot() + 
      geom_tile(aes(plotr$expl[!vars.too.far], plotr$intr[!vars.too.far], fill = plotr$fit[!vars.too.far])) +
      #geom_point(aes(caps$intrest.covar, caps$itneracting), size = cap.point.size, color = "white")+
      # geom_rug(aes(x = rug), alpha = 1/2,sides = "b",size = 0.2)+
      labs(y = interacting.lab, x = expl.lab , z = "Target Sex Capture Probability") +
      geom_contour(aes(plotr$expl[!vars.too.far], plotr$intr[!vars.too.far], z = plotr$fit[!vars.too.far]), color = "grey")+
      #  scale_fill_viridis(begin = col.scale.start, end = col.scale.end)+
      scale_fill_viridis(name = "Effect size") +
      scale_y_continuous(breaks= c(1/12, 3/12, 5/12, 7/12, 9/12,11/12)*365, labels= c("Feb", "Apr", "June", "Aug","Oct", "Dec"))+
      scale_x_continuous(breaks=seq(0,x.max, by = breaks.by ), labels= seq(0,x.max, by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
      theme_bw() +
      theme(plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color = 'black'),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12)
      )
  # tensor effect of Pasture ----
  intrest.covar = "pasture.mk2"
  expl.lab =  'Pasture area within 1.5km ('~km^2*')'
  intrest.term = "ti(julian,pasture.mk2)"
  itneracting = "julian"
  breaks.by = 1000000
  break.divider = 1000000
  interacting.lab = "Day of year"
  dist.too.far = 0.05
  divider = 1
  
  interested.length = the.interested.length
  min.inter = min(unsex.trap[,names(unsex.trap)==itneracting ])
  max.inter = max(unsex.trap[,names(unsex.trap)==itneracting ])
  min.expl = min(unsex.trap[,names(unsex.trap)==intrest.covar ])
  max.expl = max(unsex.trap[,names(unsex.trap)==intrest.covar ])
  
  my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap, interested.covar = intrest.covar, 
                                       interested.length = interested.length, by.fac = F, interesting.factor.name = NULL)
  # my.data1$exp.expl = seq(min.expl, max.expl, length = dim(my.data1)[1])
  # my.data1[, which(names(my.data1)== intrest.covar)] = my.data1$exp.expl
  my.data1[, which(names(my.data1) %in% c(intrest.covar ))] = seq(min.expl, max.expl, length = dim(my.data1)[1])
  
  
  my.data = my.data1
  my.data = my.data1[rep(1:(dim(my.data1)[1]), times = interested.length ),]
  my.data[, which(names(my.data)== itneracting)] = rep (seq(min.inter, max.inter, length = interested.length), each = interested.length)
  
  
  
  #my.data[,which(names(my.data)==intrest.covar)]  = log(1+ seq(0, exp(max(unsex.trap[,which(names(unsex.trap)==intrest.covar)])), length = dim(my.data)[1] ))
  pred.bam.st.habseason.hab.density.2.0.dum = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "terms")
  
  pred.bam.st.habseason.hab.density.2.0 = pred.bam.st.habseason.hab.density.2.0.dum
  
  my.data$fit = as.numeric( pred.bam.st.habseason.hab.density.2.0.dum$fit[,which(attributes(pred.bam.st.habseason.hab.density.2.0.dum$fit)$dimnames[[2]] == intrest.term)])
  my.data$expl = my.data[,which(names(my.data)==intrest.covar)]  
  x.max = max(my.data$expl)
  
  rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])
  ylims = c(0 , round(max(my.data$upper.95), digits = 3)+0.001)
  
  
  plotr = my.data[, which(names(my.data) %in% c(intrest.covar, itneracting ))]
  plotr$fit = my.data$fit
  
  plotr$expl = as.numeric(my.data[,which(names(my.data) %in% c(intrest.covar))]/divider)
  plotr$intr = as.numeric(my.data[,which(names(my.data) == (itneracting))])
  
  og.dat = unsex.trap[unsex.trap[,which(names(unsex.trap) %in% c( intrest.covar))]<max.expl,]
  og.dat$intr = og.dat[,which(names(og.dat) %in% c( itneracting))]
  og.dat$expl = og.dat[,which(names(og.dat) %in% c( intrest.covar))]
  rug = unique(og.dat[,which(names(og.dat) == intrest.covar)])/divider
  
  caps= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c(intrest.covar, itneracting))]
  caps$itneracting= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c( itneracting))]
  caps$intrest.covar= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c( intrest.covar))]/divider
  
  vars.too.far<- exclude.too.far(plotr$expl , plotr$intr,
                                 og.dat$expl/divider , og.dat$intr,
                                 dist.too.far)
  cap.point.size = 0.01
  
  # # # #
  #col.scale
  if(abs(range(plotr$fit[!vars.too.far])[2])> abs(range(plotr$fit[!vars.too.far])[1])){
      col.scale.end = 1
      col.scale.start = 0.5 - (abs(range(plotr$fit[!vars.too.far])[1])/abs(range(plotr$fit[!vars.too.far])[2]))
  } else{
      col.scale.start = 0
      col.scale.end = 0.5 + (abs(range(plotr$fit[!vars.too.far])[2])/abs(range(plotr$fit[!vars.too.far])[1]))
  }
  # # # #
  ggplot() + 
      geom_tile(aes(plotr$expl[!vars.too.far], plotr$intr[!vars.too.far], fill = plotr$fit[!vars.too.far])) +
      #geom_point(aes(caps$intrest.covar, caps$itneracting), size = cap.point.size, color = "white")+
      # geom_rug(aes(x = rug), alpha = 1/2,sides = "b",size = 0.2)+
      labs(y = interacting.lab, x = expl.lab , z = "Target Sex Capture Probability") +
      geom_contour(aes(plotr$expl[!vars.too.far], plotr$intr[!vars.too.far], z = plotr$fit[!vars.too.far]), color = "grey")+
      #  scale_fill_viridis(begin = col.scale.start, end = col.scale.end)+
      scale_fill_viridis(name = "Effect size") +
      scale_y_continuous(breaks= c(1/12, 3/12, 5/12, 7/12, 9/12,11/12)*365, labels= c("Feb", "Apr", "June", "Aug","Oct", "Dec"))+
      scale_x_continuous(breaks=seq(0,x.max, by = breaks.by ), labels= seq(0,x.max, by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
      theme_bw() +
      theme(plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color = 'black'),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12)
      )
  # tensor effect of rough grass ----
  intrest.covar = "rough.grass.mk2"
  expl.lab =  'Rough grass area within 1.5km ('~km^2*')'
  intrest.term = "ti(julian,rough.grass.mk2)"
  itneracting = "julian"
  breaks.by = 1000000
  break.divider = 1000000
  interacting.lab = "Day of year"
  dist.too.far = 0.05
  divider = 1
  
  interested.length = the.interested.length
  min.inter = min(unsex.trap[,names(unsex.trap)==itneracting ])
  max.inter = max(unsex.trap[,names(unsex.trap)==itneracting ])
  min.expl = min(unsex.trap[,names(unsex.trap)==intrest.covar ])
  max.expl = max(unsex.trap[,names(unsex.trap)==intrest.covar ])
  
  my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap, interested.covar = intrest.covar, 
                                       interested.length = interested.length, by.fac = F, interesting.factor.name = NULL)
  # my.data1$exp.expl = seq(min.expl, max.expl, length = dim(my.data1)[1])
  # my.data1[, which(names(my.data1)== intrest.covar)] = my.data1$exp.expl
  my.data1[, which(names(my.data1) %in% c(intrest.covar ))] = seq(min.expl, max.expl, length = dim(my.data1)[1])
  
  
  my.data = my.data1
  my.data = my.data1[rep(1:(dim(my.data1)[1]), times = interested.length ),]
  my.data[, which(names(my.data)== itneracting)] = rep (seq(min.inter, max.inter, length = interested.length), each = interested.length)
  
  
  
  #my.data[,which(names(my.data)==intrest.covar)]  = log(1+ seq(0, exp(max(unsex.trap[,which(names(unsex.trap)==intrest.covar)])), length = dim(my.data)[1] ))
  pred.bam.st.habseason.hab.density.2.0.dum = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "terms")
  
  pred.bam.st.habseason.hab.density.2.0 = pred.bam.st.habseason.hab.density.2.0.dum
  
  my.data$fit = as.numeric( pred.bam.st.habseason.hab.density.2.0.dum$fit[,which(attributes(pred.bam.st.habseason.hab.density.2.0.dum$fit)$dimnames[[2]] == intrest.term)])
  my.data$expl = my.data[,which(names(my.data)==intrest.covar)]  
  x.max = max(my.data$expl)
  
  rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])
  ylims = c(0 , round(max(my.data$upper.95), digits = 3)+0.001)
  
  
  plotr = my.data[, which(names(my.data) %in% c(intrest.covar, itneracting ))]
  plotr$fit = my.data$fit
  
  plotr$expl = as.numeric(my.data[,which(names(my.data) %in% c(intrest.covar))]/divider)
  plotr$intr = as.numeric(my.data[,which(names(my.data) == (itneracting))])
  
  og.dat = unsex.trap[unsex.trap[,which(names(unsex.trap) %in% c( intrest.covar))]<max.expl,]
  og.dat$intr = og.dat[,which(names(og.dat) %in% c( itneracting))]
  og.dat$expl = og.dat[,which(names(og.dat) %in% c( intrest.covar))]
  rug = unique(og.dat[,which(names(og.dat) == intrest.covar)])/divider
  
  caps= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c(intrest.covar, itneracting))]
  caps$itneracting= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c( itneracting))]
  caps$intrest.covar= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c( intrest.covar))]/divider
  
  vars.too.far<- exclude.too.far(plotr$expl , plotr$intr,
                                 og.dat$expl/divider , og.dat$intr,
                                 dist.too.far)
  cap.point.size = 0.01
  
  # # # #
  #col.scale
  if(abs(range(plotr$fit[!vars.too.far])[2])> abs(range(plotr$fit[!vars.too.far])[1])){
      col.scale.end = 1
      col.scale.start = 0.5 - (abs(range(plotr$fit[!vars.too.far])[1])/abs(range(plotr$fit[!vars.too.far])[2]))
  } else{
      col.scale.start = 0
      col.scale.end = 0.5 + (abs(range(plotr$fit[!vars.too.far])[2])/abs(range(plotr$fit[!vars.too.far])[1]))
  }
  # # # #
  ggplot() + 
      geom_tile(aes(plotr$expl[!vars.too.far], plotr$intr[!vars.too.far], fill = plotr$fit[!vars.too.far])) +
      #geom_point(aes(caps$intrest.covar, caps$itneracting), size = cap.point.size, color = "white")+
      # geom_rug(aes(x = rug), alpha = 1/2,sides = "b",size = 0.2)+
      labs(y = interacting.lab, x = expl.lab , z = "Target Sex Capture Probability") +
      geom_contour(aes(plotr$expl[!vars.too.far], plotr$intr[!vars.too.far], z = plotr$fit[!vars.too.far]), color = "grey")+
      #  scale_fill_viridis(begin = col.scale.start, end = col.scale.end)+
      scale_fill_viridis(name = "Effect size") +
      scale_y_continuous(breaks= c(1/12, 3/12, 5/12, 7/12, 9/12,11/12)*365, labels= c("Feb", "Apr", "June", "Aug","Oct", "Dec"))+
      scale_x_continuous(breaks=seq(0,x.max, by = breaks.by ), labels= seq(0,x.max, by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
      theme_bw() +
      theme(plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color = 'black'),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12)
      )
  # tensor effect of Water edge ----
  intrest.covar = "water.edge.length"
  expl.lab =  'Water edge within 1.5km km'
  intrest.term = "ti(julian,water.edge.length)"
  itneracting = "julian"
  breaks.by = 10000
  break.divider = 1000
  interacting.lab = "Day of year"
  dist.too.far = 0.05
  divider = 1
  
  interested.length = the.interested.length
  min.inter = min(unsex.trap[,names(unsex.trap)==itneracting ])
  max.inter = max(unsex.trap[,names(unsex.trap)==itneracting ])
  min.expl = min(unsex.trap[,names(unsex.trap)==intrest.covar ])
  max.expl = max(unsex.trap[,names(unsex.trap)==intrest.covar ])
  
  my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap, interested.covar = intrest.covar, 
                                       interested.length = interested.length, by.fac = F, interesting.factor.name = NULL)
  # my.data1$exp.expl = seq(min.expl, max.expl, length = dim(my.data1)[1])
  # my.data1[, which(names(my.data1)== intrest.covar)] = my.data1$exp.expl
  my.data1[, which(names(my.data1) %in% c(intrest.covar ))] = seq(min.expl, max.expl, length = dim(my.data1)[1])
  
  
  my.data = my.data1
  my.data = my.data1[rep(1:(dim(my.data1)[1]), times = interested.length ),]
  my.data[, which(names(my.data)== itneracting)] = rep (seq(min.inter, max.inter, length = interested.length), each = interested.length)
  
  
  
  #my.data[,which(names(my.data)==intrest.covar)]  = log(1+ seq(0, exp(max(unsex.trap[,which(names(unsex.trap)==intrest.covar)])), length = dim(my.data)[1] ))
  pred.bam.st.habseason.hab.density.2.0.dum = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "terms")
  
  pred.bam.st.habseason.hab.density.2.0 = pred.bam.st.habseason.hab.density.2.0.dum
  
  my.data$fit = as.numeric( pred.bam.st.habseason.hab.density.2.0.dum$fit[,which(attributes(pred.bam.st.habseason.hab.density.2.0.dum$fit)$dimnames[[2]] == intrest.term)])
  my.data$expl = my.data[,which(names(my.data)==intrest.covar)]  
  x.max = max(my.data$expl)
  
  rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])
  ylims = c(0 , round(max(my.data$upper.95), digits = 3)+0.001)
  
  
  plotr = my.data[, which(names(my.data) %in% c(intrest.covar, itneracting ))]
  plotr$fit = my.data$fit
  
  plotr$expl = as.numeric(my.data[,which(names(my.data) %in% c(intrest.covar))]/divider)
  plotr$intr = as.numeric(my.data[,which(names(my.data) == (itneracting))])
  
  og.dat = unsex.trap[unsex.trap[,which(names(unsex.trap) %in% c( intrest.covar))]<max.expl,]
  og.dat$intr = og.dat[,which(names(og.dat) %in% c( itneracting))]
  og.dat$expl = og.dat[,which(names(og.dat) %in% c( intrest.covar))]
  rug = unique(og.dat[,which(names(og.dat) == intrest.covar)])/divider
  
  caps= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c(intrest.covar, itneracting))]
  caps$itneracting= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c( itneracting))]
  caps$intrest.covar= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c( intrest.covar))]/divider
  
  vars.too.far<- exclude.too.far(plotr$expl , plotr$intr,
                                 og.dat$expl/divider , og.dat$intr,
                                 dist.too.far)
  cap.point.size = 0.01
  
  # # # #
  #col.scale
  if(abs(range(plotr$fit[!vars.too.far])[2])> abs(range(plotr$fit[!vars.too.far])[1])){
      col.scale.end = 1
      col.scale.start = 0.5 - (abs(range(plotr$fit[!vars.too.far])[1])/abs(range(plotr$fit[!vars.too.far])[2]))
  } else{
      col.scale.start = 0
      col.scale.end = 0.5 + (abs(range(plotr$fit[!vars.too.far])[2])/abs(range(plotr$fit[!vars.too.far])[1]))
  }
  # # # #
  ggplot() + 
      geom_tile(aes(plotr$expl[!vars.too.far], plotr$intr[!vars.too.far], fill = plotr$fit[!vars.too.far])) +
      #geom_point(aes(caps$intrest.covar, caps$itneracting), size = cap.point.size, color = "white")+
      # geom_rug(aes(x = rug), alpha = 1/2,sides = "b",size = 0.2)+
      labs(y = interacting.lab, x = expl.lab , z = "Target Sex Capture Probability") +
      geom_contour(aes(plotr$expl[!vars.too.far], plotr$intr[!vars.too.far], z = plotr$fit[!vars.too.far]), color = "grey")+
      #  scale_fill_viridis(begin = col.scale.start, end = col.scale.end)+
      scale_fill_viridis(name = "Effect size") +
      scale_y_continuous(breaks= c(1/12, 3/12, 5/12, 7/12, 9/12,11/12)*365, labels= c("Feb", "Apr", "June", "Aug","Oct", "Dec"))+
      scale_x_continuous(breaks=seq(0,x.max, by = breaks.by ), labels= seq(0,x.max, by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
      theme_bw() +
      theme(plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color = 'black'),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12)
      )
  ##########################################################################################

##########################################################################################
#  SEASONAL TENSOR PREDICTIONS ----
##########################################################################################
# predicted tensor effect of altitude ----
intrest.covar = "ALT3"
itneracting = "julian"
expl.lab =   "Altitude (m)"
interacting.lab = "Day of year"
interested.length = the.interested.length
divider = 1
dist.too.far = 0.05
max.expl = 300 # max(unsex.trap[,names(unsex.trap)==intrest.covar ]) # 2 traps over 250m, 5 sessions

#min.expl = exp(min(unsex.trap[,names(unsex.trap)==intrest.covar ]))-1
#max.expl = exp(max(unsex.trap[,names(unsex.trap)==intrest.covar ]))-1
min.expl = min(unsex.trap[,names(unsex.trap)==intrest.covar ])
#max.expl = max(unsex.trap[,names(unsex.trap)==intrest.covar ])

min.inter = min(unsex.trap[,names(unsex.trap)==itneracting ])
max.inter = max(unsex.trap[,names(unsex.trap)==itneracting ])


my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap, interested.covar = intrest.covar, 
                                     interested.length = interested.length, by.fac = F, interesting.factor.name = NULL)
# my.data1$exp.expl = seq(min.expl, max.expl, length = dim(my.data1)[1])
# my.data1[, which(names(my.data1)== intrest.covar)] = my.data1$exp.expl
my.data1[, which(names(my.data1) %in% c(intrest.covar ))] = seq(min.expl, max.expl, length = dim(my.data1)[1])


my.data = my.data1
my.data = my.data1[rep(1:(dim(my.data1)[1]), times = interested.length ),]
my.data[, which(names(my.data)== itneracting)] = rep (seq(min.inter, max.inter, length = interested.length), each = interested.length)

pred.bam.st.habseason.hab.density.2.0 =  predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")

my.data$fit = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit))
my.data$upper.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit + 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))
my.data$lower.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit - 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))
my.data$expl = my.data[,which(names(my.data)==intrest.covar)]
my.data$itneracting = my.data[,which(names(my.data)==itneracting)]

rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])

ylims = c(0 , round(max(my.data$upper.95), digits = 3)+0.001)

plotr = my.data[, which(names(my.data) %in% c(intrest.covar, itneracting ))]
plotr$fit = my.data$fit
plotr$upper.95 = my.data$upper.95
plotr$lower.95 = my.data$lower.95

plotr$expl = as.numeric(my.data[,which(names(my.data) %in% c(intrest.covar))]/divider)
plotr$intr = as.numeric(my.data[,which(names(my.data) == (itneracting))])

og.dat = unsex.trap[unsex.trap[,which(names(unsex.trap) %in% c( intrest.covar))]<max.expl,]
og.dat$intr = og.dat[,which(names(og.dat) %in% c( itneracting))]
og.dat$expl = og.dat[,which(names(og.dat) %in% c( intrest.covar))]
rug = unique(og.dat[,which(names(og.dat) == intrest.covar)])/divider

caps= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c(intrest.covar, itneracting))]
caps$itneracting= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c( itneracting))]
caps$intrest.covar= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c( intrest.covar))]/divider

vars.too.far<- exclude.too.far(plotr$expl , plotr$intr,
                               og.dat$expl , og.dat$intr,
                               dist.too.far)
cap.point.size = 0.05

# # # #
#col.scale
if(abs(range(plotr$fit[!vars.too.far])[2])> abs(range(plotr$fit[!vars.too.far])[1])){
    col.scale.end = 1
    col.scale.start = 0.5 - (abs(range(plotr$fit[!vars.too.far])[1])/abs(range(plotr$fit[!vars.too.far])[2]))
} else{
    col.scale.start = 0
    col.scale.end = 0.5 + (abs(range(plotr$fit[!vars.too.far])[2])/abs(range(plotr$fit[!vars.too.far])[1]))
}
# # # #
ggplot() + 
    geom_tile(aes(plotr$expl[!vars.too.far], plotr$intr[!vars.too.far], fill = plotr$fit[!vars.too.far])) +
    # geom_rug(aes(x = rug), alpha = 1/2,sides = "b",size = 0.2)+
    labs(y = interacting.lab, x = expl.lab , z = "Target Sex Capture Probability") +
    geom_contour(aes(plotr$expl[!vars.too.far], plotr$intr[!vars.too.far], z = plotr$fit[!vars.too.far]), color = "grey")+
    #  scale_fill_viridis(begin = col.scale.start, end = col.scale.end)+
    scale_fill_viridis(name = "Trapping rate") +
    scale_y_continuous(breaks= c(1/12, 3/12, 5/12, 7/12, 9/12,11/12)*365, labels= c("Feb", "Apr", "June", "Aug","Oct", "Dec"))+
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )

# predicted effect of pasture area ----
intrest.covar = "pasture.mk2"
itneracting = "julian"
expl.lab =   'Pasture area within 1.5km ('~km^2*')'
interacting.lab = "Day of year"
interested.length = the.interested.length
divider = 1000000
dist.too.far = 0.05

#min.expl = exp(min(unsex.trap[,names(unsex.trap)==intrest.covar ]))-1
#max.expl = exp(max(unsex.trap[,names(unsex.trap)==intrest.covar ]))-1
min.expl = min(unsex.trap[,names(unsex.trap)==intrest.covar ])
max.expl = max(unsex.trap[,names(unsex.trap)==intrest.covar ])
#max.expl = 250 # 2 traps over 250m, 5 sessions

min.inter = min(unsex.trap[,names(unsex.trap)==itneracting ])
max.inter = max(unsex.trap[,names(unsex.trap)==itneracting ])

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap, interested.covar = intrest.covar, 
                                     interested.length = interested.length, by.fac = F, interesting.factor.name = NULL)
# my.data1$exp.expl = seq(min.expl, max.expl, length = dim(my.data1)[1])
# my.data1[, which(names(my.data1)== intrest.covar)] = my.data1$exp.expl
my.data1[, which(names(my.data1) %in% c(intrest.covar ))] = seq(min.expl, max.expl, length = dim(my.data1)[1])

my.data = my.data1
my.data = my.data1[rep(1:(dim(my.data1)[1]), times = interested.length ),]
my.data[, which(names(my.data)== itneracting)] = rep (seq(min.inter, max.inter, length = interested.length), each = interested.length)


pred.bam.st.habseason.hab.density.2.0 =  predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")

my.data$fit = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit))
my.data$upper.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit + 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))
my.data$lower.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit - 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))

my.data$expl = my.data[,which(names(my.data)==intrest.covar)]
my.data$itneracting = my.data[,which(names(my.data)==itneracting)]

rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])

ylims = c(0 , round(max(my.data$upper.95), digits = 3)+0.001)

plotr = my.data[, which(names(my.data) %in% c(intrest.covar, itneracting ))]
plotr$fit = my.data$fit
plotr$upper.95 = my.data$upper.95
plotr$lower.95 = my.data$lower.95

plotr$expl = as.numeric(my.data[,which(names(my.data) %in% c(intrest.covar))]/divider)
plotr$intr = as.numeric(my.data[,which(names(my.data) == (itneracting))])

og.dat = unsex.trap[unsex.trap[,which(names(unsex.trap) %in% c( intrest.covar))]<max.expl,]
og.dat$intr = og.dat[,which(names(og.dat) %in% c( itneracting))]
og.dat$expl = og.dat[,which(names(og.dat) %in% c( intrest.covar))]/divider
rug = unique(og.dat[,which(names(og.dat) == intrest.covar)])/divider

caps= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c(intrest.covar, itneracting))]
caps$itneracting= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c( itneracting))]
caps$intrest.covar= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c( intrest.covar))]/divider

vars.too.far<- exclude.too.far(plotr$expl , plotr$intr,
                               og.dat$expl , og.dat$intr,
                               dist.too.far)
cap.point.size = 0.05

# # # #
#col.scale
if(abs(range(plotr$fit[!vars.too.far])[2])> abs(range(plotr$fit[!vars.too.far])[1])){
    col.scale.end = 1
    col.scale.start = 0.5 - (abs(range(plotr$fit[!vars.too.far])[1])/abs(range(plotr$fit[!vars.too.far])[2]))
} else{
    col.scale.start = 0
    col.scale.end = 0.5 + (abs(range(plotr$fit[!vars.too.far])[2])/abs(range(plotr$fit[!vars.too.far])[1]))
}
# # # #
ggplot() + 
    geom_tile(aes(plotr$expl[!vars.too.far], plotr$intr[!vars.too.far], fill = plotr$fit[!vars.too.far])) +
    # geom_rug(aes(x = rug), alpha = 1/2,sides = "b",size = 0.2)+
    labs(y = interacting.lab, x = expl.lab , z = "Target Sex Capture Probability") +
    geom_contour(aes(plotr$expl[!vars.too.far], plotr$intr[!vars.too.far], z = plotr$fit[!vars.too.far]), color = "grey")+
    #  scale_fill_viridis(begin = col.scale.start, end = col.scale.end)+
    scale_fill_viridis(name = "Trapping rate") +
    scale_y_continuous(breaks= c(1/12, 3/12, 5/12, 7/12, 9/12,11/12)*365, labels= c("Feb", "Apr", "June", "Aug","Oct", "Dec"))+
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )
# date of peak at low bog
my.data$itneracting[my.data$expl==0][my.data$fit[my.data$expl==0]==max( my.data$fit[my.data$expl==0]) ] # 1st sept = day 244
#lowest
my.data[my.data$fit == min(my.data$fit ),] # 1st march and 3.499km2

# predicted effect of rough grass area ----
intrest.covar = "rough.grass.mk2"
itneracting = "julian"
expl.lab =   'Rough grass area within 1.5km ('~km^2*')'
interacting.lab = "Day of year"
interested.length = the.interested.length
divider = 1000000
dist.too.far = 0.05

#min.expl = exp(min(unsex.trap[,names(unsex.trap)==intrest.covar ]))-1
#max.expl = exp(max(unsex.trap[,names(unsex.trap)==intrest.covar ]))-1
min.expl = min(unsex.trap[,names(unsex.trap)==intrest.covar ])
max.expl = max(unsex.trap[,names(unsex.trap)==intrest.covar ])
#max.expl = 250 # 2 traps over 250m, 5 sessions

min.inter = min(unsex.trap[,names(unsex.trap)==itneracting ])
max.inter = max(unsex.trap[,names(unsex.trap)==itneracting ])

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap, interested.covar = intrest.covar, 
                                     interested.length = interested.length, by.fac = F, interesting.factor.name = NULL)
# my.data1$exp.expl = seq(min.expl, max.expl, length = dim(my.data1)[1])
# my.data1[, which(names(my.data1)== intrest.covar)] = my.data1$exp.expl
my.data1[, which(names(my.data1) %in% c(intrest.covar ))] = seq(min.expl, max.expl, length = dim(my.data1)[1])

my.data = my.data1
my.data = my.data1[rep(1:(dim(my.data1)[1]), times = interested.length ),]
my.data[, which(names(my.data)== itneracting)] = rep (seq(min.inter, max.inter, length = interested.length), each = interested.length)


pred.bam.st.habseason.hab.density.2.0 =  predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")

my.data$fit = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit))
my.data$upper.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit + 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))
my.data$lower.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit - 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))

my.data$expl = my.data[,which(names(my.data)==intrest.covar)]
my.data$itneracting = my.data[,which(names(my.data)==itneracting)]

rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])

ylims = c(0 , round(max(my.data$upper.95), digits = 3)+0.001)

plotr = my.data[, which(names(my.data) %in% c(intrest.covar, itneracting ))]
plotr$fit = my.data$fit
plotr$upper.95 = my.data$upper.95
plotr$lower.95 = my.data$lower.95

plotr$expl = as.numeric(my.data[,which(names(my.data) %in% c(intrest.covar))]/divider)
plotr$intr = as.numeric(my.data[,which(names(my.data) == (itneracting))])

og.dat = unsex.trap[unsex.trap[,which(names(unsex.trap) %in% c( intrest.covar))]<max.expl,]
og.dat$intr = og.dat[,which(names(og.dat) %in% c( itneracting))]
og.dat$expl = og.dat[,which(names(og.dat) %in% c( intrest.covar))]/divider
rug = unique(og.dat[,which(names(og.dat) == intrest.covar)])/divider

caps= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c(intrest.covar, itneracting))]
caps$itneracting= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c( itneracting))]
caps$intrest.covar= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c( intrest.covar))]/divider

vars.too.far<- exclude.too.far(plotr$expl , plotr$intr,
                               og.dat$expl , og.dat$intr,
                               dist.too.far)
cap.point.size = 0.05

# # # #
#col.scale
if(abs(range(plotr$fit[!vars.too.far])[2])> abs(range(plotr$fit[!vars.too.far])[1])){
    col.scale.end = 1
    col.scale.start = 0.5 - (abs(range(plotr$fit[!vars.too.far])[1])/abs(range(plotr$fit[!vars.too.far])[2]))
} else{
    col.scale.start = 0
    col.scale.end = 0.5 + (abs(range(plotr$fit[!vars.too.far])[2])/abs(range(plotr$fit[!vars.too.far])[1]))
}
# # # #
ggplot() + 
    geom_tile(aes(plotr$expl[!vars.too.far], plotr$intr[!vars.too.far], fill = plotr$fit[!vars.too.far])) +
    # geom_rug(aes(x = rug), alpha = 1/2,sides = "b",size = 0.2)+
    labs(y = interacting.lab, x = expl.lab , z = "Target Sex Capture Probability") +
    geom_contour(aes(plotr$expl[!vars.too.far], plotr$intr[!vars.too.far], z = plotr$fit[!vars.too.far]), color = "grey")+
    #  scale_fill_viridis(begin = col.scale.start, end = col.scale.end)+
    scale_fill_viridis(name = "Trapping rate") +
    scale_y_continuous(breaks= c(1/12, 3/12, 5/12, 7/12, 9/12,11/12)*365, labels= c("Feb", "Apr", "June", "Aug","Oct", "Dec"))+
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )
# date of peak at low bog
my.data$itneracting[my.data$expl==0][my.data$fit[my.data$expl==0]==max( my.data$fit[my.data$expl==0]) ] # 1st sept = day 244
#lowest
my.data[my.data$fit == min(my.data$fit ),] # 1st march and 3.499km2

# predicted tensor effect of Water edge ----
intrest.covar = "water.edge.length"
itneracting = "julian"
expl.lab =   "Water edge within 1.5 km (km)"
interacting.lab = "Day of year"
interested.length = the.interested.length
divider = 10000
dist.too.far = 0.05
max.expl =  max(unsex.trap[,names(unsex.trap)==intrest.covar ]) # 2 traps over 250m, 5 sessions

#min.expl = exp(min(unsex.trap[,names(unsex.trap)==intrest.covar ]))-1
#max.expl = exp(max(unsex.trap[,names(unsex.trap)==intrest.covar ]))-1
min.expl = min(unsex.trap[,names(unsex.trap)==intrest.covar ])
#max.expl = max(unsex.trap[,names(unsex.trap)==intrest.covar ])

min.inter = min(unsex.trap[,names(unsex.trap)==itneracting ])
max.inter = max(unsex.trap[,names(unsex.trap)==itneracting ])


my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap, interested.covar = intrest.covar, 
                                     interested.length = interested.length, by.fac = F, interesting.factor.name = NULL)
# my.data1$exp.expl = seq(min.expl, max.expl, length = dim(my.data1)[1])
# my.data1[, which(names(my.data1)== intrest.covar)] = my.data1$exp.expl
my.data1[, which(names(my.data1) %in% c(intrest.covar ))] = seq(min.expl, max.expl, length = dim(my.data1)[1])


my.data = my.data1
my.data = my.data1[rep(1:(dim(my.data1)[1]), times = interested.length ),]
my.data[, which(names(my.data)== itneracting)] = rep (seq(min.inter, max.inter, length = interested.length), each = interested.length)

pred.bam.st.habseason.hab.density.2.0 =  predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")

my.data$fit = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit))
my.data$upper.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit + 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))
my.data$lower.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit - 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))
my.data$expl = my.data[,which(names(my.data)==intrest.covar)]
my.data$itneracting = my.data[,which(names(my.data)==itneracting)]

rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])

ylims = c(0 , round(max(my.data$upper.95), digits = 3)+0.001)

plotr = my.data[, which(names(my.data) %in% c(intrest.covar, itneracting ))]
plotr$fit = my.data$fit
plotr$upper.95 = my.data$upper.95
plotr$lower.95 = my.data$lower.95

plotr$expl = as.numeric(my.data[,which(names(my.data) %in% c(intrest.covar))]/divider)
plotr$intr = as.numeric(my.data[,which(names(my.data) == (itneracting))])

og.dat = unsex.trap[unsex.trap[,which(names(unsex.trap) %in% c( intrest.covar))]<max.expl,]
og.dat$intr = og.dat[,which(names(og.dat) %in% c( itneracting))]
og.dat$expl = og.dat[,which(names(og.dat) %in% c( intrest.covar))]
rug = unique(og.dat[,which(names(og.dat) == intrest.covar)])/divider

caps= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c(intrest.covar, itneracting))]
caps$itneracting= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c( itneracting))]
caps$intrest.covar= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c( intrest.covar))]/divider

vars.too.far<- exclude.too.far(plotr$expl , plotr$intr,
                               og.dat$expl/ divider, og.dat$intr,
                               dist.too.far)
cap.point.size = 0.05

# # # #
#col.scale
if(abs(range(plotr$fit[!vars.too.far])[2])> abs(range(plotr$fit[!vars.too.far])[1])){
    col.scale.end = 1
    col.scale.start = 0.5 - (abs(range(plotr$fit[!vars.too.far])[1])/abs(range(plotr$fit[!vars.too.far])[2]))
} else{
    col.scale.start = 0
    col.scale.end = 0.5 + (abs(range(plotr$fit[!vars.too.far])[2])/abs(range(plotr$fit[!vars.too.far])[1]))
}
# # # #
ggplot() + 
    geom_tile(aes(plotr$expl[!vars.too.far], plotr$intr[!vars.too.far], fill = plotr$fit[!vars.too.far])) +
    # geom_rug(aes(x = rug), alpha = 1/2,sides = "b",size = 0.2)+
    labs(y = interacting.lab, x = expl.lab , z = "Target Sex Capture Probability") +
    geom_contour(aes(plotr$expl[!vars.too.far], plotr$intr[!vars.too.far], z = plotr$fit[!vars.too.far]), color = "grey")+
    #  scale_fill_viridis(begin = col.scale.start, end = col.scale.end)+
    scale_fill_viridis(name = "Trapping rate") +
    scale_y_continuous(breaks= c(1/12, 3/12, 5/12, 7/12, 9/12,11/12)*365, labels= c("Feb", "Apr", "June", "Aug","Oct", "Dec"))+
    scale_x_continuous(breaks= seq(0,70, by = 10)/10, labels= seq(0,70, by = 10))+
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )


##########################################################################################
#
##########################################################################################
# TENSOR OF DENSITY -----
##########################################################################################
# altitude ----
intrest.covar = "ALT3"
expl.lab =   "Altitude (m)"
intrest.term = "ti(log.dens,ALT3)"
itneracting = "log.dens"
breaks.by = 100
break.divider = 1
interacting.lab = "Regional abundance (N)"
dist.too.far = 0.06
divider = 1
dens.to.axis = c(100,250,500,1000, 2000)

interested.length = the.interested.length
min.inter = min(unsex.trap[,names(unsex.trap)==itneracting ])
max.inter = max(unsex.trap[,names(unsex.trap)==itneracting ])
min.expl = min(unsex.trap[,names(unsex.trap)==intrest.covar ])
max.expl = max(unsex.trap[,names(unsex.trap)==intrest.covar ])

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap, interested.covar = intrest.covar, 
                                     interested.length = interested.length, by.fac = F, interesting.factor.name = NULL)
# my.data1$exp.expl = seq(min.expl, max.expl, length = dim(my.data1)[1])
# my.data1[, which(names(my.data1)== intrest.covar)] = my.data1$exp.expl
my.data1[, which(names(my.data1) %in% c(intrest.covar ))] = seq(min.expl, max.expl, length = dim(my.data1)[1])

my.data = my.data1
my.data = my.data1[rep(1:(dim(my.data1)[1]), times = interested.length ),]
my.data[, which(names(my.data)== itneracting)] = rep (seq(min.inter, max.inter, length = interested.length), each = interested.length)

#my.data[,which(names(my.data)==intrest.covar)]  = log(1+ seq(0, exp(max(unsex.trap[,which(names(unsex.trap)==intrest.covar)])), length = dim(my.data)[1] ))
pred.bam.st.habseason.hab.density.2.0.dum = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "terms")

pred.bam.st.habseason.hab.density.2.0 = pred.bam.st.habseason.hab.density.2.0.dum

my.data$fit = as.numeric( pred.bam.st.habseason.hab.density.2.0.dum$fit[,which(attributes(pred.bam.st.habseason.hab.density.2.0.dum$fit)$dimnames[[2]] == intrest.term)])
my.data$expl = my.data[,which(names(my.data)==intrest.covar)]  
x.max = 300#max(my.data$expl)

rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])
ylims = c(0 , round(max(my.data$upper.95), digits = 3)+0.001)

plotr = my.data[, which(names(my.data) %in% c(intrest.covar, itneracting ))]
plotr$fit = my.data$fit

plotr$expl = as.numeric(my.data[,which(names(my.data) %in% c(intrest.covar))]/divider)
plotr$intr = as.numeric(my.data[,which(names(my.data) == (itneracting))])

og.dat = unsex.trap[unsex.trap[,which(names(unsex.trap) %in% c( intrest.covar))]<max.expl,]
og.dat$intr = og.dat[,which(names(og.dat) %in% c( itneracting))]
og.dat$expl = og.dat[,which(names(og.dat) %in% c( intrest.covar))]
rug = unique(og.dat[,which(names(og.dat) == intrest.covar)])/divider

caps= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c(intrest.covar, itneracting))]
caps$itneracting= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c( itneracting))]
caps$intrest.covar= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c( intrest.covar))]/divider

vars.too.far<- exclude.too.far(plotr$expl , plotr$intr,
                               og.dat$expl/divider , og.dat$intr,
                               dist.too.far)
cap.point.size = 0.01

# # # #
#col.scale
if(abs(range(plotr$fit[!vars.too.far])[2])> abs(range(plotr$fit[!vars.too.far])[1])){
    col.scale.end = 1
    col.scale.start = 0.5 - (abs(range(plotr$fit[!vars.too.far])[1])/abs(range(plotr$fit[!vars.too.far])[2]))
} else{
    col.scale.start = 0
    col.scale.end = 0.5 + (abs(range(plotr$fit[!vars.too.far])[2])/abs(range(plotr$fit[!vars.too.far])[1]))
}
# # # #
ggplot() + 
    geom_tile(aes(plotr$expl[!vars.too.far], plotr$intr[!vars.too.far], fill = plotr$fit[!vars.too.far])) +
    #geom_point(aes(caps$intrest.covar, caps$itneracting), size = cap.point.size, color = "white")+
    # geom_rug(aes(x = rug), alpha = 1/2,sides = "b",size = 0.2)+
    labs(y = interacting.lab, x = expl.lab , z = "Target Sex Capture Probability") +
    geom_contour(aes(plotr$expl[!vars.too.far], plotr$intr[!vars.too.far], z = plotr$fit[!vars.too.far]), color = "grey")+
    #  scale_fill_viridis(begin = col.scale.start, end = col.scale.end)+
    scale_fill_viridis(name = "Effect size") +
    scale_y_continuous(breaks = log(dens.to.axis + 1) , labels = dens.to.axis )+
    scale_x_continuous(breaks = seq(0,x.max, by = breaks.by ), labels= seq(0,x.max, by = breaks.by )/break.divider ,  limits=c(0,max(my.data$expl))) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )

##########################################################################################
# DENSITY TENSOR PREDICTION
##########################################################################################
# tensor prediction of altitude ----
intrest.covar = "ALT3"
itneracting = "log.dens"
expl.lab =   "Altitude (m)"
interacting.lab = "Regional abundance (N)"
interested.length = the.interested.length
divider = 1
dens.to.axis = c(100,250,500,1000, 2000)


min.expl = min(unsex.trap[,names(unsex.trap)==intrest.covar ])
max.expl = max(unsex.trap[,names(unsex.trap)==intrest.covar ])
max.expl = 300 

#min.interact = min(unsex.trap$block.dens)
#max.interact = max(unsex.trap$block.dens)
min.interact = min(unsex.trap$log.dens)
max.interact = max(unsex.trap$log.dens)


my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap, interested.covar = intrest.covar, 
                                     interested.length = interested.length, by.fac = F, interesting.factor.name = NULL)
# my.data1$exp.expl = seq(min.expl, max.expl, length = dim(my.data1)[1])
# my.data1[, which(names(my.data1)== intrest.covar)] = my.data1$exp.expl
my.data1[, which(names(my.data1) %in% c(intrest.covar ))] = seq(min.expl, max.expl, length = dim(my.data1)[1])


my.data = my.data1
my.data = my.data1[rep(1:(dim(my.data1)[1]), times = interested.length ),]
#my.data[, which(names(my.data)== itneracting)] = log(rep (seq(min.interact, max.interact, length = interested.length), each = interested.length)+1)
my.data[, which(names(my.data)== itneracting)] = rep (seq(min.interact, max.interact, length = interested.length), each = interested.length)


pred.bam.st.habseason.hab.density.2.0 = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")
# extra stuff

my.data$fit = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit))
my.data$upper.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit + 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))
my.data$lower.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit - 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))
# end of hack
my.data$expl = my.data[,which(names(my.data)==intrest.covar)]
#my.data$itneracting =exp( my.data[,which(names(my.data)==itneracting)])-1
my.data$itneracting =my.data[,which(names(my.data)==itneracting)]


rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])

ylims = c(0 , round(max(my.data$upper.95), digits = 3)+0.001)

plotr = my.data[, which(names(my.data) %in% c(intrest.covar, itneracting ))]
plotr$fit = my.data$fit
plotr$upper.95 = my.data$upper.95
plotr$lower.95 = my.data$lower.95

plotr$expl =  my.data$expl
plotr$intr = my.data$itneracting

og.dat = unsex.trap[unsex.trap[,which(names(unsex.trap) %in% c( intrest.covar))]<max.expl,]
og.dat$intr = og.dat[,which(names(og.dat) %in% c( itneracting))]
#og.dat$intr = exp(og.dat[,which(names(og.dat) %in% c( itneracting))])-1

og.dat$expl = og.dat[,which(names(og.dat) %in% c( intrest.covar))]
rug = unique(og.dat[,which(names(og.dat) == intrest.covar)])/divider

caps= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c(intrest.covar, itneracting))]
#caps$itneracting= exp(og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c( itneracting))])-1
caps$itneracting= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c( itneracting))]
caps$intrest.covar= og.dat[og.dat$sex.cap==1,which(names(og.dat) %in% c( intrest.covar))]/divider

dist.too.far = 0.1
vars.too.far<- exclude.too.far(plotr$expl , plotr$intr,
                               og.dat$expl , og.dat$intr,
                               dist.too.far)
cap.point.size = 0.05

# # # #
#col.scale
if(abs(range(plotr$fit[!vars.too.far])[2])> abs(range(plotr$fit[!vars.too.far])[1])){
    col.scale.end = 1
    col.scale.start = 0.5 - (abs(range(plotr$fit[!vars.too.far])[1])/abs(range(plotr$fit[!vars.too.far])[2]))
} else{
    col.scale.start = 0
    col.scale.end = 0.5 + (abs(range(plotr$fit[!vars.too.far])[2])/abs(range(plotr$fit[!vars.too.far])[1]))
}
# # # #
ggplot() + 
    geom_tile(aes(plotr$expl[!vars.too.far], plotr$intr[!vars.too.far], fill = plotr$fit[!vars.too.far])) +
 #   geom_point(aes(caps$intrest.covar, caps$itneracting), size = 0.1, color = "white")+
    # geom_rug(aes(x = rug), alpha = 1/2,sides = "b",size = 0.2)+
  labs(y = interacting.lab, x = expl.lab , z = "Target Sex Capture Probability") +
  geom_contour(aes(plotr$expl[!vars.too.far], plotr$intr[!vars.too.far], z = plotr$fit[!vars.too.far]), color = "white")+
    #  scale_fill_viridis(begin = col.scale.start, end = col.scale.end)+
  scale_fill_viridis(name = "Trapping rate") +
  scale_y_continuous(breaks=log(dens.to.axis+1), labels=dens.to.axis ,  limits=log(range(dens.to.axis+1)))+
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )

# lowest abundance highest quality ----
my.data$expl[my.data$itneracting == min(my.data$itneracting )][my.data$fit[my.data$itneracting == min(my.data$itneracting )] == max(my.data$fit[my.data$itneracting == min(my.data$itneracting )])]
##########################################################################################
