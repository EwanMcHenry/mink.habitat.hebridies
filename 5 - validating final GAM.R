# 
# validation
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
library(data.table)
library(MuMIn)
library(TTR)

library(smooth)
library(Mcomp)

library(corrplot)


##################################################################################################################################################################

##################################################################################################################################################################
# load .csvs ----
par(mfrow = c(1,1))
def.par = par()
unsex.trap = read.csv("unsex.trap.v7.01.csv") 
#unsex.trap$sesid = as.factor(unsex.trap$sesid)
##################################################################################################################################################################
# run functions code -- needs to be done ----
source("functions.code.R")

# load model object
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.2.0")

# plot residuals by space

grid.size = 1000
east.sq = (floor(unsex.trap$east/grid.size) + 0.5) * grid.size
north.sq = (floor(unsex.trap$north/grid.size) + 0.5)* grid.size

kim.sq = paste(east.sq, north.sq)

resid.plot = data.frame(resids = residuals(bam.st.habseason.hab.density.2.0, type= "pearson"),
                        kim.sq = paste(east.sq, north.sq))

km.mean.resids = aggregate(resid.plot$resids, by = list(east.sq, north.sq), mean )

ggplot(km.mean.resids)+
    geom_point(aes(x = Group.1, y = Group.2, color = x)) +
    xlab("Eastings (m)") +
    ylab("Northings (m)") +
    scale_color_viridis(name = "mean residual")+
    coord_equal()+
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12, face = "bold")
    )

hist(km.mean.resids$x, xlim = c(-0.8,0.8), breaks = 200)
#######################################################################
# plot residuals by fitted values
resids = rstudent(bam.st.habseason.hab.density.2.0)

    residuals(bam.st.habseason.hab.density.2.0, type= "pearson")
rstandard(bam.st.habseason.hab.density.2.0)

my.data= unsex.trap
fit = pred.bam.st.habseason.hab.density.2.0 = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=F, type = "link")

plot(resids~ fit)


#######################################################################
# plot residuals by julian

num.data.points = 5000
obs.per.point = ceiling(length(unsex.trap$julian)/num.data.points)

resids = residuals(bam.st.habseason.hab.density.2.0, type= "pearson" )
#    - mean(residuals(bam.st.habseason.hab.density.2.0)))/ sd(residuals(bam.st.habseason.hab.density.2.0))

julian.validation = data.frame(julian = sort(unsex.trap$julian),
                               resids = resids[order(unsex.trap$julian)])
julian.validation$bin.num = rep(1:num.data.points, each = obs.per.point, length =  length(unsex.trap$julian))

julian.validation$mean.julian = NA
julian.validation$mean.resid = NA
for (i in 1: num.data.points){
    julian.validation$mean.julian[julian.validation$bin.num == i] = mean(julian.validation$julian[julian.validation$bin.num == i])
    julian.validation$mean.resid[julian.validation$bin.num == i] = mean(julian.validation$resids[julian.validation$bin.num == i])
}

plot.julian.resids = julian.validation[!duplicated(julian.validation$bin.num),]

ggplot(plot.julian.resids)+
    geom_point(aes(x = mean.julian, y = mean.resid)) +
    xlab("julian date") +
    ylab("residual") +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    )

