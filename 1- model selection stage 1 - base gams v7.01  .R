##################################################################################################################################################################
# base bam analaysis, altered from "base bams 7.5.18.R" 
# changes to reflect superviror meeting on 14.5.18
###     loaded unsex.trap.mk6 instead of mk5
###     include scent effects
###     include trap competition
###     include local population depletion
###     changed all bam() to bam()



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

##################################################################################################################################################################
##################################################################################################################################################################
# run functions code -- needs to be done ----
source("functions.code.R")

##################################################################################################################################################################
# load .csvs ----
par(mfrow = c(1,1))
def.par = par()
unsex.trap = read.csv("unsex.trap.v7.01.csv") 


##################################################################################################################################################################


##################################################################################################################################################################
# run base models ----

time0 = Sys.time()
bam.logdens.base.1  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                      s(log.dens) + 
                      sub.ses.num +
                      s(sub.sesh.night.num, k = 4) +

                      s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                  data = unsex.trap , family = binomial, method = "ML")
time.bam.logdens.base.1 = Sys.time() - time0
save(bam.logdens.base.1, file = paste(sep = "", model.outputs.path, "bam.logdens.base.1"))
rm(bam.logdens.base.1)
gc()

time0 = Sys.time()
bam.logdens.base.2  = bam(unsex.capt ~ s(julian, bs = "cc") +
                      s(log.dens) + 
                      sub.ses.num +
                      s(sub.sesh.night.num, k = 4) +

                    mink.scent +
                          
                      s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                  data = unsex.trap , family = binomial, method = "ML")
time.bam.logdens.base.2 = Sys.time() - time0
save(bam.logdens.base.2, file = paste(sep = "", model.outputs.path, "bam.logdens.base.2"))
rm(bam.logdens.base.2)
gc()


###################################################################################################################################################################
 # load models and compare aics ----

load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.logdens.base.1")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.logdens.base.2")

AIC(bam.logdens.base.1, bam.logdens.base.2)

summary(bam.logdens.base.26)

###################################################################################################################################################################
