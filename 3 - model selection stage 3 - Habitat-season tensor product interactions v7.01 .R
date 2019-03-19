# examination of habitat-seasonality interactions

#################################################################################################################################################################
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
library(GGally)

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
# lvl 1 models ----

time0 = Sys.time()
bam.s.t.hab.season.1.0 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                  s(log.dens) + 
                                  sub.ses.num +
                                  s(sub.sesh.night.num, k = 4) +
                                  
                                  mink.scent +
                                  
                                  s(salat.sum.means) +
                                  coastal +
                                  s(combo.heather_all.mk2) + 
                                  s(ALT3) +
                                  s(builtup.mk2) + 
                                  s(rough.grass.mk2) + 
                                  s(water.edge.length) +
                                  s(combo.beach_saltmarsh.mk2) +
                                  s(pasture.mk2) + 
                                  s(bog.mk2) +
                                  
                                  s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                              data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.1.0 = Sys.time() - time0
save(bam.s.t.hab.season.1.0, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.1.0"))
rm(bam.s.t.hab.season.1.0)
gc()


time0 = Sys.time()
bam.s.t.hab.season.1.01 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   ti(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   s(ALT3) +
                                   s(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   s(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, salat.sum.means, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.1.01 = Sys.time() - time0
save(bam.s.t.hab.season.1.01, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.1.01"))
rm(bam.s.t.hab.season.1.01)
gc()

time0 = Sys.time()
bam.s.t.hab.season.1.02 =  bam(unsex.capt ~ ti(julian, by = coastal, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   s(ALT3) +
                                   s(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   s(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.1.02 = Sys.time() - time0
save(bam.s.t.hab.season.1.02, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.1.02"))
rm(bam.s.t.hab.season.1.02)
gc()

time0 = Sys.time()
bam.s.t.hab.season.1.03 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   ti(combo.heather_all.mk2) + 
                                   s(ALT3) +
                                   s(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   s(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, combo.heather_all.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.1.03 = Sys.time() - time0
save(bam.s.t.hab.season.1.03, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.1.03"))
rm(bam.s.t.hab.season.1.03)
gc()

time0 = Sys.time()
bam.s.t.hab.season.1.04 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   s(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.1.04 = Sys.time() - time0
save(bam.s.t.hab.season.1.04, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.1.04"))
rm(bam.s.t.hab.season.1.04)
gc()

time0 = Sys.time()
bam.s.t.hab.season.1.05 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   s(ALT3) +
                                   ti(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   s(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, builtup.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.1.05 = Sys.time() - time0
save(bam.s.t.hab.season.1.05, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.1.05"))
rm(bam.s.t.hab.season.1.05)
gc()

time0 = Sys.time()
bam.s.t.hab.season.1.06 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   s(ALT3) +
                                   s(builtup.mk2) + 
                                   ti(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   s(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, rough.grass.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.1.06 = Sys.time() - time0
save(bam.s.t.hab.season.1.06, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.1.06"))
rm(bam.s.t.hab.season.1.06)
gc()

time0 = Sys.time()
bam.s.t.hab.season.1.07 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   s(ALT3) +
                                   s(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   ti(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   s(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, water.edge.length, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.1.07 = Sys.time() - time0
save(bam.s.t.hab.season.1.07, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.1.07"))
rm(bam.s.t.hab.season.1.07)
gc()

time0 = Sys.time()
bam.s.t.hab.season.1.08 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   s(ALT3) +
                                   s(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   ti(combo.beach_saltmarsh.mk2) +
                                   s(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, combo.beach_saltmarsh.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.1.08 = Sys.time() - time0
save(bam.s.t.hab.season.1.08, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.1.08"))
rm(bam.s.t.hab.season.1.08)
gc()

time0 = Sys.time()
bam.s.t.hab.season.1.09 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   s(ALT3) +
                                   s(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.1.09 = Sys.time() - time0
save(bam.s.t.hab.season.1.09, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.1.09"))
rm(bam.s.t.hab.season.1.09)
gc()

time0 = Sys.time()
bam.s.t.hab.season.1.10 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   s(ALT3) +
                                   s(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   s(pasture.mk2) + 
                                   ti(bog.mk2) +
                                   
                                   ti(julian, bog.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.1.10 = Sys.time() - time0
save(bam.s.t.hab.season.1.10, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.1.10"))
rm(bam.s.t.hab.season.1.10)
gc()

##################################################################################################################################################################
# lvl 1 aic comparisions ----

load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.1.0")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.1.01")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.1.02")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.1.03")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.1.04")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.1.05")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.1.06")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.1.07")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.1.08")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.1.09")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.1.10")

AIC(bam.s.t.hab.season.1.0, bam.s.t.hab.season.1.01, 
    bam.s.t.hab.season.1.02, bam.s.t.hab.season.1.03, 
    bam.s.t.hab.season.1.04, 
    bam.s.t.hab.season.1.05, 
    bam.s.t.hab.season.1.06, 
    bam.s.t.hab.season.1.07, 
    bam.s.t.hab.season.1.08, 
    bam.s.t.hab.season.1.09, 
    bam.s.t.hab.season.1.10
)

# df      AIC
# bam.s.t.hab.season.1.0  121.9077 15001.90
# bam.s.t.hab.season.1.01 126.5229 14976.00
# bam.s.t.hab.season.1.02 121.2483 14971.37
# bam.s.t.hab.season.1.03 125.2391 15005.64
# bam.s.t.hab.season.1.04 124.5100 14952.16 * * ALT3
# bam.s.t.hab.season.1.05 125.8361 15000.06
# bam.s.t.hab.season.1.06 133.0078 14984.59
# bam.s.t.hab.season.1.07 129.0782 15001.67
# bam.s.t.hab.season.1.08 125.3632 14993.44
# bam.s.t.hab.season.1.09 133.2610 14967.36
# bam.s.t.hab.season.1.10 124.4744 14978.29
##################################################################################################################################################################
# lvl 2 models ----

time0 = Sys.time()
bam.s.t.hab.season.2.0 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                  s(log.dens) + 
                                  sub.ses.num +
                                  s(sub.sesh.night.num, k = 4) +
                                  
                                  mink.scent +
                                  
                                  s(salat.sum.means) +
                                  coastal +
                                  s(combo.heather_all.mk2) + 
                                  ti(ALT3) +
                                  s(builtup.mk2) + 
                                  s(rough.grass.mk2) + 
                                  s(water.edge.length) +
                                  s(combo.beach_saltmarsh.mk2) +
                                  s(pasture.mk2) + 
                                  s(bog.mk2) +
                                  
                                  ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                  
                                  s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                              data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.2.0 = Sys.time() - time0
save(bam.s.t.hab.season.2.0, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.2.0"))
rm(bam.s.t.hab.season.2.0)
gc()

time0 = Sys.time()
bam.s.t.hab.season.2.01 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   ti(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   s(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, salat.sum.means, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.2.01 = Sys.time() - time0
save(bam.s.t.hab.season.2.01, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.2.01"))
rm(bam.s.t.hab.season.2.01)
gc()

bam.s.t.hab.season.2.02 =  bam(unsex.capt ~ ti(julian, by = coastal, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   s(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +

                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.2.02 = Sys.time() - time0
save(bam.s.t.hab.season.2.02, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.2.02"))
rm(bam.s.t.hab.season.2.02)
gc()

bam.s.t.hab.season.2.03 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   ti(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   s(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, combo.heather_all.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.2.03 = Sys.time() - time0
save(bam.s.t.hab.season.2.03, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.2.03"))
rm(bam.s.t.hab.season.2.03)
gc()

bam.s.t.hab.season.2.04 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   ti(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   s(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, builtup.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.2.04 = Sys.time() - time0
save(bam.s.t.hab.season.2.04, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.2.04"))
rm(bam.s.t.hab.season.2.04)
gc()

bam.s.t.hab.season.2.05 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   ti(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   s(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, rough.grass.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.2.05 = Sys.time() - time0
save(bam.s.t.hab.season.2.05, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.2.05"))
rm(bam.s.t.hab.season.2.05)
gc()

bam.s.t.hab.season.2.06 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   ti(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   s(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, water.edge.length, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.2.06 = Sys.time() - time0
save(bam.s.t.hab.season.2.06, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.2.06"))
rm(bam.s.t.hab.season.2.06)
gc()

bam.s.t.hab.season.2.07 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   ti(combo.beach_saltmarsh.mk2) +
                                   s(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, combo.beach_saltmarsh.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.2.07 = Sys.time() - time0
save(bam.s.t.hab.season.2.07, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.2.07"))
rm(bam.s.t.hab.season.2.07)
gc()

bam.s.t.hab.season.2.08 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.2.08 = Sys.time() - time0
save(bam.s.t.hab.season.2.08, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.2.08"))
rm(bam.s.t.hab.season.2.08)
gc()

bam.s.t.hab.season.2.09 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   s(pasture.mk2) + 
                                   ti(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, bog.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.2.09 = Sys.time() - time0
save(bam.s.t.hab.season.2.09, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.2.09"))
rm(bam.s.t.hab.season.2.09)
gc()

##################################################################################################################################################################
# lvl 2 aic comparisions ----

load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.2.0")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.2.01")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.2.02")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.2.03")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.2.04")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.2.05")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.2.06")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.2.07")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.2.08")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.2.09")

AIC(bam.s.t.hab.season.2.0, bam.s.t.hab.season.2.01, 
    bam.s.t.hab.season.2.02, bam.s.t.hab.season.2.03, 
    bam.s.t.hab.season.2.04, bam.s.t.hab.season.2.05, 
    bam.s.t.hab.season.2.06, bam.s.t.hab.season.2.07 , 
    bam.s.t.hab.season.2.08, 
    bam.s.t.hab.season.2.09
)

# df      AIC
# bam.s.t.hab.season.2.0  124.5100 14952.16
# bam.s.t.hab.season.2.01 124.1767 14952.37
# bam.s.t.hab.season.2.02 124.9890 14949.87
# bam.s.t.hab.season.2.03 125.3691 14953.92
# bam.s.t.hab.season.2.04 120.3874 14949.30
# bam.s.t.hab.season.2.05 130.3818 14950.80
# bam.s.t.hab.season.2.06 128.0656 14951.30
# bam.s.t.hab.season.2.07 128.2824 14954.87
# bam.s.t.hab.season.2.08 133.1309 14936.26 ** pasture
# bam.s.t.hab.season.2.09 124.7576 14949.24

##################################################################################################################################################################
# lvl 3.0 models ----


time0 = Sys.time()
bam.s.t.hab.season.3.0 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                  s(log.dens) + 
                                  sub.ses.num +
                                  s(sub.sesh.night.num, k = 4) +
                                  
                                  mink.scent +
                                  
                                  s(salat.sum.means) +
                                  coastal +
                                  s(combo.heather_all.mk2) + 
                                  ti(ALT3) +
                                  s(builtup.mk2) + 
                                  s(rough.grass.mk2) + 
                                  s(water.edge.length) +
                                  s(combo.beach_saltmarsh.mk2) +
                                  ti(pasture.mk2) + 
                                  s(bog.mk2) +
                                  
                                  ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                  ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                  
                                  s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                              data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.3.0 = Sys.time() - time0
save(bam.s.t.hab.season.3.0, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.3.0"))
rm(bam.s.t.hab.season.3.0)
gc()

time0 = Sys.time()
bam.s.t.hab.season.3.01 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   ti(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, salat.sum.means, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.3.01 = Sys.time() - time0
save(bam.s.t.hab.season.3.01, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.3.01"))
rm(bam.s.t.hab.season.3.01)
gc()

time0 = Sys.time()
bam.s.t.hab.season.3.02 =  bam(unsex.capt ~ ti(julian, by = coastal,  bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +

                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.3.02 = Sys.time() - time0
save(bam.s.t.hab.season.3.02, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.3.02"))
rm(bam.s.t.hab.season.3.02)
gc()

time0 = Sys.time()
bam.s.t.hab.season.3.03 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   ti(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, combo.heather_all.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.3.03 = Sys.time() - time0
save(bam.s.t.hab.season.3.03, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.3.03"))
rm(bam.s.t.hab.season.3.03)
gc()

time0 = Sys.time()
bam.s.t.hab.season.3.04 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   ti(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, builtup.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.3.04 = Sys.time() - time0
save(bam.s.t.hab.season.3.04, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.3.04"))
rm(bam.s.t.hab.season.3.04)
gc()

time0 = Sys.time()
bam.s.t.hab.season.3.05 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   ti(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, rough.grass.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.3.05 = Sys.time() - time0
save(bam.s.t.hab.season.3.05, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.3.05"))
rm(bam.s.t.hab.season.3.05)
gc()

time0 = Sys.time()
bam.s.t.hab.season.3.06 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   ti(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, water.edge.length, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.3.06 = Sys.time() - time0
save(bam.s.t.hab.season.3.06, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.3.06"))
rm(bam.s.t.hab.season.3.06)
gc()

time0 = Sys.time()
bam.s.t.hab.season.3.07 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   ti(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, combo.beach_saltmarsh.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.3.07 = Sys.time() - time0
save(bam.s.t.hab.season.3.07, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.3.07"))
rm(bam.s.t.hab.season.3.07)
gc()

time0 = Sys.time()
bam.s.t.hab.season.3.08 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   s(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   ti(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, bog.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.3.08 = Sys.time() - time0
save(bam.s.t.hab.season.3.08, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.3.08"))
rm(bam.s.t.hab.season.3.08)
gc()

##################################################################################################################################################################
# lvl 3 aic comparisions ----

load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.3.0")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.3.01")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.3.02")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.3.03")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.3.04")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.3.05")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.3.06")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.3.07")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.3.08")

AIC(bam.s.t.hab.season.3.0, bam.s.t.hab.season.3.01, 
    bam.s.t.hab.season.3.02, bam.s.t.hab.season.3.03, bam.s.t.hab.season.3.04, 
    bam.s.t.hab.season.3.05, bam.s.t.hab.season.3.06, bam.s.t.hab.season.3.07, bam.s.t.hab.season.3.08
)

# df      AIC
# bam.s.t.hab.season.3.0  133.1309 14936.26
# bam.s.t.hab.season.3.01 136.1264 14943.27
# bam.s.t.hab.season.3.02 131.3095 14931.85
# bam.s.t.hab.season.3.03 131.1783 14932.38
# bam.s.t.hab.season.3.04 147.3104 14971.31
# bam.s.t.hab.season.3.05 135.2021 14931.35 ** rough.grass.mk2
# bam.s.t.hab.season.3.06 135.5922 14934.88
# bam.s.t.hab.season.3.07 134.4256 14956.60
# bam.s.t.hab.season.3.08 137.6997 14941.63

##################################################################################################################################################################
# lvl 4.0 models ----

time0 = Sys.time()
bam.s.t.hab.season.4.0 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                  s(log.dens) + 
                                  sub.ses.num +
                                  s(sub.sesh.night.num, k = 4) +
                                  
                                  mink.scent +
                                  
                                  s(salat.sum.means) +
                                  coastal +
                                  s(combo.heather_all.mk2) + 
                                  ti(ALT3) +
                                  s(builtup.mk2) + 
                                  ti(rough.grass.mk2) + 
                                  s(water.edge.length) +
                                  s(combo.beach_saltmarsh.mk2) +
                                  ti(pasture.mk2) + 
                                  s(bog.mk2) +
                                  
                                  ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                  ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                  ti(julian, rough.grass.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                  
                                  s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                              data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.4.0 = Sys.time() - time0
save(bam.s.t.hab.season.4.0, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.4.0"))
rm(bam.s.t.hab.season.4.0)
gc()

time0 = Sys.time()
bam.s.t.hab.season.4.01 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   ti(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   ti(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, rough.grass.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, salat.sum.means, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.4.01 = Sys.time() - time0
save(bam.s.t.hab.season.4.01, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.4.01"))
rm(bam.s.t.hab.season.4.01)
gc()

time0 = Sys.time()
bam.s.t.hab.season.4.02 =  bam(unsex.capt ~ ti(julian, by = coastal, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   ti(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, rough.grass.mk2, bs = c("cc", "cr"), d = c(1,1)) +

                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.4.02 = Sys.time() - time0
save(bam.s.t.hab.season.4.02, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.4.02"))
rm(bam.s.t.hab.season.4.02)
gc()

time0 = Sys.time()
bam.s.t.hab.season.4.03 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   ti(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   ti(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, rough.grass.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, combo.heather_all.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.4.03 = Sys.time() - time0
save(bam.s.t.hab.season.4.03, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.4.03"))
rm(bam.s.t.hab.season.4.03)
gc()

time0 = Sys.time()
bam.s.t.hab.season.4.04 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   ti(builtup.mk2) + 
                                   ti(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, rough.grass.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, builtup.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.4.04 = Sys.time() - time0
save(bam.s.t.hab.season.4.04, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.4.04"))
rm(bam.s.t.hab.season.4.04)
gc()

time0 = Sys.time()
bam.s.t.hab.season.4.05 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   ti(rough.grass.mk2) + 
                                   ti(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, rough.grass.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, water.edge.length, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.4.05 = Sys.time() - time0
save(bam.s.t.hab.season.4.05, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.4.05"))
rm(bam.s.t.hab.season.4.05)
gc()

time0 = Sys.time()
bam.s.t.hab.season.4.06 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   ti(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   ti(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, rough.grass.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, combo.beach_saltmarsh.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.4.06 = Sys.time() - time0
save(bam.s.t.hab.season.4.06, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.4.06"))
rm(bam.s.t.hab.season.4.06)
gc()

time0 = Sys.time()
bam.s.t.hab.season.4.07 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   ti(rough.grass.mk2) + 
                                   s(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   ti(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, rough.grass.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, bog.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.4.07 = Sys.time() - time0
save(bam.s.t.hab.season.4.07, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.4.07"))
rm(bam.s.t.hab.season.4.07)
gc()


##################################################################################################################################################################
# lvl 4 aic comparisions ----

load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.4.0")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.4.01")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.4.02")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.4.03")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.4.04")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.4.05")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.4.06")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.4.07")

AIC(bam.s.t.hab.season.4.0, bam.s.t.hab.season.4.01, 
    bam.s.t.hab.season.4.02, 
    bam.s.t.hab.season.4.03,
    bam.s.t.hab.season.4.04, 
    bam.s.t.hab.season.4.05, bam.s.t.hab.season.4.06, 
    bam.s.t.hab.season.4.07)

# df      AIC
# bam.s.t.hab.season.4.0  135.2021 14931.35
# bam.s.t.hab.season.4.01 135.4641 14932.15
# bam.s.t.hab.season.4.02 134.5728 14928.42
# bam.s.t.hab.season.4.03 134.0626 14929.17
# bam.s.t.hab.season.4.04 132.1104 14932.43
# bam.s.t.hab.season.4.05 138.5051 14927.94 *** water edge
# bam.s.t.hab.season.4.06 138.3233 14934.30
# bam.s.t.hab.season.4.07 155.9146 14968.30


##################################################################################################################################################################
# lvl 5.0 models ----

time0 = Sys.time()
bam.s.t.hab.season.5.0 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                  s(log.dens) + 
                                  sub.ses.num +
                                  s(sub.sesh.night.num, k = 4) +
                                  
                                  mink.scent +
                                  
                                  s(salat.sum.means) +
                                  coastal +
                                  s(combo.heather_all.mk2) + 
                                  ti(ALT3) +
                                  s(builtup.mk2) + 
                                  ti(rough.grass.mk2) + 
                                  ti(water.edge.length) +
                                  s(combo.beach_saltmarsh.mk2) +
                                  ti(pasture.mk2) + 
                                  s(bog.mk2) +
                                  
                                  ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                  ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                  ti(julian, rough.grass.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                  ti(julian, water.edge.length, bs = c("cc", "cr"), d = c(1,1)) +
                                  
                                  s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                              data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.5.0 = Sys.time() - time0
save(bam.s.t.hab.season.5.0, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.5.0"))
rm(bam.s.t.hab.season.5.0)
gc()

time0 = Sys.time()
bam.s.t.hab.season.5.01 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   ti(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   ti(rough.grass.mk2) + 
                                   ti(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, rough.grass.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, water.edge.length, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, salat.sum.means, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.5.01 = Sys.time() - time0
save(bam.s.t.hab.season.5.01, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.5.01"))
rm(bam.s.t.hab.season.5.01)
gc()

time0 = Sys.time()
bam.s.t.hab.season.5.02 =  bam(unsex.capt ~ ti(julian, by = coastal, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   ti(rough.grass.mk2) + 
                                   ti(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, rough.grass.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, water.edge.length, bs = c("cc", "cr"), d = c(1,1)) +

                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.5.02 = Sys.time() - time0
save(bam.s.t.hab.season.5.02, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.5.02"))
rm(bam.s.t.hab.season.5.02)
gc()

time0 = Sys.time()
bam.s.t.hab.season.5.03 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   ti(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   ti(rough.grass.mk2) + 
                                   ti(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, rough.grass.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, water.edge.length, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, combo.heather_all.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.5.03 = Sys.time() - time0
save(bam.s.t.hab.season.5.03, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.5.03"))
rm(bam.s.t.hab.season.5.03)
gc()

time0 = Sys.time()
bam.s.t.hab.season.5.04 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   ti(builtup.mk2) + 
                                   ti(rough.grass.mk2) + 
                                   ti(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, rough.grass.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, water.edge.length, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, builtup.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.5.04 = Sys.time() - time0
save(bam.s.t.hab.season.5.04, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.5.04"))
rm(bam.s.t.hab.season.5.04)
gc()

time0 = Sys.time()
bam.s.t.hab.season.5.05 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   ti(rough.grass.mk2) + 
                                   ti(water.edge.length) +
                                   ti(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   s(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, rough.grass.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, water.edge.length, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, combo.beach_saltmarsh.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.5.05 = Sys.time() - time0
save(bam.s.t.hab.season.5.05, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.5.05"))
rm(bam.s.t.hab.season.5.05)
gc()

time0 = Sys.time()
bam.s.t.hab.season.5.06 =  bam(unsex.capt ~ ti(julian, bs = "cc") +
                                   s(log.dens) + 
                                   sub.ses.num +
                                   s(sub.sesh.night.num, k = 4) +
                                   
                                   mink.scent +
                                   
                                   s(salat.sum.means) +
                                   coastal +
                                   s(combo.heather_all.mk2) + 
                                   ti(ALT3) +
                                   s(builtup.mk2) + 
                                   ti(rough.grass.mk2) + 
                                   ti(water.edge.length) +
                                   s(combo.beach_saltmarsh.mk2) +
                                   ti(pasture.mk2) + 
                                   ti(bog.mk2) +
                                   
                                   ti(julian, ALT3, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, pasture.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, rough.grass.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, water.edge.length, bs = c("cc", "cr"), d = c(1,1)) +
                                   ti(julian, bog.mk2, bs = c("cc", "cr"), d = c(1,1)) +
                                   
                                   s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                               data = unsex.trap , family = binomial, method = "ML")
time.bam.s.t.hab.season.5.06 = Sys.time() - time0
save(bam.s.t.hab.season.5.06, file = paste(sep = "", model.outputs.path, "bam.s.t.hab.season.5.06"))
rm(bam.s.t.hab.season.5.06)
gc()

##################################################################################################################################################################
# lvl 5 aic comparisions ----

load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.5.0")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.5.01")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.5.02")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.5.03")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.5.04")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.5.05")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.t.hab.season.5.06")

AIC(bam.s.t.hab.season.5.0, bam.s.t.hab.season.5.01, 
    bam.s.t.hab.season.5.02, bam.s.t.hab.season.5.03, bam.s.t.hab.season.5.04, 
    bam.s.t.hab.season.5.05, bam.s.t.hab.season.5.06
    )

# df      AIC
# bam.s.t.hab.season.5.0  138.5051 14927.94 ***
# bam.s.t.hab.season.5.01 137.8010 14928.11 - salat
# bam.s.t.hab.season.5.02 140.1037 14930.43 - coastal
# bam.s.t.hab.season.5.03 139.7813 14930.56 - hether
# bam.s.t.hab.season.5.04 136.8774 14931.75 - infrastucture
# bam.s.t.hab.season.5.05 142.9812 14933.33 - littoral sed
# bam.s.t.hab.season.5.06 146.4659 14941.13 -  bog



summary(bam.s.t.hab.season.5.0)
# Family: binomial 
# Link function: logit 
# 
# Formula:
#     unsex.capt ~ ti(julian, bs = "cc") + s(log.dens) + sub.ses.num + 
#     s(sub.sesh.night.num, k = 4) + mink.scent + s(salat.sum.means) + 
#     coastal + s(combo.heather_all.mk2) + ti(ALT3) + s(builtup.mk2) + 
#     ti(rough.grass.mk2) + ti(water.edge.length) + s(combo.beach_saltmarsh.mk2) + 
#     ti(pasture.mk2) + s(bog.mk2) + ti(julian, ALT3, bs = c("cc", 
#                                                            "cr"), d = c(1, 1)) + ti(julian, pasture.mk2, bs = c("cc", 
#                                                                                                                 "cr"), d = c(1, 1)) + ti(julian, rough.grass.mk2, bs = c("cc", 
#                                                                                                                                                                          "cr"), d = c(1, 1)) + ti(julian, water.edge.length, bs = c("cc", 
#                                                                                                                                                                                                                                     "cr"), d = c(1, 1)) + s(two.month.mega, bs = "re") + s(fx.trapper, 
#                                                                                                                                                                                                                                                                                            bs = "re")
# 
# Parametric coefficients:
#     Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -4.94968    0.15705 -31.516  < 2e-16 ***
#     sub.ses.num -0.47013    0.08268  -5.686  1.3e-08 ***
#     mink.scent   1.11804    0.08741  12.791  < 2e-16 ***
#     coastal     -0.03850    0.10919  -0.353    0.724    
# ---
#     Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Approximate significance of smooth terms:
#     edf  Ref.df  Chi.sq  p-value    
# ti(julian)                    1.930   3.000 231.517  0.00040 ***
#     s(log.dens)                   2.905   3.311 113.865  < 2e-16 ***
#     s(sub.sesh.night.num)         1.006   1.012  68.756  < 2e-16 ***
#     s(salat.sum.means)            2.664   3.421  27.090 2.40e-05 ***
#     s(combo.heather_all.mk2)      1.007   1.015  36.705 1.59e-09 ***
#     ti(ALT3)                      3.070   3.491  30.532 3.15e-06 ***
#     s(builtup.mk2)                5.563   6.672  44.034 1.46e-07 ***
#     ti(rough.grass.mk2)           1.018   1.035  30.270 4.21e-08 ***
#     ti(water.edge.length)         1.008   1.015   7.680  0.00596 ** 
#     s(combo.beach_saltmarsh.mk2)  1.007   1.014   8.859  0.00298 ** 
#     ti(pasture.mk2)               2.758   3.031  14.416  0.00240 ** 
#     s(bog.mk2)                    1.024   1.048   2.207  0.15378    
# ti(julian,ALT3)               2.372  12.000  61.115 4.96e-09 ***
#     ti(julian,pasture.mk2)        5.717  12.000  28.821  0.00670 ** 
#     ti(julian,rough.grass.mk2)    2.878  12.000  14.749  0.00278 ** 
#     ti(julian,water.edge.length)  2.799  12.000   8.868  0.05764 .  
# s(two.month.mega)            72.040 119.000 334.314  < 2e-16 ***
#     s(fx.trapper)                13.347  22.000  89.099 5.97e-10 ***
#     ---
#     Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# R-sq.(adj) =  0.0276   Deviance explained = 43.4%
# -ML = 2.3883e+05  Scale est. = 1         n = 175461
intrest.covar = "ALT3"
interacting = "julian" 
percentiles.inter = c(0.3, 0.6) # for julien 103 and 225

my.data = make.dummy.dat.for.2way.inter.model (model.name =bam.s.t.hab.season.2.0 , original.df = unsex.trap,
                                               interested.covar = intrest.covar,interacting.covar = interacting, 
                                               interacting.percentiles=percentiles.inter , interested.length = 100)
pred.bam.s.t.hab.season.2.0 = predict(bam.s.t.hab.season.2.0, newdata= my.data , se=TRUE, type = "response")

plot.interaciton ( model.used = bam.s.t.hab.season.2.0 , interested.covar = intrest.covar , interacting.covar = interacting,
                   interested.length = 100,prediction.object = pred.bam.s.t.hab.season.2.0, dummy.prediciton.data = my.data, 
                   orignal.data = unsex.trap , main = "alt:dens", resp.name = "TRAP.PROB",  ylimit = c(0, 0.005), interacting.percentiles = percentiles.inter)


layout(matrix(1:4,ncol=2))
plot(bam.s.t.hab.season.2.0,scheme=2)
layout(1)

