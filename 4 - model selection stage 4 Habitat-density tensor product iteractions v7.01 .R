# density effect on habtiat effects
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
##################################################################################################################################################################
# lvl 1 models -----

time0 = Sys.time()
bam.st.habseason.hab.density.1.0 = bam(unsex.capt ~ ti(julian, bs = "cc") +
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
time.bam.st.habseason.hab.density.1.0 = Sys.time() - time0
save(bam.st.habseason.hab.density.1.0, file = paste(sep = "", model.outputs.path, "bam.st.habseason.hab.density.1.0"))
rm(bam.st.habseason.hab.density.1.0)
gc()

time0 = Sys.time()
bam.st.habseason.hab.density.1.01 = bam(unsex.capt ~ ti(julian, bs = "cc") +
                                            ti(log.dens) + 
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
                                            
                                            ti(log.dens, salat.sum.means, bs = c("cr", "cr"), d = c(1,1)) +
                                            
                                            s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                                        data = unsex.trap , family = binomial, method = "ML")
time.bam.st.habseason.hab.density.1.01 = Sys.time() - time0
save(bam.st.habseason.hab.density.1.01, file = paste(sep = "", model.outputs.path, "bam.st.habseason.hab.density.1.01"))
rm(bam.st.habseason.hab.density.1.01)
gc()

time0 = Sys.time()
bam.st.habseason.hab.density.1.02 = bam(unsex.capt ~ ti(julian, bs = "cc") +
                                            ti(log.dens, by = coastal ) + 
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
time.bam.st.habseason.hab.density.1.02 = Sys.time() - time0
save(bam.st.habseason.hab.density.1.02, file = paste(sep = "", model.outputs.path, "bam.st.habseason.hab.density.1.02"))
rm(bam.st.habseason.hab.density.1.02)
gc()

time0 = Sys.time()
bam.st.habseason.hab.density.1.03 = bam(unsex.capt ~ ti(julian, bs = "cc") +
                                            ti(log.dens) + 
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
                                            
                                            ti(log.dens, combo.heather_all.mk2, bs = c("cr", "cr"), d = c(1,1)) +
                                            
                                            s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                                        data = unsex.trap , family = binomial, method = "ML")
time.bam.st.habseason.hab.density.1.03 = Sys.time() - time0
save(bam.st.habseason.hab.density.1.03, file = paste(sep = "", model.outputs.path, "bam.st.habseason.hab.density.1.03"))
rm(bam.st.habseason.hab.density.1.03)
gc()

time0 = Sys.time()
bam.st.habseason.hab.density.1.04 = bam(unsex.capt ~ ti(julian, bs = "cc") +
                                            ti(log.dens) + 
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
                                            
                                            ti(log.dens, ALT3, bs = c("cr", "cr"), d = c(1,1)) +
                                            
                                            s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                                        data = unsex.trap , family = binomial, method = "ML")
time.bam.st.habseason.hab.density.1.04 = Sys.time() - time0
save(bam.st.habseason.hab.density.1.04, file = paste(sep = "", model.outputs.path, "bam.st.habseason.hab.density.1.04"))
rm(bam.st.habseason.hab.density.1.04)
gc()

time0 = Sys.time()
bam.st.habseason.hab.density.1.05 = bam(unsex.capt ~ ti(julian, bs = "cc") +
                                            ti(log.dens) + 
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
                                            
                                            ti(log.dens, builtup.mk2, bs = c("cr", "cr"), d = c(1,1)) +
                                            
                                            s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                                        data = unsex.trap , family = binomial, method = "ML")
time.bam.st.habseason.hab.density.1.05 = Sys.time() - time0
save(bam.st.habseason.hab.density.1.05, file = paste(sep = "", model.outputs.path, "bam.st.habseason.hab.density.1.05"))
rm(bam.st.habseason.hab.density.1.05)
gc()

time0 = Sys.time()
bam.st.habseason.hab.density.1.06 = bam(unsex.capt ~ ti(julian, bs = "cc") +
                                            ti(log.dens) + 
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
                                            
                                            ti(log.dens, rough.grass.mk2, bs = c("cr", "cr"), d = c(1,1)) +
                                            
                                            s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                                        data = unsex.trap , family = binomial, method = "ML")
time.bam.st.habseason.hab.density.1.06 = Sys.time() - time0
save(bam.st.habseason.hab.density.1.06, file = paste(sep = "", model.outputs.path, "bam.st.habseason.hab.density.1.06"))
rm(bam.st.habseason.hab.density.1.06)
gc()

time0 = Sys.time()
bam.st.habseason.hab.density.1.07 = bam(unsex.capt ~ ti(julian, bs = "cc") +
                                            ti(log.dens) + 
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
                                            
                                            ti(log.dens, water.edge.length, bs = c("cr", "cr"), d = c(1,1)) +
                                            
                                            s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                                        data = unsex.trap , family = binomial, method = "ML")
time.bam.st.habseason.hab.density.1.07 = Sys.time() - time0
save(bam.st.habseason.hab.density.1.07, file = paste(sep = "", model.outputs.path, "bam.st.habseason.hab.density.1.07"))
rm(bam.st.habseason.hab.density.1.07)
gc()

time0 = Sys.time()
bam.st.habseason.hab.density.1.08 = bam(unsex.capt ~ ti(julian, bs = "cc") +
                                            ti(log.dens) + 
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
                                            
                                            ti(log.dens, combo.beach_saltmarsh.mk2, bs = c("cr", "cr"), d = c(1,1)) +
                                            
                                            s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                                        data = unsex.trap , family = binomial, method = "ML")
time.bam.st.habseason.hab.density.1.08 = Sys.time() - time0
save(bam.st.habseason.hab.density.1.08, file = paste(sep = "", model.outputs.path, "bam.st.habseason.hab.density.1.08"))
rm(bam.st.habseason.hab.density.1.08)
gc()

time0 = Sys.time()
bam.st.habseason.hab.density.1.09 = bam(unsex.capt ~ ti(julian, bs = "cc") +
                                            ti(log.dens) + 
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
                                            
                                            ti(log.dens, pasture.mk2, bs = c("cr", "cr"), d = c(1,1)) +
                                            
                                            s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                                        data = unsex.trap , family = binomial, method = "ML")
time.bam.st.habseason.hab.density.1.09 = Sys.time() - time0
save(bam.st.habseason.hab.density.1.09, file = paste(sep = "", model.outputs.path, "bam.st.habseason.hab.density.1.09"))
rm(bam.st.habseason.hab.density.1.09)
gc()

time0 = Sys.time()
bam.st.habseason.hab.density.1.10 = bam(unsex.capt ~ ti(julian, bs = "cc") +
                                            ti(log.dens) + 
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
                                            
                                            ti(log.dens, bog.mk2, bs = c("cr", "cr"), d = c(1,1)) +
                                            
                                            s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                                        data = unsex.trap , family = binomial, method = "ML")
time.bam.st.habseason.hab.density.1.10 = Sys.time() - time0
save(bam.st.habseason.hab.density.1.10, file = paste(sep = "", model.outputs.path, "bam.st.habseason.hab.density.1.10"))
rm(bam.st.habseason.hab.density.1.10)
gc()


##################################################################################################################################################################
# lvl 1 aic comparisions ----

load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.1.0")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.1.01")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.1.02")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.1.03")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.1.04")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.1.05")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.1.06")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.1.07")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.1.08")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.1.09")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.1.10")

AIC(bam.st.habseason.hab.density.1.0,
    bam.st.habseason.hab.density.1.01, 
    bam.st.habseason.hab.density.1.02, bam.st.habseason.hab.density.1.03, 
    bam.st.habseason.hab.density.1.04, bam.st.habseason.hab.density.1.05, 
    bam.st.habseason.hab.density.1.06, bam.st.habseason.hab.density.1.07, 
    bam.st.habseason.hab.density.1.08, bam.st.habseason.hab.density.1.09, bam.st.habseason.hab.density.1.10
)


# df      AIC
# bam.st.habseason.hab.density.1.0  138.5051 14927.94
# bam.st.habseason.hab.density.1.01 144.0011 14934.45
# bam.st.habseason.hab.density.1.02 146.2773 14954.66
# bam.st.habseason.hab.density.1.03 140.6487 14932.34
# bam.st.habseason.hab.density.1.04 141.6000 14915.57 ** ALT
# bam.st.habseason.hab.density.1.05 142.6583 14939.43
# bam.st.habseason.hab.density.1.06 142.8342 14930.36
# bam.st.habseason.hab.density.1.07 149.3966 14932.75
# bam.st.habseason.hab.density.1.08 142.4566 14934.97
# bam.st.habseason.hab.density.1.09 142.9670 14933.49
# bam.st.habseason.hab.density.1.10 140.0491 14931.66



##################################################################################################################################################################
# lvl 2 models ----

time0 = Sys.time()
bam.st.habseason.hab.density.2.0 = bam(unsex.capt ~ ti(julian, bs = "cc") +
                                           ti(log.dens) + 
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
                                           
                                           ti(log.dens, ALT3, bs = c("cr", "cr"), d = c(1,1)) +
                                           
                                           s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                                       data = unsex.trap , family = binomial, method = "ML")
time.bam.st.habseason.hab.density.2.0 = Sys.time() - time0
save(bam.st.habseason.hab.density.2.0, file = paste(sep = "", model.outputs.path, "bam.st.habseason.hab.density.2.0"))
rm(bam.st.habseason.hab.density.2.0)
gc()

time0 = Sys.time()
bam.st.habseason.hab.density.2.01 = bam(unsex.capt ~ ti(julian, bs = "cc") +
                                            ti(log.dens) + 
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
                                            
                                            ti(log.dens, ALT3, bs = c("cr", "cr"), d = c(1,1)) +
                                            ti(log.dens, salat.sum.means, bs = c("cr", "cr"), d = c(1,1)) +
                                            
                                            s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                                        data = unsex.trap , family = binomial, method = "ML")
time.bam.st.habseason.hab.density.2.01 = Sys.time() - time0
save(bam.st.habseason.hab.density.2.01, file = paste(sep = "", model.outputs.path, "bam.st.habseason.hab.density.2.01"))
rm(bam.st.habseason.hab.density.2.01)
gc()

time0 = Sys.time()
bam.st.habseason.hab.density.2.02 = bam(unsex.capt ~ ti(julian , bs = "cc") +
                                            ti(log.dens, by = coastal) + 
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
                                            
                                            ti(log.dens, ALT3, bs = c("cr", "cr"), d = c(1,1)) +

                                            s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                                        data = unsex.trap , family = binomial, method = "ML")
time.bam.st.habseason.hab.density.2.02 = Sys.time() - time0
save(bam.st.habseason.hab.density.2.02, file = paste(sep = "", model.outputs.path, "bam.st.habseason.hab.density.2.02"))
rm(bam.st.habseason.hab.density.2.02)
gc()

time0 = Sys.time()
bam.st.habseason.hab.density.2.03 = bam(unsex.capt ~ ti(julian, bs = "cc") +
                                            ti(log.dens) + 
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
                                            
                                            ti(log.dens, ALT3, bs = c("cr", "cr"), d = c(1,1)) +
                                            ti(log.dens, combo.heather_all.mk2, bs = c("cr", "cr"), d = c(1,1)) +
                                            
                                            s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                                        data = unsex.trap , family = binomial, method = "ML")
time.bam.st.habseason.hab.density.2.03 = Sys.time() - time0
save(bam.st.habseason.hab.density.2.03, file = paste(sep = "", model.outputs.path, "bam.st.habseason.hab.density.2.03"))
rm(bam.st.habseason.hab.density.2.03)
gc()

time0 = Sys.time()
bam.st.habseason.hab.density.2.04 = bam(unsex.capt ~ ti(julian, bs = "cc") +
                                            ti(log.dens) + 
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
                                            
                                            ti(log.dens, ALT3, bs = c("cr", "cr"), d = c(1,1)) +
                                            ti(log.dens, builtup.mk2, bs = c("cr", "cr"), d = c(1,1)) +
                                            
                                            s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                                        data = unsex.trap , family = binomial, method = "ML")
time.bam.st.habseason.hab.density.2.04 = Sys.time() - time0
save(bam.st.habseason.hab.density.2.04, file = paste(sep = "", model.outputs.path, "bam.st.habseason.hab.density.2.04"))
rm(bam.st.habseason.hab.density.2.04)
gc()

time0 = Sys.time()
bam.st.habseason.hab.density.2.05 = bam(unsex.capt ~ ti(julian, bs = "cc") +
                                            ti(log.dens) + 
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
                                            
                                            ti(log.dens, ALT3, bs = c("cr", "cr"), d = c(1,1)) +
                                            ti(log.dens, rough.grass.mk2, bs = c("cr", "cr"), d = c(1,1)) +
                                            
                                            s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                                        data = unsex.trap , family = binomial, method = "ML")
time.bam.st.habseason.hab.density.2.05 = Sys.time() - time0
save(bam.st.habseason.hab.density.2.05, file = paste(sep = "", model.outputs.path, "bam.st.habseason.hab.density.2.05"))
rm(bam.st.habseason.hab.density.2.05)
gc()

time0 = Sys.time()
bam.st.habseason.hab.density.2.06 = bam(unsex.capt ~ ti(julian, bs = "cc") +
                                            ti(log.dens) + 
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
                                            
                                            ti(log.dens, ALT3, bs = c("cr", "cr"), d = c(1,1)) +
                                            ti(log.dens, water.edge.length, bs = c("cr", "cr"), d = c(1,1)) +
                                            
                                            s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                                        data = unsex.trap , family = binomial, method = "ML")
time.bam.st.habseason.hab.density.2.06 = Sys.time() - time0
save(bam.st.habseason.hab.density.2.06, file = paste(sep = "", model.outputs.path, "bam.st.habseason.hab.density.2.06"))
rm(bam.st.habseason.hab.density.2.06)
gc()

time0 = Sys.time()
bam.st.habseason.hab.density.2.07 = bam(unsex.capt ~ ti(julian, bs = "cc") +
                                            ti(log.dens) + 
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
                                            
                                            ti(log.dens, ALT3, bs = c("cr", "cr"), d = c(1,1)) +
                                            ti(log.dens, combo.beach_saltmarsh.mk2, bs = c("cr", "cr"), d = c(1,1)) +
                                            
                                            s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                                        data = unsex.trap , family = binomial, method = "ML")
time.bam.st.habseason.hab.density.2.07 = Sys.time() - time0
save(bam.st.habseason.hab.density.2.07, file = paste(sep = "", model.outputs.path, "bam.st.habseason.hab.density.2.07"))
rm(bam.st.habseason.hab.density.2.07)
gc()

time0 = Sys.time()
bam.st.habseason.hab.density.2.08 = bam(unsex.capt ~ ti(julian, bs = "cc") +
                                            ti(log.dens) + 
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
                                            
                                            ti(log.dens, ALT3, bs = c("cr", "cr"), d = c(1,1)) +
                                            ti(log.dens, pasture.mk2, bs = c("cr", "cr"), d = c(1,1)) +
                                            
                                            s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                                        data = unsex.trap , family = binomial, method = "ML")
time.bam.st.habseason.hab.density.2.08 = Sys.time() - time0
save(bam.st.habseason.hab.density.2.08, file = paste(sep = "", model.outputs.path, "bam.st.habseason.hab.density.2.08"))
rm(bam.st.habseason.hab.density.2.08)
gc()

time0 = Sys.time()
bam.st.habseason.hab.density.2.09 = bam(unsex.capt ~ ti(julian, bs = "cc") +
                                            ti(log.dens) + 
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
                                            
                                            ti(log.dens, ALT3, bs = c("cr", "cr"), d = c(1,1)) +
                                            ti(log.dens, bog.mk2, bs = c("cr", "cr"), d = c(1,1)) +
                                            
                                            s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                                        data = unsex.trap , family = binomial, method = "ML")
time.bam.st.habseason.hab.density.2.09 = Sys.time() - time0
save(bam.st.habseason.hab.density.2.09, file = paste(sep = "", model.outputs.path, "bam.st.habseason.hab.density.2.09"))
rm(bam.st.habseason.hab.density.2.09)
gc()

##################################################################################################################################################################
# lvl 2 aic comparisions ----

load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.2.0")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.2.01")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.2.02")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.2.03")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.2.04")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.2.05")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.2.06")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.2.07")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.2.08")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.st.habseason.hab.density.2.09")

AIC(bam.st.habseason.hab.density.2.0, bam.st.habseason.hab.density.2.01, 
    bam.st.habseason.hab.density.2.02, 
    bam.st.habseason.hab.density.2.03, bam.st.habseason.hab.density.2.04, 
    bam.st.habseason.hab.density.2.05, bam.st.habseason.hab.density.2.06, 
    bam.st.habseason.hab.density.2.07, 
    bam.st.habseason.hab.density.2.08, bam.st.habseason.hab.density.2.09
)

# df      AIC
# bam.st.habseason.hab.density.2.0  141.6000 14915.57 *
# bam.st.habseason.hab.density.2.01 141.2259 14915.75 * salat
# bam.st.habseason.hab.density.2.02 141.7194 14914.81 * coastal - right?
# bam.st.habseason.hab.density.2.03 141.8132 14915.34 * heather
# bam.st.habseason.hab.density.2.04 144.4447 14926.83 - builtup
# bam.st.habseason.hab.density.2.05 149.7874 14918.05 - rough.grass
# bam.st.habseason.hab.density.2.06 149.6789 14919.78 - water.edge.length
# bam.st.habseason.hab.density.2.07 141.7534 14915.75 * .beach_saltmarsh
# bam.st.habseason.hab.density.2.08 141.3869 14914.91 * pasture
# bam.st.habseason.hab.density.2.09 142.3761 14917.31 - bog


##################################################################################################################################################################
# fianl model as a gam




time0 = Sys.time()
gam.st.habseason.hab.density.2.0 = gam(unsex.capt ~ ti(julian, bs = "cc") +
                                           ti(log.dens) + 
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
                                           
                                           ti(log.dens, ALT3, bs = c("cr", "cr"), d = c(1,1)) +
                                           
                                           s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                                       data = unsex.trap , family = binomial, method = "ML")
time.gam.st.habseason.hab.density.2.0 = Sys.time() - time0
save(gam.st.habseason.hab.density.2.0, file = paste(sep = "", model.outputs.path, "gam.st.habseason.hab.density.2.0"))
rm(gam.st.habseason.hab.density.2.0)
gc()
##################################################################################################################################################################
plot(bam.st.habseason.hab.density.2.0)

# exploring best model
habdens.sum = summary (bam.st.habseason.hab.density.2.0) #----
# 
# 
# Family: binomial 
# Link function: logit 
# 
# Formula:
#     sex.cap ~ ti(julian, bs = "cc") + female + s(day.since.start) + 
#     s(block.dens) + sub.ses.num + s(sub.sesh.night.num, k = 4) + 
#     s(depletion.sex) + s(traps.open.within.buffer) + tern.col + 
#     boat.access + othersex.scent + samesex.scent + s(salat.sum.means) + 
#     s(combo.heather_all.mk2) + ti(coast.dist) + ti(ALT3) + ti(rough.grass.mk2) + 
#     s(builtup.mk2) + ti(combo.bog_all.mk2) + ti(pasture.mk2) + 
#     s(woodland.mk2) + ti(julian, combo.bog_all.mk2, bs = c("cc", 
#       "cr"), d = c(1, 1)) + ti(julian, pasture.mk2, bs = c("cc", 
#       "cr"), d = c(1, 1)) + ti(julian, coast.dist, bs = c("cc", 
# 
# Parametric coefficients:
#     Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -5.67855    0.13604 -41.741  < 2e-16 ***
#     female          0.09460    0.05258   1.799  0.07198 .  
# sub.ses.num    -0.55948    0.08610  -6.498 8.15e-11 ***
#     tern.col        1.02718    0.35714   2.876  0.00403 ** 
#     boat.access     0.02215    0.09463   0.234  0.81490    
# othersex.scent  1.00388    0.10568   9.500  < 2e-16 ***
#     samesex.scent   0.63657    0.11837   5.378 7.54e-08 ***
#     ---
#     Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Approximate significance of smooth terms:
#     edf  Ref.df  Chi.sq  p-value    
# ti(julian)                        2.3280   3.000 370.221 6.93e-09 ***
#     s(day.since.start)                1.0123   1.017  65.095 7.18e-16 ***
#     s(block.dens)                     1.0018   1.003   4.548 0.033295 *  
#     s(sub.sesh.night.num)             1.0036   1.007  70.143  < 2e-16 ***
#     s(depletion.sex)                  4.7436   5.666  42.760 2.00e-07 ***
#     s(traps.open.within.buffer)       3.1973   3.980  52.878 7.95e-11 ***
#     s(salat.sum.means)                2.6796   3.447  22.161 9.46e-05 ***
#     s(combo.heather_all.mk2)          1.0044   1.009  28.710 8.98e-08 ***
#     ti(coast.dist)                    2.9985   3.261  28.647 1.68e-05 ***
#     ti(ALT3)                          2.6715   3.136  15.497 0.001726 ** 
#     ti(rough.grass.mk2)               3.4686   3.780  35.186 4.52e-07 ***
#     s(builtup.mk2)                    5.5826   6.689  38.872 1.36e-06 ***
#     ti(combo.bog_all.mk2)             2.9732   3.466  19.291 0.000540 ***
#     ti(pasture.mk2)                   2.3726   2.662  18.575 0.000521 ***
#     s(woodland.mk2)                   1.0054   1.011   6.725 0.009876 ** 
#     ti(julian,combo.bog_all.mk2)      0.0267  12.000   0.018 0.640527    
# ti(julian,pasture.mk2)            5.0150  12.000  21.941 0.002343 ** 
#     ti(julian,coast.dist)             4.5828  12.000  22.943 0.011033 *  
#     ti(julian,ALT3)                   2.4082  12.000  39.338 4.55e-08 ***
#     ti(julian,rough.grass.mk2)        2.5714  12.000   9.183 0.009881 ** 
#     ti(block.dens,combo.bog_all.mk2)  1.9551   2.248   5.943 0.062029 .  
# s(two.month.mega)                62.3980 119.000 209.922  < 2e-16 ***
#     s(fx.trapper)                    12.8246  22.000  77.899 7.67e-11 ***
#     ---
#     Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# R-sq.(adj) =  0.0154   Deviance explained =   52%
# -ML = 4.7756e+05  Scale est. = 1         n = 350922

####
habdens.CC = concurvity(bam.st.habseason.hab.density.2.0) # ----
# para ti(julian) s(day.since.start) s(block.dens) s(sub.sesh.night.num) s(depletion.sex) s(traps.open.within.buffer) s(salat.sum.means) s(combo.heather_all.mk2) ti(coast.dist)
# worst       1  0.9720355          0.9995001     0.9940860            0.06433633        0.3576757                   0.3840586          0.7351537                0.6875697      0.8342401
# observed    1  0.9694773          0.9994397     0.9354596            0.02717903        0.2006180                   0.3415670          0.6806074                0.5706993      0.6148741
# estimate    1  0.9363928          0.9989275     0.9318035            0.02667047        0.3428416                   0.3121726          0.3648490                0.4428247      0.4514291

# ti(ALT3) ti(rough.grass.mk2) s(builtup.mk2) ti(combo.bog_all.mk2) ti(pasture.mk2) s(woodland.mk2) ti(julian,combo.bog_all.mk2) ti(julian,pasture.mk2) ti(julian,coast.dist)
# worst    0.6394971           0.6961806      0.6784027             0.8294449       0.7448941       0.2761109                    0.7734467              0.6154710             0.8354126
# observed 0.6229077           0.5042185      0.4182191             0.5635095       0.6932692       0.2192885                    0.5301555              0.4357916             0.4095448
# estimate 0.3478668           0.3809252      0.3188010             0.4453794       0.6201699       0.1748253                    0.3543839              0.3730359             0.3714189

# ti(julian,ALT3) ti(julian,rough.grass.mk2) ti(block.dens,combo.bog_all.mk2) s(two.month.mega) s(fx.trapper)
# worst          0.5598269                  0.6699591                        0.6196601         1.0000000    1.00000000
# observed       0.4456562                  0.5487273                        0.4738271         0.5134346    0.09161245
# estimate       0.2851568                  0.3168585                        0.2251995         0.3648566    0.19024206
####
vcov.habdens = vcov(bam.st.habseason.hab.density.2.0) # ----
####
cov2cor(vcov.habdens) # ------
# no corr >0.9 or -0.9 not form same spline
#####
# plots ----
# main habitat effects interactions ----
#     s(salat.sum.means)                2.6796   3.447  22.161 9.46e-05 *** ----
intrest.covar = "salat.sum.means"
expl.lab =   "Sheltered Coast Index"
ylims = c(0 , 0.015)

interested.length = 200

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                     interested.covar = intrest.covar, interested.length = interested.length, 
                                     by.fac = F, interesting.factor.name = NULL)
my.data = my.data1
pred.bam.st.habseason.hab.density.2.0 = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")
my.data$fit = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit))
my.data$upper.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit + 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))
my.data$lower.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit - 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))

my.data$expl = my.data[,which(names(my.data)==intrest.covar )]
rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])
##ylims =  c(0,max( my.data$upper.95 * 1.1))
 
ggplot()+
    geom_ribbon(aes(ymin = my.data$lower.95, ymax = my.data$upper.95, x = my.data$expl),
                alpha = 0.4, fill = "#afcbff") +
    geom_line(aes(y = my.data$fit, x = my.data$expl), size = 1, col = "black") +
    geom_rug( aes(x = rugs), alpha = 0.5 ,sides = "b",# position = "jitter",
              size = 0.1) +
    labs(x = expl.lab, y = "Target Sex Capture Probability") +
    scale_y_continuous(limits=ylims) +
    scale_x_continuous(limits=c(0,max(my.data$expl))) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18, face = "bold")
    )
 
#     s(combo.heather_all.mk2)          1.0044   1.009  28.710 8.98e-08 *** ----
intrest.covar = "combo.heather_all.mk2"
expl.lab =   "Heather Area (km2)"
divider = 1000000
#ylims = NA

interested.length = 200

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                     interested.covar = intrest.covar, interested.length = interested.length, 
                                     by.fac = F, interesting.factor.name = NULL)
my.data = my.data1
pred.bam.st.habseason.hab.density.2.0 = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")
my.data$fit = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit))
my.data$upper.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit + 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))
my.data$lower.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit - 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))

my.data$expl = my.data[,which(names(my.data)==intrest.covar )]/divider
rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])/divider
#ylims =  c(0,max( my.data$upper.95 * 1.1))

ggplot()+
    geom_ribbon(aes(ymin = my.data$lower.95, ymax = my.data$upper.95, x = my.data$expl),
                alpha = 0.4, fill = "#afcbff") +
    geom_line(aes(y = my.data$fit, x = my.data$expl), size = 1, col = "black") +
    geom_rug( aes(x = rugs), alpha = 0.5 ,sides = "b",# position = "jitter",
              size = 0.2) +
    labs(x = expl.lab, y = "Target Sex Capture Probability") +
    scale_y_continuous(limits=ylims) +
    scale_x_continuous(limits=c(0,max(my.data$expl))) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18, face = "bold")
    )


#     s(builtup.mk2)                    5.5826   6.689  38.872 1.36e-06 *** ----
intrest.covar = "builtup.mk2"
expl.lab =   "Urban Area (km2)"
divider = 1000000
#ylims = NA

interested.length = 200

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                     interested.covar = intrest.covar, interested.length = interested.length, 
                                     by.fac = F, interesting.factor.name = NULL)
my.data = my.data1
pred.bam.st.habseason.hab.density.2.0 = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")
my.data$fit = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit))
my.data$upper.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit + 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))
my.data$lower.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit - 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))

my.data$expl = my.data[,which(names(my.data)==intrest.covar )]/divider
rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])/divider
#ylims =  c(0,max( my.data$upper.95 * 1.1))

ggplot()+
    geom_ribbon(aes(ymin = my.data$lower.95, ymax = my.data$upper.95, x = my.data$expl),
                alpha = 0.4, fill = "#afcbff") +
    geom_line(aes(y = my.data$fit, x = my.data$expl), size = 1, col = "black") +
    geom_rug( aes(x = rugs), alpha = 0.5 ,sides = "b",# position = "jitter",
              size = 0.2) +
    labs(x = expl.lab, y = "Target Sex Capture Probability") +
    scale_y_continuous(limits=ylims) +
    scale_x_continuous(limits=c(0,max(my.data$expl))) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18, face = "bold")
    )



#     s(woodland.mk2)                   1.0054   1.011   6.725 0.009876 **  -----
intrest.covar = "woodland.mk2"
expl.lab =   "Woodland Area (km2)"
divider = 1000000
#ylims = NA

interested.length = 200

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                     interested.covar = intrest.covar, interested.length = interested.length, 
                                     by.fac = F, interesting.factor.name = NULL)
my.data = my.data1
pred.bam.st.habseason.hab.density.2.0 = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")
my.data$fit = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit))
my.data$upper.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit + 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))
my.data$lower.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit - 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))

my.data$expl = my.data[,which(names(my.data)==intrest.covar )]/divider
rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])/divider
#ylims =  c(0,max( my.data$upper.95 * 1.1))

ggplot()+
    geom_ribbon(aes(ymin = my.data$lower.95, ymax = my.data$upper.95, x = my.data$expl),
                alpha = 0.4, fill = "#afcbff") +
    geom_line(aes(y = my.data$fit, x = my.data$expl), size = 1, col = "black") +
    geom_rug( aes(x = rugs), alpha = 0.5 ,sides = "b",# position = "jitter",
              size = 0.2) +
    labs(x = expl.lab, y = "Target Sex Capture Probability") +
    scale_y_continuous(limits=ylims) +
    scale_x_continuous(limits=c(0,max(my.data$expl))) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18, face = "bold")
    )


###
#####################################################################################
# seasonal habitat interactions ----
#ti(julian,combo.bog_all.mk2)      0.0267  12.000   0.018 0.640527     ----

intrest.covar = "combo.bog_all.mk2"
itneracting = "julian"
y.lab = "Day of Year"
x.lab = "Bog Area (km2)"
divider = 1000000

interested.length = 100
min.inter = min (unsex.trap[,names(unsex.trap)==itneracting ])
max.inter = max (unsex.trap[,names(unsex.trap)==itneracting ])

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                     interested.covar = intrest.covar, interested.length = interested.length, 
                                     by.fac = F, interesting.factor.name = NULL)
my.data = my.data1[rep(1:(dim(my.data1)[1]), times = interested.length ),]
my.data[, which(names(my.data)== itneracting)] = rep (seq(min.inter, max.inter, length = interested.length), each = interested.length)

pred.bam.st.habseason.hab.density.2.0 = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")

plotr =my.data[, which(names(my.data) %in% c(intrest.covar, itneracting ))]
plotr$fit = as.numeric( plogis(pred.bam.st.habseason.hab.density.2.0$fit))
plotr$upper.95 = as.numeric( plogis(pred.bam.st.habseason.hab.density.2.0$fit + 1.96 * pred.bam.st.habseason.hab.density.2.0$se.fit))
plotr$lower.95 = as.numeric( plogis(pred.bam.st.habseason.hab.density.2.0$fit - 1.96 * pred.bam.st.habseason.hab.density.2.0$se.fit))
plotr$expl = as.numeric(my.data[,which(names(my.data) %in% c(intrest.covar))]/divider)
plotr$intr = as.numeric(my.data[,which(names(my.data) == (itneracting))])
rug = unique(unsex.trap[,which(names(unsex.trap) == intrest.covar)])/divider
caps= unsex.trap[unsex.trap$sex.cap==1,which(names(unsex.trap) %in% c(intrest.covar, itneracting))]
caps$itneracting= unsex.trap[unsex.trap$sex.cap==1,which(names(unsex.trap) %in% c( itneracting))]
caps$intrest.covar= unsex.trap[unsex.trap$sex.cap==1,which(names(unsex.trap) %in% c( intrest.covar))]/divider

#main.plot = 
ggplot() + 
    geom_tile(aes(plotr$expl, plotr$intr, fill = plotr$fit)) +
    geom_point(aes(caps$intrest.covar, caps$itneracting), color = "white", alpha = 0.5)+
    geom_rug(aes(x = rug), alpha = 1/2,sides = "b",size = 0.2)+
    labs(y = y.lab, x = x.lab , z = "Target Sex Capture Probability") +
    geom_contour(aes(plotr$expl, plotr$intr, z = plotr$fit))+
    scale_y_continuous(breaks= c(1/12, 3/12, 5/12, 7/12, 9/12,11/12)*365, labels= c("Jan", "March", "May", "July","Sept", "Nov"))+
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18, face = "bold")
    )
####
# ti(julian,pasture.mk2)            5.0150  12.000  21.941 0.002343 **  ----

intrest.covar = "pasture.mk2"
itneracting = "julian"
y.lab = "Day of Year"
x.lab = "Pasture Area (km2)"
divider = 1000000

interested.length = 100
min.inter = min (unsex.trap[,names(unsex.trap)==itneracting ])
max.inter = max (unsex.trap[,names(unsex.trap)==itneracting ])

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                     interested.covar = intrest.covar, interested.length = interested.length, 
                                     by.fac = F, interesting.factor.name = NULL)
my.data = my.data1[rep(1:(dim(my.data1)[1]), times = interested.length ),]
my.data[, which(names(my.data)== itneracting)] = rep (seq(min.inter, max.inter, length = interested.length), each = interested.length)

pred.bam.st.habseason.hab.density.2.0 = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")

plotr =my.data[, which(names(my.data) %in% c(intrest.covar, itneracting ))]
plotr$fit = as.numeric( plogis(pred.bam.st.habseason.hab.density.2.0$fit))
plotr$upper.95 = as.numeric( plogis(pred.bam.st.habseason.hab.density.2.0$fit + 1.96 * pred.bam.st.habseason.hab.density.2.0$se.fit))
plotr$lower.95 = as.numeric( plogis(pred.bam.st.habseason.hab.density.2.0$fit - 1.96 * pred.bam.st.habseason.hab.density.2.0$se.fit))
plotr$expl = as.numeric(my.data[,which(names(my.data) %in% c(intrest.covar))]/divider)
plotr$intr = as.numeric(my.data[,which(names(my.data) == (itneracting))])
rug = unique(unsex.trap[,which(names(unsex.trap) == intrest.covar)])/divider
caps= unsex.trap[unsex.trap$sex.cap==1,which(names(unsex.trap) %in% c(intrest.covar, itneracting))]
caps$itneracting= unsex.trap[unsex.trap$sex.cap==1,which(names(unsex.trap) %in% c( itneracting))]
caps$intrest.covar= unsex.trap[unsex.trap$sex.cap==1,which(names(unsex.trap) %in% c( intrest.covar))]/divider

#main.plot = 
ggplot() + 
    geom_tile(aes(plotr$expl, plotr$intr, fill = plotr$fit)) +
    geom_point(aes(caps$intrest.covar, caps$itneracting), color = "white", alpha = 0.5)+
    geom_rug(aes(x = rug), alpha = 1/2,sides = "b",size = 0.2)+
    labs(y = y.lab, x = x.lab , z = "Target Sex Capture Probability") +
    geom_contour(aes(plotr$expl, plotr$intr, z = plotr$fit))+
    scale_y_continuous(breaks= c(1/12, 3/12, 5/12, 7/12, 9/12,11/12)*365, labels= c("Jan", "March", "May", "July","Sept", "Nov"))+
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18, face = "bold")
    )
#     ti(julian,coast.dist)             4.5828  12.000  22.943 0.011033 * ----  

intrest.covar = "coast.dist"
itneracting = "julian"
y.lab = "Day of Year"
x.lab = "Distance to Coast (km)"
divider = 1000

interested.length = 100
min.inter = min (unsex.trap[,names(unsex.trap)==itneracting ])
max.inter = max (unsex.trap[,names(unsex.trap)==itneracting ])

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                     interested.covar = intrest.covar, interested.length = interested.length, 
                                     by.fac = F, interesting.factor.name = NULL)
my.data = my.data1[rep(1:(dim(my.data1)[1]), times = interested.length ),]
my.data[, which(names(my.data)== itneracting)] = rep (seq(min.inter, max.inter, length = interested.length), each = interested.length)

pred.bam.st.habseason.hab.density.2.0 = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")

plotr =my.data[, which(names(my.data) %in% c(intrest.covar, itneracting ))]
plotr$fit = as.numeric( plogis(pred.bam.st.habseason.hab.density.2.0$fit))
plotr$upper.95 = as.numeric( plogis(pred.bam.st.habseason.hab.density.2.0$fit + 1.96 * pred.bam.st.habseason.hab.density.2.0$se.fit))
plotr$lower.95 = as.numeric( plogis(pred.bam.st.habseason.hab.density.2.0$fit - 1.96 * pred.bam.st.habseason.hab.density.2.0$se.fit))
plotr$expl = as.numeric(my.data[,which(names(my.data) %in% c(intrest.covar))]/divider)
plotr$intr = as.numeric(my.data[,which(names(my.data) == (itneracting))])
rug = unique(unsex.trap[,which(names(unsex.trap) == intrest.covar)])/divider
caps= unsex.trap[unsex.trap$sex.cap==1,which(names(unsex.trap) %in% c(intrest.covar, itneracting))]
caps$itneracting= unsex.trap[unsex.trap$sex.cap==1,which(names(unsex.trap) %in% c( itneracting))]
caps$intrest.covar= unsex.trap[unsex.trap$sex.cap==1,which(names(unsex.trap) %in% c( intrest.covar))]/divider

#main.plot = 
ggplot() + 
    geom_tile(aes(plotr$expl, plotr$intr, fill = plotr$fit)) +
    geom_point(aes(caps$intrest.covar, caps$itneracting), color = "white", alpha = 0.5)+
    geom_rug(aes(x = rug), alpha = 1/2,sides = "b",size = 0.2)+
    labs(y = y.lab, x = x.lab , z = "Target Sex Capture Probability") +
    geom_contour(aes(plotr$expl, plotr$intr, z = plotr$fit))+
    scale_y_continuous(breaks= c(1/12, 3/12, 5/12, 7/12, 9/12,11/12)*365, labels= c("Jan", "March", "May", "July","Sept", "Nov"))+
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18, face = "bold")
    )
#     ti(julian,ALT3)                   2.4082  12.000  39.338 4.55e-08 *** ----

intrest.covar = "ALT3"
itneracting = "julian"
y.lab = "Day of Year"
x.lab = "Altitude (m)"
divider = 1

interested.length = 100
min.inter = min (unsex.trap[,names(unsex.trap)==itneracting ])
max.inter = max (unsex.trap[,names(unsex.trap)==itneracting ])

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                     interested.covar = intrest.covar, interested.length = interested.length, 
                                     by.fac = F, interesting.factor.name = NULL)
my.data = my.data1[rep(1:(dim(my.data1)[1]), times = interested.length ),]
my.data[, which(names(my.data)== itneracting)] = rep (seq(min.inter, max.inter, length = interested.length), each = interested.length)

pred.bam.st.habseason.hab.density.2.0 = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")

plotr =my.data[, which(names(my.data) %in% c(intrest.covar, itneracting ))]
plotr$fit = as.numeric( plogis(pred.bam.st.habseason.hab.density.2.0$fit))
plotr$upper.95 = as.numeric( plogis(pred.bam.st.habseason.hab.density.2.0$fit + 1.96 * pred.bam.st.habseason.hab.density.2.0$se.fit))
plotr$lower.95 = as.numeric( plogis(pred.bam.st.habseason.hab.density.2.0$fit - 1.96 * pred.bam.st.habseason.hab.density.2.0$se.fit))
plotr$expl = as.numeric(my.data[,which(names(my.data) %in% c(intrest.covar))]/divider)
plotr$intr = as.numeric(my.data[,which(names(my.data) == (itneracting))])
rug = unique(unsex.trap[,which(names(unsex.trap) == intrest.covar)])/divider
caps= unsex.trap[unsex.trap$sex.cap==1,which(names(unsex.trap) %in% c(intrest.covar, itneracting))]
caps$itneracting= unsex.trap[unsex.trap$sex.cap==1,which(names(unsex.trap) %in% c( itneracting))]
caps$intrest.covar= unsex.trap[unsex.trap$sex.cap==1,which(names(unsex.trap) %in% c( intrest.covar))]/divider

#main.plot = 
ggplot() + 
    geom_tile(aes(plotr$expl, plotr$intr, fill = plotr$fit)) +
    geom_point(aes(caps$intrest.covar, caps$itneracting), color = "white", alpha = 0.5)+
    geom_rug(aes(x = rug), alpha = 1/2,sides = "b",size = 0.2)+
    labs(y = y.lab, x = x.lab , z = "Target Sex Capture Probability") +
    geom_contour(aes(plotr$expl, plotr$intr, z = plotr$fit))+
    scale_y_continuous(breaks= c(1/12, 3/12, 5/12, 7/12, 9/12,11/12)*365, labels= c("Jan", "March", "May", "July","Sept", "Nov"))+
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18, face = "bold")
    )
#     ti(julian,rough.grass.mk2)        2.5714  12.000   9.183 0.009881 **  ----

intrest.covar = "rough.grass.mk2"
itneracting = "julian"
y.lab = "Day of Year"
x.lab = "Rough Grass Area (km2)"
divider = 1000000

interested.length = 100
min.inter = min (unsex.trap[,names(unsex.trap)==itneracting ])
max.inter = max (unsex.trap[,names(unsex.trap)==itneracting ])

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                     interested.covar = intrest.covar, interested.length = interested.length, 
                                     by.fac = F, interesting.factor.name = NULL)
my.data = my.data1[rep(1:(dim(my.data1)[1]), times = interested.length ),]
my.data[, which(names(my.data)== itneracting)] = rep (seq(min.inter, max.inter, length = interested.length), each = interested.length)

pred.bam.st.habseason.hab.density.2.0 = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")

plotr =my.data[, which(names(my.data) %in% c(intrest.covar, itneracting ))]
plotr$fit = as.numeric( plogis(pred.bam.st.habseason.hab.density.2.0$fit))
plotr$upper.95 = as.numeric( plogis(pred.bam.st.habseason.hab.density.2.0$fit + 1.96 * pred.bam.st.habseason.hab.density.2.0$se.fit))
plotr$lower.95 = as.numeric( plogis(pred.bam.st.habseason.hab.density.2.0$fit - 1.96 * pred.bam.st.habseason.hab.density.2.0$se.fit))
plotr$expl = as.numeric(my.data[,which(names(my.data) %in% c(intrest.covar))]/divider)
plotr$intr = as.numeric(my.data[,which(names(my.data) == (itneracting))])
rug = unique(unsex.trap[,which(names(unsex.trap) == intrest.covar)])/divider
caps= unsex.trap[unsex.trap$sex.cap==1,which(names(unsex.trap) %in% c(intrest.covar, itneracting))]
caps$itneracting= unsex.trap[unsex.trap$sex.cap==1,which(names(unsex.trap) %in% c( itneracting))]
caps$intrest.covar= unsex.trap[unsex.trap$sex.cap==1,which(names(unsex.trap) %in% c( intrest.covar))]/divider

#main.plot = 
ggplot() + 
    geom_tile(aes(plotr$expl, plotr$intr, fill = plotr$fit)) +
    geom_point(aes(caps$intrest.covar, caps$itneracting), color = "white", alpha = 0.5)+
    geom_rug(aes(x = rug), alpha = 1/2,sides = "b",size = 0.2)+
    labs(y = y.lab, x = x.lab , z = "Target Sex Capture Probability") +
    geom_contour(aes(plotr$expl, plotr$intr, z = plotr$fit))+
    scale_y_continuous(breaks= c(1/12, 3/12, 5/12, 7/12, 9/12,11/12)*365, labels= c("Jan", "March", "May", "July","Sept", "Nov"))+
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18, face = "bold")
    )
#####################################################################################
# density habitat interactions ----
#     ti(block.dens,combo.bog_all.mk2)  1.9551   2.248   5.943 0.062029 .  
# plot heatmap -----
intrest.covar = "combo.bog_all.mk2"
itneracting = "block.dens"
y.lab = "Estimated Population Size"
x.lab = "Bog Area (km2)"
divider = 1000000

interested.length = 100
min.inter = min (unsex.trap[,names(unsex.trap)==itneracting ])
max.inter = max (unsex.trap[,names(unsex.trap)==itneracting ])

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                     interested.covar = intrest.covar, interested.length = interested.length, 
                                     by.fac = F, interesting.factor.name = NULL)
my.data = my.data1[rep(1:(dim(my.data1)[1]), times = interested.length ),]
my.data[, which(names(my.data)== itneracting)] = rep (seq(min.inter, max.inter, length = interested.length), each = interested.length)
#my.data$day.since.start = rep (seq(max(unsex.trap$day.since.start), min(unsex.trap$day.since.start), length = interested.length), each = interested.length)


pred.bam.st.habseason.hab.density.2.0 = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")

plotr =my.data[, which(names(my.data) %in% c(intrest.covar, itneracting ))]
plotr$fit = as.numeric( plogis(pred.bam.st.habseason.hab.density.2.0$fit))
plotr$upper.95 = as.numeric( plogis(pred.bam.st.habseason.hab.density.2.0$fit + 1.96 * pred.bam.st.habseason.hab.density.2.0$se.fit))
plotr$lower.95 = as.numeric( plogis(pred.bam.st.habseason.hab.density.2.0$fit - 1.96 * pred.bam.st.habseason.hab.density.2.0$se.fit))
plotr$expl = as.numeric(my.data[,which(names(my.data) %in% c(intrest.covar))]/divider)
plotr$intr = as.numeric(my.data[,which(names(my.data) == (itneracting))])
rug = unique(unsex.trap[,which(names(unsex.trap) == intrest.covar)])/divider
caps= unsex.trap[unsex.trap$sex.cap==1,which(names(unsex.trap) %in% c(intrest.covar, itneracting))]
caps$itneracting= unsex.trap[unsex.trap$sex.cap==1,which(names(unsex.trap) %in% c( itneracting))]
caps$intrest.covar= unsex.trap[unsex.trap$sex.cap==1,which(names(unsex.trap) %in% c( intrest.covar))]/divider

#main.plot = 
ggplot() + 
    geom_tile(aes(plotr$expl, plotr$intr, fill = plotr$fit)) +
    geom_point(aes(caps$intrest.covar, caps$itneracting), color = "white", alpha = 0.5)+
    geom_rug(aes(x = rug), alpha = 1/2,sides = "b",size = 0.2)+
    labs(y = y.lab, x = x.lab , z = "Target Sex Capture Probability") +
    geom_contour(aes(plotr$expl, plotr$intr, z = plotr$fit))+
    scale_y_continuous(breaks= seq(0,2000,500), labels= seq(0,2000,500))+
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18, face = "bold")
    )
# plot line for high and low density -----
plotr$link.fit = as.numeric( (pred.bam.st.habseason.hab.density.2.0$fit))
plotr$link.upper.95 = as.numeric( (pred.bam.st.habseason.hab.density.2.0$fit + 1.96 * pred.bam.st.habseason.hab.density.2.0$se.fit))
plotr$link.lower.95 = as.numeric( (pred.bam.st.habseason.hab.density.2.0$fit - 1.96 * pred.bam.st.habseason.hab.density.2.0$se.fit))

high.dens.plotr = plotr [plotr [,which(names(plotr)==itneracting )] == max(plotr [,which(names(plotr)==itneracting) ]),]
low.dens.plotr = plotr [plotr [,names(plotr)==itneracting ] == min(plotr [,names(plotr)==itneracting ]),]

hilo.plotr = rbind(high.dens.plotr, low.dens.plotr)

ylims = c(0 , 0.04)

interested.length = 200

my.data1 = make.dummy.dat.for.model (model.name =bam.st.habseason.hab.density.2.0 , original.df = unsex.trap,
                                     interested.covar = intrest.covar, interested.length = interested.length, 
                                     by.fac = F, interesting.factor.name = NULL)
my.data = my.data1
pred.bam.st.habseason.hab.density.2.0 = predict(bam.st.habseason.hab.density.2.0, newdata= my.data , se=TRUE, type = "link")
my.data$fit = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit))
my.data$upper.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit + 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))
my.data$lower.95 = plogis(as.numeric( pred.bam.st.habseason.hab.density.2.0$fit - 1.96* pred.bam.st.habseason.hab.density.2.0$se.fit))

hilo.plotr$dens.col = hilo.plotr$block.dens
hilo.plotr$dens.col[hilo.plotr$block.dens== max(hilo.plotr$block.dens)] = "red"
hilo.plotr$dens.col[hilo.plotr$block.dens== min(hilo.plotr$block.dens)] = "deepskyblue4"


my.data$expl = my.data[,which(names(my.data)==intrest.covar )]
rugs = unique(unsex.trap[,which(names(unsex.trap)==intrest.covar)])/divider
# both plots on one ----
ggplot()+
    geom_ribbon(aes(ymin = hilo.plotr$lower.95, ymax = hilo.plotr$upper.95, x = hilo.plotr$expl),
                alpha = 0.4, group = hilo.plotr$block.dens, fill = hilo.plotr$dens.col) +
    geom_line(aes(y = high.dens.plotr$fit, x = high.dens.plotr$expl), size = 1, col = "black") +
    geom_line(aes(y = low.dens.plotr$fit, x = low.dens.plotr$expl), size = 1, col = "black") +
    geom_rug( aes(x = rugs), alpha = 0.5 ,sides = "b",# position = "jitter",
              size = 0.1) +
    labs(x = x.lab, y = "Target Sex Capture Probability") +
    scale_y_continuous(limits=c(0,0.015)) +
    scale_x_continuous(limits=c(0,max(plotr$expl))) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18, face = "bold")
    )
# both plots seperate ----
ggplot()+
    geom_ribbon(aes(ymin = high.dens.plotr$lower.95, ymax = high.dens.plotr$upper.95, x = high.dens.plotr$expl),
                alpha = 0.4) +
    geom_line(aes(y = high.dens.plotr$fit, x = high.dens.plotr$expl), size = 1, col = "black") +
    geom_rug( aes(x = rugs), alpha = 0.5 ,sides = "b",# position = "jitter",
              size = 0.1) +
    labs(x = x.lab, y = "Target Sex Capture Probability") +
    scale_y_continuous(limits=ylims) +
    scale_x_continuous(limits=c(0,max(plotr$expl))) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18, face = "bold")
    )

ggplot()+
    geom_ribbon(aes(ymin = low.dens.plotr$lower.95, ymax = low.dens.plotr$upper.95, x = low.dens.plotr$expl),
                alpha = 0.4) +
    geom_line(aes(y = low.dens.plotr$fit, x = low.dens.plotr$expl), size = 1, col = "black") +
    geom_rug( aes(x = rugs), alpha = 0.5 ,sides = "b",# position = "jitter",
              size = 0.1) +
    labs(x = x.lab, y = "Target Sex Capture Probability") +
    #scale_y_continuous(limits=ylims) +
    scale_x_continuous(limits=c(0,max(plotr$expl))) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18, face = "bold")
    )

# both plots on one (link) ----
ggplot()+
    geom_ribbon(aes(ymin = hilo.plotr$link.lower.95, ymax = hilo.plotr$link.upper.95, x = hilo.plotr$expl),
                alpha = 0.4, group = hilo.plotr$block.dens) +
    geom_line(aes(y = high.dens.plotr$link.fit, x = high.dens.plotr$expl), size = 1, col = "black") +
    geom_line(aes(y = low.dens.plotr$link.fit, x = low.dens.plotr$expl), size = 1, col = "black") +
    geom_rug( aes(x = rugs), alpha = 0.5 ,sides = "b",# position = "jitter",
              size = 0.1) +
    labs(x = x.lab, y = "Logit(Target Sex Capture Probability)") +
  #  scale_y_continuous(limits=ylims) +
    scale_x_continuous(limits=c(0,max(plotr$expl))) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = 'black'),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18, face = "bold")
    )
#####

