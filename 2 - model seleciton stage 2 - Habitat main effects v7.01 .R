# habtiat variables includign spline effects
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

## corrplot 0.84 loaded
corr.habs = unsex.trap[, which(names(unsex.trap) %in% c("salat.sum.means", "coastal", "ALT3 ", 
                                                  "water.edge.length", "combo.heather_all.mk2", "combo.dunes_all.mk2", 
                                                  "coast.rock.mk2", "acid.grass.mk2", "woodland.mk2", 
                                                  "combo.beach_saltmarsh.mk2", "bog.mk2", "pasture.mk2", 
                                                  "builtup.mk2", "rough.grass.mk2", "inland.rock.mk2"
))]
names(corr.habs)[names(corr.habs)== "combo.heather_all.mk2"] = "heather"
names(corr.habs)[names(corr.habs)== "combo.dunes_all.mk2"] = "dunes"
names(corr.habs)[names(corr.habs)== "combo.beach_saltmarsh.mk2"] = "littoral sed"
names(corr.habs)[names(corr.habs)== "bog.mk2"] = "bog"
names(corr.habs)[names(corr.habs)== "salat.sum.means"] = "coastal shelter"
names(corr.habs)[names(corr.habs)== "water.edge.length"] = "water edge"
names(corr.habs)[names(corr.habs)== "builtup.mk2"] = "infrastructure"
names(corr.habs)[names(corr.habs)== "woodland.mk2"] = "woodland"
names(corr.habs)[names(corr.habs)== "coast.rock.mk2"] = "coastal rock"
names(corr.habs)[names(corr.habs)== "rough.grass.mk2"] = "rough grass"
names(corr.habs)[names(corr.habs)== "acid.grass.mk2"] = "Acid grass"
names(corr.habs)[names(corr.habs)== "pasture.mk2"] = "Pasture"
names(corr.habs)[names(corr.habs)== "inland.rock.mk2"] = "inland rock"


M <- cor(corr.habs)
round (M, digits = 2)
corrplot(M, method = "circle")



##################################################################################################################################################################
# hab models building up step1 ----

time0 = Sys.time()
bam.s.hab1.0  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                        s(log.dens) + 
                        sub.ses.num +
                        s(sub.sesh.night.num, k = 4) +
                      
                      mink.scent +
                      
                        s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                    data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab1.0 = Sys.time() - time0
save(bam.s.hab1.0, file = paste(sep = "", model.outputs.path, "bam.s.hab1.0"))
rm(bam.s.hab1.0)
gc()

# single habitat effects

time0 = Sys.time()
bam.s.hab1.01  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                       
                       mink.scent +
                       
                         s(salat.sum.means) +
                         coastal +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab1.01 = Sys.time() - time0
save(bam.s.hab1.01, file = paste(sep = "", model.outputs.path, "bam.s.hab1.01"))
rm(bam.s.hab1.01)
gc()

time0 = Sys.time()
bam.s.hab1.02  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                       
                       mink.scent +
                       
                         coastal +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab1.02 = Sys.time() - time0
save(bam.s.hab1.02, file = paste(sep = "", model.outputs.path, "bam.s.hab1.02"))
rm(bam.s.hab1.02)
gc()

time0 = Sys.time()
bam.s.hab1.03  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                       
                       mink.scent +
                       
                         s(ALT3) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab1.03 = Sys.time() - time0
save(bam.s.hab1.03, file = paste(sep = "", model.outputs.path, "bam.s.hab1.03"))
rm(bam.s.hab1.03)
gc()

time0 = Sys.time()
bam.s.hab1.04  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                       
                       mink.scent +
                       
                         s(water.edge.length) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab1.04 = Sys.time() - time0
save(bam.s.hab1.04, file = paste(sep = "", model.outputs.path, "bam.s.hab1.04"))
rm(bam.s.hab1.04)
gc()

time0 = Sys.time()
bam.s.hab1.05  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                       
                       mink.scent +
                       
                         s(combo.heather_all.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab1.05 = Sys.time() - time0
save(bam.s.hab1.05, file = paste(sep = "", model.outputs.path, "bam.s.hab1.05"))
rm(bam.s.hab1.05)
gc()

time0 = Sys.time()
bam.s.hab1.06  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                       
                       mink.scent +
                       
                         s(combo.dunes_all.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab1.06 = Sys.time() - time0
save(bam.s.hab1.06, file = paste(sep = "", model.outputs.path, "bam.s.hab1.06"))
rm(bam.s.hab1.06)
gc()

time0 = Sys.time()
bam.s.hab1.07  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                       
                       mink.scent +
                       
                         s(coast.rock.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab1.07 = Sys.time() - time0
save(bam.s.hab1.07, file = paste(sep = "", model.outputs.path, "bam.s.hab1.07"))
rm(bam.s.hab1.07)
gc()


time0 = Sys.time()
bam.s.hab1.08  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                       
                       mink.scent +
                       
                         s(acid.grass.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab1.08 = Sys.time() - time0
save(bam.s.hab1.08, file = paste(sep = "", model.outputs.path, "bam.s.hab1.08"))
rm(bam.s.hab1.08)
gc()

time0 = Sys.time()
bam.s.hab1.09  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                       
                       mink.scent +
                       
                         s(woodland.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab1.09 = Sys.time() - time0
save(bam.s.hab1.09, file = paste(sep = "", model.outputs.path, "bam.s.hab1.09"))
rm(bam.s.hab1.09)
gc()


time0 = Sys.time()
bam.s.hab1.10  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                       
                       mink.scent +
                       
                         coastal +
                         s(combo.beach_saltmarsh.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab1.10 = Sys.time() - time0
save(bam.s.hab1.10, file = paste(sep = "", model.outputs.path, "bam.s.hab1.10"))
rm(bam.s.hab1.10)
gc()


time0 = Sys.time()
bam.s.hab1.11  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                       
                       mink.scent +
                       
                         s(bog.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab1.11 = Sys.time() - time0
save(bam.s.hab1.11, file = paste(sep = "", model.outputs.path, "bam.s.hab1.11"))
rm(bam.s.hab1.11)
gc()


time0 = Sys.time()
bam.s.hab1.12  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                       
                       mink.scent +
                       
                         s(pasture.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab1.12 = Sys.time() - time0
save(bam.s.hab1.12, file = paste(sep = "", model.outputs.path, "bam.s.hab1.12"))
rm(bam.s.hab1.12)
gc()

time0 = Sys.time()
bam.s.hab1.13  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                       
                       mink.scent +
                       
                         s(builtup.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab1.13 = Sys.time() - time0
save(bam.s.hab1.13, file = paste(sep = "", model.outputs.path, "bam.s.hab1.13"))
rm(bam.s.hab1.13)
gc()

time0 = Sys.time()
bam.s.hab1.14  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                       
                       mink.scent +
                       
                         s(rough.grass.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab1.14 = Sys.time() - time0
save(bam.s.hab1.14, file = paste(sep = "", model.outputs.path, "bam.s.hab1.14"))
rm(bam.s.hab1.14)
gc()

time0 = Sys.time()
bam.s.hab1.15  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                       
                       mink.scent +
                       
                         s(inland.rock.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab1.15 = Sys.time() - time0
save(bam.s.hab1.15, file = paste(sep = "", model.outputs.path, "bam.s.hab1.15"))
rm(bam.s.hab1.15)
gc()


##################################################################################################################################################################
#   aic comparisons ----
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab1.0")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab1.01")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab1.02")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab1.03")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab1.04")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab1.05")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab1.06")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab1.07")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab1.08")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab1.09")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab1.10")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab1.11")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab1.12")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab1.13")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab1.14")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab1.15")

AIC(bam.s.hab1.0, bam.s.hab1.01, 
    bam.s.hab1.02, bam.s.hab1.03, bam.s.hab1.04, bam.s.hab1.05, 
    bam.s.hab1.06, bam.s.hab1.07, bam.s.hab1.08, bam.s.hab1.09, 
    bam.s.hab1.10, bam.s.hab1.11, bam.s.hab1.12, bam.s.hab1.13, bam.s.hab1.14, bam.s.hab1.15)

#               df      AIC
# bam.s.hab1.0  105.6535 15302.78
# bam.s.hab1.01 116.7914 15168.33 **** salat + coastal
# bam.s.hab1.02 125.5508 15285.92
# bam.s.hab1.03 131.1641 15187.17
# bam.s.hab1.04 104.9485 15246.23
# bam.s.hab1.05 105.7536 15226.82
# bam.s.hab1.06 109.4947 15301.09
# bam.s.hab1.07 108.4400 15179.91
# bam.s.hab1.08 105.4085 15288.93
# bam.s.hab1.09 106.7355 15304.79
# bam.s.hab1.10 126.2605 15285.34
# bam.s.hab1.11 107.6495 15289.58
# bam.s.hab1.12 107.3574 15303.70
# bam.s.hab1.13 110.4489 15259.97
# bam.s.hab1.14 109.8202 15299.73
# bam.s.hab1.15 110.6275 15277.04


##################################################################################################################################################################
# lvl2 models ----

time0 = Sys.time()
bam.s.hab2.0  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                        s(log.dens) + 
                        sub.ses.num +
                        s(sub.sesh.night.num, k = 4) +
                        
                        mink.scent +
                        
                        s(salat.sum.means) +
                        coastal +
                        
                        s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                    data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab2.0 = Sys.time() - time0
save(bam.s.hab2.0, file = paste(sep = "", model.outputs.path, "bam.s.hab2.0"))
rm(bam.s.hab2.0)
gc()

time0 = Sys.time()
bam.s.hab2.01  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(ALT3) +

                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab2.01 = Sys.time() - time0
save(bam.s.hab2.01, file = paste(sep = "", model.outputs.path, "bam.s.hab2.01"))
rm(bam.s.hab2.01)
gc()
    
time0 = Sys.time()
bam.s.hab2.02  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(water.edge.length) +

                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")

time.bam.s.hab2.02 = Sys.time() - time0
save(bam.s.hab2.02, file = paste(sep = "", model.outputs.path, "bam.s.hab2.02"))
rm(bam.s.hab2.02)
gc()

time0 = Sys.time()
bam.s.hab2.03  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 

                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")

time.bam.s.hab2.03 = Sys.time() - time0
save(bam.s.hab2.03, file = paste(sep = "", model.outputs.path, "bam.s.hab2.03"))
rm(bam.s.hab2.03)
gc()

time0 = Sys.time()
bam.s.hab2.04  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.dunes_all.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")

time.bam.s.hab2.04 = Sys.time() - time0
save(bam.s.hab2.04, file = paste(sep = "", model.outputs.path, "bam.s.hab2.04"))
rm(bam.s.hab2.04)
gc()

time0 = Sys.time()
bam.s.hab2.05  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(coast.rock.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")

time.bam.s.hab2.05 = Sys.time() - time0
save(bam.s.hab2.05, file = paste(sep = "", model.outputs.path, "bam.s.hab2.05"))
rm(bam.s.hab2.05)
gc()

time0 = Sys.time()
bam.s.hab2.06  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(acid.grass.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")

time.bam.s.hab2.06 = Sys.time() - time0
save(bam.s.hab2.06, file = paste(sep = "", model.outputs.path, "bam.s.hab2.06"))
rm(bam.s.hab2.06)
gc()

time0 = Sys.time()
bam.s.hab2.07  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(woodland.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")

time.bam.s.hab2.07 = Sys.time() - time0
save(bam.s.hab2.07, file = paste(sep = "", model.outputs.path, "bam.s.hab2.07"))
rm(bam.s.hab2.07)
gc()

time0 = Sys.time()
bam.s.hab2.08  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.beach_saltmarsh.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")

time.bam.s.hab2.08 = Sys.time() - time0
save(bam.s.hab2.08, file = paste(sep = "", model.outputs.path, "bam.s.hab2.08"))
rm(bam.s.hab2.08)
gc()

time0 = Sys.time()
bam.s.hab2.09  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(bog.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab2.09 = Sys.time() - time0
save(bam.s.hab2.09, file = paste(sep = "", model.outputs.path, "bam.s.hab2.09"))
rm(bam.s.hab2.09)
gc()

time0 = Sys.time()
bam.s.hab2.10  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(pasture.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab2.10 = Sys.time() - time0
save(bam.s.hab2.10, file = paste(sep = "", model.outputs.path, "bam.s.hab2.10"))
rm(bam.s.hab2.10)
gc()

time0 = Sys.time()
bam.s.hab2.11  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(builtup.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab2.11 = Sys.time() - time0
save(bam.s.hab2.11, file = paste(sep = "", model.outputs.path, "bam.s.hab2.11"))
rm(bam.s.hab2.11)
gc()

time0 = Sys.time()
bam.s.hab2.12  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(rough.grass.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab2.12 = Sys.time() - time0
save(bam.s.hab2.12, file = paste(sep = "", model.outputs.path, "bam.s.hab2.12"))
rm(bam.s.hab2.12)
gc()

time0 = Sys.time()
bam.s.hab2.13  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(inland.rock.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab2.13 = Sys.time() - time0
save(bam.s.hab2.13, file = paste(sep = "", model.outputs.path, "bam.s.hab2.13"))
rm(bam.s.hab2.13)
gc()


##################################################################################################################################################################
# lvl2 AICs ----

load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab2.0")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab2.01")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab2.02")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab2.03")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab2.04")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab2.05")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab2.06")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab2.07")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab2.08")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab2.09")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab2.10")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab2.11")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab2.12")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab2.13")

AIC(bam.s.hab2.0, bam.s.hab2.01, 
  bam.s.hab2.02, bam.s.hab2.03, bam.s.hab2.04, bam.s.hab2.05, bam.s.hab2.06, bam.s.hab2.07, bam.s.hab2.08, bam.s.hab2.09, 
    bam.s.hab2.10, bam.s.hab2.11, bam.s.hab2.12, bam.s.hab2.13)

# df      AIC
# bam.s.hab2.0  116.7914 15168.33
# bam.s.hab2.01 113.2401 15125.28
# bam.s.hab2.02 108.6916 15115.10
# bam.s.hab2.03 110.4179 15106.26 ** heather
# bam.s.hab2.04 109.6742 15155.93
# bam.s.hab2.05 112.6811 15146.81
# bam.s.hab2.06 109.5607 15153.74
# bam.s.hab2.07 109.7560 15154.66
# bam.s.hab2.08 109.1740 15150.58
# bam.s.hab2.09 111.3489 15150.71
# bam.s.hab2.10 110.0809 15154.23
# bam.s.hab2.11 112.9301 15132.20
# bam.s.hab2.12 110.0801 15129.83
# bam.s.hab2.13 113.5236 15152.86

##################################################################################################################################################################
# lvl3 models ----
time0 = Sys.time()
bam.s.hab3.0  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                      s(log.dens) + 
                      sub.ses.num +
                      s(sub.sesh.night.num, k = 4) +
                      
                      mink.scent +
                      
                      s(salat.sum.means) +
                      coastal +
                      s(combo.heather_all.mk2) + 
                      
                      s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                    data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab3.0 = Sys.time() - time0
save(bam.s.hab3.0, file = paste(sep = "", model.outputs.path, "bam.s.hab3.0"))
rm(bam.s.hab3.0)
gc()

time0 = Sys.time()
bam.s.hab3.01  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(ALT3) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab3.01 = Sys.time() - time0
save(bam.s.hab3.01, file = paste(sep = "", model.outputs.path, "bam.s.hab3.01"))
rm(bam.s.hab3.01)
gc()

time0 = Sys.time()
bam.s.hab3.02  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(water.edge.length) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab3.02 = Sys.time() - time0
save(bam.s.hab3.02, file = paste(sep = "", model.outputs.path, "bam.s.hab3.02"))
rm(bam.s.hab3.02)
gc()

time0 = Sys.time()
bam.s.hab3.03  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) +                          
                         s(combo.dunes_all.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab3.03 = Sys.time() - time0
save(bam.s.hab3.03, file = paste(sep = "", model.outputs.path, "bam.s.hab3.03"))
rm(bam.s.hab3.03)
gc()

time0 = Sys.time()
bam.s.hab3.04  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) +
                         s(coast.rock.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab3.04 = Sys.time() - time0
save(bam.s.hab3.04, file = paste(sep = "", model.outputs.path, "bam.s.hab3.04"))
rm(bam.s.hab3.04)
gc()

time0 = Sys.time()
bam.s.hab3.05  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) +
                         s(acid.grass.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab3.05 = Sys.time() - time0
save(bam.s.hab3.05, file = paste(sep = "", model.outputs.path, "bam.s.hab3.05"))
rm(bam.s.hab3.05)
gc()

time0 = Sys.time()
bam.s.hab3.06  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(woodland.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab3.06 = Sys.time() - time0
save(bam.s.hab3.06, file = paste(sep = "", model.outputs.path, "bam.s.hab3.06"))
rm(bam.s.hab3.06)
gc()

time0 = Sys.time()
bam.s.hab3.07  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(combo.beach_saltmarsh.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab3.07 = Sys.time() - time0
save(bam.s.hab3.07, file = paste(sep = "", model.outputs.path, "bam.s.hab3.07"))
rm(bam.s.hab3.07)
gc()

time0 = Sys.time()
bam.s.hab3.08  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(bog.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab3.08 = Sys.time() - time0
save(bam.s.hab3.08, file = paste(sep = "", model.outputs.path, "bam.s.hab3.08"))
rm(bam.s.hab3.08)
gc()

time0 = Sys.time()
bam.s.hab3.09  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(pasture.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab3.09 = Sys.time() - time0
save(bam.s.hab3.09, file = paste(sep = "", model.outputs.path, "bam.s.hab3.09"))
rm(bam.s.hab3.09)
gc()

time0 = Sys.time()
bam.s.hab3.10  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(builtup.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab3.10 = Sys.time() - time0
save(bam.s.hab3.10, file = paste(sep = "", model.outputs.path, "bam.s.hab3.10"))
rm(bam.s.hab3.10)
gc()

time0 = Sys.time()
bam.s.hab3.11  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(rough.grass.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab3.11 = Sys.time() - time0
save(bam.s.hab3.11, file = paste(sep = "", model.outputs.path, "bam.s.hab3.11"))
rm(bam.s.hab3.11)
gc()

time0 = Sys.time()
bam.s.hab3.12  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(inland.rock.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab3.12 = Sys.time() - time0
save(bam.s.hab3.12, file = paste(sep = "", model.outputs.path, "bam.s.hab3.12"))
rm(bam.s.hab3.12)
gc()


##################################################################################################################################################################
# lvl3 AICs ----

load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab3.0")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab3.01")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab3.02")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab3.03")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab3.04")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab3.05")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab3.06")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab3.07")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab3.08")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab3.09")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab3.10")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab3.11")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab3.12")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab3.13")

AIC(bam.s.hab3.0, bam.s.hab3.01, 
    bam.s.hab3.02, bam.s.hab3.03, bam.s.hab3.04, bam.s.hab3.05, bam.s.hab3.06, bam.s.hab3.07, bam.s.hab3.08, bam.s.hab3.09, 
    bam.s.hab3.10, bam.s.hab3.11, bam.s.hab3.12 )

# df      AIC
# bam.s.hab3.0  110.4179 15106.26
# bam.s.hab3.01 114.6952 15079.64 ** ALT
# bam.s.hab3.02 108.7876 15086.48
# bam.s.hab3.03 114.0832 15106.64
# bam.s.hab3.04 111.4528 15103.00
# bam.s.hab3.05 111.5820 15107.39
# bam.s.hab3.06 111.4279 15105.69
# bam.s.hab3.07 110.8627 15096.99
# bam.s.hab3.08 111.3493 15105.43
# bam.s.hab3.09 112.0470 15101.17
# bam.s.hab3.10 115.5981 15087.08
# bam.s.hab3.11 112.2809 15091.63
# bam.s.hab3.12 115.9376 15106.38

##################################################################################################################################################################
# lvl4 models ---- 

time0 = Sys.time()
bam.s.hab4.0  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                        s(log.dens) + 
                        sub.ses.num +
                        s(sub.sesh.night.num, k = 4) +
                        
                        mink.scent +
                        
                        s(salat.sum.means) +
                        coastal +
                        s(combo.heather_all.mk2) + 
                        s(ALT3) +
                        
                        s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                    data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab4.0 = Sys.time() - time0
save(bam.s.hab4.0, file = paste(sep = "", model.outputs.path, "bam.s.hab4.0"))
rm(bam.s.hab4.0)
gc()

time0 = Sys.time()
bam.s.hab4.01  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(ALT3) +
                         s(water.edge.length) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab4.01 = Sys.time() - time0
save(bam.s.hab4.01, file = paste(sep = "", model.outputs.path, "bam.s.hab4.01"))
rm(bam.s.hab4.01)
gc()

time0 = Sys.time()
bam.s.hab4.02  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(ALT3) +
                         s(combo.dunes_all.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab4.02 = Sys.time() - time0
save(bam.s.hab4.02, file = paste(sep = "", model.outputs.path, "bam.s.hab4.02"))
rm(bam.s.hab4.02)
gc()

time0 = Sys.time()
bam.s.hab4.03  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(ALT3) +
                         s(coast.rock.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab4.03 = Sys.time() - time0
save(bam.s.hab4.03, file = paste(sep = "", model.outputs.path, "bam.s.hab4.03"))
rm(bam.s.hab4.03)
gc()

time0 = Sys.time()
bam.s.hab4.04  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(ALT3) +
                         s(acid.grass.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab4.04 = Sys.time() - time0
save(bam.s.hab4.04, file = paste(sep = "", model.outputs.path, "bam.s.hab4.04"))
rm(bam.s.hab4.04)
gc()

time0 = Sys.time()
bam.s.hab4.05  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(ALT3) +
                         s(woodland.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab4.05 = Sys.time() - time0
save(bam.s.hab4.05, file = paste(sep = "", model.outputs.path, "bam.s.hab4.05"))
rm(bam.s.hab4.05)
gc()

time0 = Sys.time()
bam.s.hab4.06  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(ALT3) +
                         s(combo.beach_saltmarsh.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab4.06 = Sys.time() - time0
save(bam.s.hab4.06, file = paste(sep = "", model.outputs.path, "bam.s.hab4.06"))
rm(bam.s.hab4.06)
gc()

time0 = Sys.time()
bam.s.hab4.07  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(ALT3) +
                         s(bog.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab4.07 = Sys.time() - time0
save(bam.s.hab4.07, file = paste(sep = "", model.outputs.path, "bam.s.hab4.07"))
rm(bam.s.hab4.07)
gc()

time0 = Sys.time()
bam.s.hab4.08  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(ALT3) +
                         s(pasture.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab4.08 = Sys.time() - time0
save(bam.s.hab4.08, file = paste(sep = "", model.outputs.path, "bam.s.hab4.08"))
rm(bam.s.hab4.08)
gc()

time0 = Sys.time()
bam.s.hab4.09  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(ALT3) +
                         s(builtup.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab4.09 = Sys.time() - time0
save(bam.s.hab4.09, file = paste(sep = "", model.outputs.path, "bam.s.hab4.09"))
rm(bam.s.hab4.09)
gc()

time0 = Sys.time()
bam.s.hab4.10  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(ALT3) +
                         s(rough.grass.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab4.10 = Sys.time() - time0
save(bam.s.hab4.10, file = paste(sep = "", model.outputs.path, "bam.s.hab4.10"))
rm(bam.s.hab4.10)
gc()

time0 = Sys.time()
bam.s.hab4.11  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(ALT3) +
                         s(inland.rock.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab4.11 = Sys.time() - time0
save(bam.s.hab4.11, file = paste(sep = "", model.outputs.path, "bam.s.hab4.11"))
rm(bam.s.hab4.11)
gc()

##################################################################################################################################################################
# lvl4 AICs ----

load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab4.0")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab4.01")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab4.02")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab4.03")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab4.04")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab4.05")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab4.06")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab4.07")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab4.08")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab4.09")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab4.10")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab4.11")

AIC(bam.s.hab4.0, bam.s.hab4.01, bam.s.hab4.02, bam.s.hab4.03, bam.s.hab4.04, bam.s.hab4.05, bam.s.hab4.06, bam.s.hab4.07, bam.s.hab4.08, bam.s.hab4.09, 
    bam.s.hab4.10, bam.s.hab4.11)

# df      AIC
# bam.s.hab4.0  114.6952 15079.64
# bam.s.hab4.01 113.1659 15067.81
# bam.s.hab4.02 117.4896 15074.39
# bam.s.hab4.03 115.6303 15080.18
# bam.s.hab4.04 116.0136 15081.08
# bam.s.hab4.05 115.5819 15078.73
# bam.s.hab4.06 114.7254 15065.55
# bam.s.hab4.07 115.9921 15081.19
# bam.s.hab4.08 116.2109 15071.17
# bam.s.hab4.09 120.0132 15062.39 **
# bam.s.hab4.10 124.3340 15078.24
# bam.s.hab4.11 120.2444 15079.77

##################################################################################################################################################################
# lvl5 models ---- 

time0 = Sys.time()
bam.s.hab5.0  =  bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(ALT3) +
                         s(builtup.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab5.0 = Sys.time() - time0
save(bam.s.hab5.0, file = paste(sep = "", model.outputs.path, "bam.s.hab5.0"))
rm(bam.s.hab5.0)
gc()

time0 = Sys.time()
bam.s.hab5.01  =  bam(unsex.capt ~ s(julian, bs = "cc")  +
                          s(log.dens) + 
                          sub.ses.num +
                          s(sub.sesh.night.num, k = 4) +
                          
                          mink.scent +
                          
                          s(salat.sum.means) +
                          coastal +
                          s(combo.heather_all.mk2) + 
                          s(ALT3) +
                          s(builtup.mk2) + 
                          s(water.edge.length) +
                          
                          s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                      data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab5.01 = Sys.time() - time0
save(bam.s.hab5.01, file = paste(sep = "", model.outputs.path, "bam.s.hab5.01"))
rm(bam.s.hab5.01)
gc()

time0 = Sys.time()
bam.s.hab5.02  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(ALT3) +
                         s(builtup.mk2) + 
                         s(combo.dunes_all.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab5.02 = Sys.time() - time0
save(bam.s.hab5.02, file = paste(sep = "", model.outputs.path, "bam.s.hab5.02"))
rm(bam.s.hab5.02)
gc()

time0 = Sys.time()
bam.s.hab5.03  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(ALT3) +
                         s(builtup.mk2) + 
                         s(coast.rock.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab5.03 = Sys.time() - time0
save(bam.s.hab5.03, file = paste(sep = "", model.outputs.path, "bam.s.hab5.03"))
rm(bam.s.hab5.03)
gc()

time0 = Sys.time()
bam.s.hab5.04  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(ALT3) +
                         s(builtup.mk2) + 
                         s(acid.grass.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab5.04 = Sys.time() - time0
save(bam.s.hab5.04, file = paste(sep = "", model.outputs.path, "bam.s.hab5.04"))
rm(bam.s.hab5.04)
gc()

time0 = Sys.time()
bam.s.hab5.05  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(ALT3) +
                         s(builtup.mk2) + 
                         s(woodland.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab5.05 = Sys.time() - time0
save(bam.s.hab5.05, file = paste(sep = "", model.outputs.path, "bam.s.hab5.05"))
rm(bam.s.hab5.05)
gc()

time0 = Sys.time()
bam.s.hab5.06  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(ALT3) +
                         s(builtup.mk2) + 
                         s(combo.beach_saltmarsh.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab5.06 = Sys.time() - time0
save(bam.s.hab5.06, file = paste(sep = "", model.outputs.path, "bam.s.hab5.06"))
rm(bam.s.hab5.06)
gc()

time0 = Sys.time()
bam.s.hab5.07  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(ALT3) +
                         s(builtup.mk2) + 
                         s(bog.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab5.07 = Sys.time() - time0
save(bam.s.hab5.07, file = paste(sep = "", model.outputs.path, "bam.s.hab5.07"))
rm(bam.s.hab5.07)
gc()

time0 = Sys.time()
bam.s.hab5.08  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(ALT3) +
                         s(builtup.mk2) + 
                         s(pasture.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab5.08 = Sys.time() - time0
save(bam.s.hab5.08, file = paste(sep = "", model.outputs.path, "bam.s.hab5.08"))
rm(bam.s.hab5.08)
gc()

time0 = Sys.time()
bam.s.hab5.09  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab5.09 = Sys.time() - time0
save(bam.s.hab5.09, file = paste(sep = "", model.outputs.path, "bam.s.hab5.09"))
rm(bam.s.hab5.09)
gc()

time0 = Sys.time()
bam.s.hab5.10  = bam(unsex.capt ~ s(julian, bs = "cc")  +
                         s(log.dens) + 
                         sub.ses.num +
                         s(sub.sesh.night.num, k = 4) +
                         
                         mink.scent +
                         
                         s(salat.sum.means) +
                         coastal +
                         s(combo.heather_all.mk2) + 
                         s(ALT3) +
                         s(builtup.mk2) + 
                         s(inland.rock.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab5.10 = Sys.time() - time0
save(bam.s.hab5.10, file = paste(sep = "", model.outputs.path, "bam.s.hab5.10"))
rm(bam.s.hab5.10)
gc()



##################################################################################################################################################################
# lvl5 AICs ----

load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab5.0")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab5.01")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab5.02")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab5.03")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab5.04")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab5.05")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab5.06")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab5.07")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab5.08")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab5.09")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab5.10")

AIC(bam.s.hab5.0, bam.s.hab5.01, bam.s.hab5.02, bam.s.hab5.03, bam.s.hab5.04, bam.s.hab5.05, 
    bam.s.hab5.06, bam.s.hab5.07, bam.s.hab5.08 , bam.s.hab5.09, 
    bam.s.hab5.10 
    )

# df       AIC
# bam.s.hab5.0  120.01320 15062.394
# bam.s.hab5.01 118.7802  15050.44
# bam.s.hab5.02 122.63205 15057.615
# bam.s.hab5.03 121.07068 15065.306
# bam.s.hab5.04 121.39886 15063.647
# bam.s.hab5.05 120.92807 15060.235
# bam.s.hab5.06 120.27703 15046.275
# bam.s.hab5.07 121.27463 15064.828
# bam.s.hab5.08 121.80037 15040.712
# bam.s.hab5.09 119.06435 15027.857 ** rough grass
# bam.s.hab5.10 125.66428 15062.671

##################################################################################################################################################################
# lvl6 models ---- 

time0 = Sys.time()
bam.s.hab6.0  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                        
                        s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                    data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab6.0 = Sys.time() - time0
save(bam.s.hab6.0, file = paste(sep = "", model.outputs.path, "bam.s.hab6.0"))
rm(bam.s.hab6.0)
gc()

time0 = Sys.time()
bam.s.hab6.01  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab6.01 = Sys.time() - time0
save(bam.s.hab6.01, file = paste(sep = "", model.outputs.path, "bam.s.hab6.01"))
rm(bam.s.hab6.01)
gc()

time0 = Sys.time()
bam.s.hab6.02  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(combo.dunes_all.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab6.02 = Sys.time() - time0
save(bam.s.hab6.02, file = paste(sep = "", model.outputs.path, "bam.s.hab6.02"))
rm(bam.s.hab6.02)
gc()

time0 = Sys.time()
bam.s.hab6.03  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(coast.rock.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab6.03 = Sys.time() - time0
save(bam.s.hab6.03, file = paste(sep = "", model.outputs.path, "bam.s.hab6.03"))
rm(bam.s.hab6.03)
gc()

time0 = Sys.time()
bam.s.hab6.04  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(acid.grass.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab6.04 = Sys.time() - time0
save(bam.s.hab6.04, file = paste(sep = "", model.outputs.path, "bam.s.hab6.04"))
rm(bam.s.hab6.04)
gc()

time0 = Sys.time()
bam.s.hab6.05  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(woodland.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab6.05 = Sys.time() - time0
save(bam.s.hab6.05, file = paste(sep = "", model.outputs.path, "bam.s.hab6.05"))
rm(bam.s.hab6.05)
gc()

time0 = Sys.time()
bam.s.hab6.06  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(combo.beach_saltmarsh.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab6.06 = Sys.time() - time0
save(bam.s.hab6.06, file = paste(sep = "", model.outputs.path, "bam.s.hab6.06"))
rm(bam.s.hab6.06)
gc()

time0 = Sys.time()
bam.s.hab6.07  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(bog.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab6.07 = Sys.time() - time0
save(bam.s.hab6.07, file = paste(sep = "", model.outputs.path, "bam.s.hab6.07"))
rm(bam.s.hab6.07)
gc()

time0 = Sys.time()
bam.s.hab6.08  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(pasture.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab6.08 = Sys.time() - time0
save(bam.s.hab6.08, file = paste(sep = "", model.outputs.path, "bam.s.hab6.08"))
rm(bam.s.hab6.08)
gc()

time0 = Sys.time()
bam.s.hab6.09  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(inland.rock.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab6.09 = Sys.time() - time0
save(bam.s.hab6.09, file = paste(sep = "", model.outputs.path, "bam.s.hab6.09"))
rm(bam.s.hab6.09)
gc()


##################################################################################################################################################################
# lvl6 AICs ----

load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab6.0")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab6.01")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab6.02")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab6.03")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab6.04")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab6.05")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab6.06")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab6.07")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab6.08")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab6.09")

AIC(bam.s.hab6.0, bam.s.hab6.01, bam.s.hab6.02, bam.s.hab6.03, bam.s.hab6.04, 
    bam.s.hab6.05, bam.s.hab6.06, bam.s.hab6.07, 
    bam.s.hab6.08, bam.s.hab6.09 
    )

# df      AIC
# bam.s.hab6.0  119.0644 15027.86
# bam.s.hab6.01 118.9679 15015.95 **water.edge.length
# bam.s.hab6.02 119.2234 15028.65
# bam.s.hab6.03 119.2665 15026.40
# bam.s.hab6.04 122.1525 15029.29
# bam.s.hab6.05 119.4119 15025.65
# bam.s.hab6.06 124.6718 15028.66
# bam.s.hab6.07 121.5792 15017.89
# bam.s.hab6.08 126.8720 15028.42
# bam.s.hab6.09 118.8378 15027.95
##################################################################################################################################################################
# LVL7 models ---- 
time0 = Sys.time()
bam.s.hab7.0  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                        
                        s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                    data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab7.0 = Sys.time() - time0
save(bam.s.hab7.0, file = paste(sep = "", model.outputs.path, "bam.s.hab7.0"))
rm(bam.s.hab7.0)
gc()

time0 = Sys.time()
bam.s.hab7.01  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(combo.dunes_all.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab7.01 = Sys.time() - time0
save(bam.s.hab7.01, file = paste(sep = "", model.outputs.path, "bam.s.hab7.01"))
rm(bam.s.hab7.01)
gc()

time0 = Sys.time()
bam.s.hab7.02  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(coast.rock.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab7.02 = Sys.time() - time0
save(bam.s.hab7.02, file = paste(sep = "", model.outputs.path, "bam.s.hab7.02"))
rm(bam.s.hab7.02)
gc()

time0 = Sys.time()
bam.s.hab7.03  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(acid.grass.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab7.03 = Sys.time() - time0
save(bam.s.hab7.03, file = paste(sep = "", model.outputs.path, "bam.s.hab7.03"))
rm(bam.s.hab7.03)
gc()

time0 = Sys.time()
bam.s.hab7.04  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(woodland.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab7.04 = Sys.time() - time0
save(bam.s.hab7.04, file = paste(sep = "", model.outputs.path, "bam.s.hab7.04"))
rm(bam.s.hab7.04)
gc()

time0 = Sys.time()
bam.s.hab7.05  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab7.05 = Sys.time() - time0
save(bam.s.hab7.05, file = paste(sep = "", model.outputs.path, "bam.s.hab7.05"))
rm(bam.s.hab7.05)
gc()

time0 = Sys.time()
bam.s.hab7.06  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(bog.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab7.06 = Sys.time() - time0
save(bam.s.hab7.06, file = paste(sep = "", model.outputs.path, "bam.s.hab7.06"))
rm(bam.s.hab7.06)
gc()

time0 = Sys.time()
bam.s.hab7.07  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(pasture.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab7.07 = Sys.time() - time0
save(bam.s.hab7.07, file = paste(sep = "", model.outputs.path, "bam.s.hab7.07"))
rm(bam.s.hab7.07)
gc()

time0 = Sys.time()
bam.s.hab7.08  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(inland.rock.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab7.08 = Sys.time() - time0
save(bam.s.hab7.08, file = paste(sep = "", model.outputs.path, "bam.s.hab7.08"))
rm(bam.s.hab7.08)
gc()

##################################################################################################################################################################
# lvl7 AICs ----

load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab7.0")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab7.01")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab7.02")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab7.03")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab7.04")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab7.05")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab7.06")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab7.07")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab7.08")

AIC(bam.s.hab7.0, bam.s.hab7.01, bam.s.hab7.02, bam.s.hab7.03, bam.s.hab7.04, 
    bam.s.hab7.05, bam.s.hab7.06, bam.s.hab7.07, bam.s.hab7.08 
     )

# df      AIC
# bam.s.hab7.0  118.9679 15015.95
# bam.s.hab7.01 120.1042 15018.34
# bam.s.hab7.02 119.9288 15017.03
# bam.s.hab7.03 123.6490 15021.74
# bam.s.hab7.04 120.2235 15014.58
# bam.s.hab7.05 118.9462 15007.32 **   sediment
# bam.s.hab7.06 119.8924 15014.35
# bam.s.hab7.07 120.5386 15007.69 * pasture
# bam.s.hab7.08 128.2826 15026.41

##################################################################################################################################################################
# LVL8 models ----

time0 = Sys.time()
bam.s.hab8.0  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                        
                        s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                    data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab8.0 = Sys.time() - time0
save(bam.s.hab8.0, file = paste(sep = "", model.outputs.path, "bam.s.hab8.0"))
rm(bam.s.hab8.0)
gc()

time0 = Sys.time()
bam.s.hab8.01  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(combo.dunes_all.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab8.01 = Sys.time() - time0
save(bam.s.hab8.01, file = paste(sep = "", model.outputs.path, "bam.s.hab8.01"))
rm(bam.s.hab8.01)
gc()

time0 = Sys.time()
bam.s.hab8.02  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(coast.rock.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab8.02 = Sys.time() - time0
save(bam.s.hab8.02, file = paste(sep = "", model.outputs.path, "bam.s.hab8.02"))
rm(bam.s.hab8.02)
gc()

time0 = Sys.time()
bam.s.hab8.03  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(acid.grass.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab8.03 = Sys.time() - time0
save(bam.s.hab8.03, file = paste(sep = "", model.outputs.path, "bam.s.hab8.03"))
rm(bam.s.hab8.03)
gc()

time0 = Sys.time()
bam.s.hab8.04  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(woodland.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab8.04 = Sys.time() - time0
save(bam.s.hab8.04, file = paste(sep = "", model.outputs.path, "bam.s.hab8.04"))
rm(bam.s.hab8.04)
gc()

time0 = Sys.time()
bam.s.hab8.05  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(bog.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab8.05 = Sys.time() - time0
save(bam.s.hab8.05, file = paste(sep = "", model.outputs.path, "bam.s.hab8.05"))
rm(bam.s.hab8.05)
gc()

time0 = Sys.time()
bam.s.hab8.06  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab8.06 = Sys.time() - time0
save(bam.s.hab8.06, file = paste(sep = "", model.outputs.path, "bam.s.hab8.06"))
rm(bam.s.hab8.06)
gc()

time0 = Sys.time()
bam.s.hab8.07  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(inland.rock.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab8.07 = Sys.time() - time0
save(bam.s.hab8.07, file = paste(sep = "", model.outputs.path, "bam.s.hab8.07"))
rm(bam.s.hab8.07)
gc()
 

##################################################################################################################################################################
# lvl8 AICs ----

load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab8.0")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab8.01")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab8.02")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab8.03")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab8.04")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab8.05")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab8.06")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab8.07")

AIC(bam.s.hab8.0, bam.s.hab8.01, 
    bam.s.hab8.02, bam.s.hab8.03, bam.s.hab8.04, 
    bam.s.hab8.05, bam.s.hab8.06, bam.s.hab8.07 
)

# df      AIC
# bam.s.hab8.0  118.9462 15007.32
# bam.s.hab8.01 120.0791 15008.27
# bam.s.hab8.02 126.8556 15013.88
# bam.s.hab8.03 123.4485 15012.64
# bam.s.hab8.04 120.1098 15005.59
# bam.s.hab8.05 130.1739 15023.18
# bam.s.hab8.06 120.0997 15000.09 ** pasture
# bam.s.hab8.07 119.6046 15009.16

##################################################################################################################################################################
# LVL9 models ----

time0 = Sys.time()
bam.s.hab9.0  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                        
                        s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                    data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab9.0 = Sys.time() - time0
save(bam.s.hab9.0, file = paste(sep = "", model.outputs.path, "bam.s.hab9.0"))
rm(bam.s.hab9.0)
gc()

time0 = Sys.time()
bam.s.hab9.01  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(combo.dunes_all.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab9.01 = Sys.time() - time0
save(bam.s.hab9.01, file = paste(sep = "", model.outputs.path, "bam.s.hab9.01"))
rm(bam.s.hab9.01)
gc()

time0 = Sys.time()
bam.s.hab9.02  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(coast.rock.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab9.02 = Sys.time() - time0
save(bam.s.hab9.02, file = paste(sep = "", model.outputs.path, "bam.s.hab9.02"))
rm(bam.s.hab9.02)
gc()

time0 = Sys.time()
bam.s.hab9.03  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(acid.grass.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab9.03 = Sys.time() - time0
save(bam.s.hab9.03, file = paste(sep = "", model.outputs.path, "bam.s.hab9.03"))
rm(bam.s.hab9.03)
gc()

time0 = Sys.time()
bam.s.hab9.04  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(woodland.mk2) + 
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab9.04 = Sys.time() - time0
save(bam.s.hab9.04, file = paste(sep = "", model.outputs.path, "bam.s.hab9.04"))
rm(bam.s.hab9.04)
gc()

time0 = Sys.time()
bam.s.hab9.05  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
time.bam.s.hab9.05 = Sys.time() - time0
save(bam.s.hab9.05, file = paste(sep = "", model.outputs.path, "bam.s.hab9.05"))
rm(bam.s.hab9.05)
gc()

time0 = Sys.time()
bam.s.hab9.06  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                         s(inland.rock.mk2) +
                         
                         s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                     data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab9.06 = Sys.time() - time0
save(bam.s.hab9.06, file = paste(sep = "", model.outputs.path, "bam.s.hab9.06"))
rm(bam.s.hab9.06)
gc()


##################################################################################################################################################################
# lvl9 AICs ----

load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab9.0")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab9.01")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab9.02")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab9.03")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab9.04")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab9.05")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab9.06")

AIC(bam.s.hab9.0, 
    bam.s.hab9.01, 
    bam.s.hab9.02, bam.s.hab9.03, bam.s.hab9.04, 
    bam.s.hab9.05, bam.s.hab9.06 
)

# df      AIC
# bam.s.hab9.0  120.0997 15000.09
# bam.s.hab9.01 121.1921 15000.05
# bam.s.hab9.02 126.0541 15000.64
# bam.s.hab9.03 121.3230 15002.15
# bam.s.hab9.04 121.2400 14999.05
# bam.s.hab9.05 123.6469 14996.36 ** bog
# bam.s.hab9.06 121.0639 15002.50

##################################################################################################################################################################
# LVL10 models ----

time0 = Sys.time()
bam.s.hab10.0  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
time.bam.s.hab10.0 = Sys.time() - time0
save(bam.s.hab10.0, file = paste(sep = "", model.outputs.path, "bam.s.hab10.0"))
rm(bam.s.hab10.0)
gc()

time0 = Sys.time()
bam.s.hab10.01  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                          s(combo.dunes_all.mk2) + 
                          
                          s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                      data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab10.01 = Sys.time() - time0
save(bam.s.hab10.01, file = paste(sep = "", model.outputs.path, "bam.s.hab10.01"))
rm(bam.s.hab10.01)
gc()

time0 = Sys.time()
bam.s.hab10.02  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                          s(coast.rock.mk2) +
                          
                          s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                      data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab10.02 = Sys.time() - time0
save(bam.s.hab10.02, file = paste(sep = "", model.outputs.path, "bam.s.hab10.02"))
rm(bam.s.hab10.02)
gc()

time0 = Sys.time()
bam.s.hab10.03  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                          s(acid.grass.mk2) + 
                          
                          s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                      data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab10.03 = Sys.time() - time0
save(bam.s.hab10.03, file = paste(sep = "", model.outputs.path, "bam.s.hab10.03"))
rm(bam.s.hab10.03)
gc()

time0 = Sys.time()
bam.s.hab10.04  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                          s(woodland.mk2) + 
                          
                          s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                      data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab10.04 = Sys.time() - time0
save(bam.s.hab10.04, file = paste(sep = "", model.outputs.path, "bam.s.hab10.04"))
rm(bam.s.hab10.04)
gc()

time0 = Sys.time()
bam.s.hab10.05  = bam(unsex.capt ~ s(julian, bs = "cc")  +
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
                          s(inland.rock.mk2) + 
                          
                          s(two.month.mega, bs = "re") + s(fx.trapper, bs = "re"),
                      data = unsex.trap , family = binomial, method = "ML")
time.bam.s.hab10.05 = Sys.time() - time0
save(bam.s.hab10.05, file = paste(sep = "", model.outputs.path, "bam.s.hab10.05"))
rm(bam.s.hab10.05)
gc() 

##################################################################################################################################################################
# lvl10 AICs ----

load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab10.0")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab10.01")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab10.02")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab10.03")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab10.04")
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab10.05")

AIC(bam.s.hab10.0, 
    bam.s.hab10.01, bam.s.hab10.02, bam.s.hab10.03, bam.s.hab10.04, 
    bam.s.hab10.05)
# df      AIC
# bam.s.hab10.0  123.6469 14996.36 **
# bam.s.hab10.01 124.9851 14997.17 -  dunes
# bam.s.hab10.02 129.7251 14998.47 -  Coastal rock
# bam.s.hab10.03 128.6051 15002.73 - acid grass
# bam.s.hab10.04 124.6187 14995.12 - woodland
# bam.s.hab10.05 124.8103 14999.49 - Inland rock 

##################################################################################################################################################################
# investigate best model ----
load("T:\\HMP\\HMP work\\Analysis\\objects.corrections\\bam.s.hab10.0")

summary (bam.s.hab10.0)
# Family: binomial 
# Link function: logit 
# 
# Formula:
#     sex.cap ~ s(julian, bs = "cc") + female + s(log.dens) + s(day.since.start) + 
#     sub.ses.num + s(sub.sesh.night.num, k = 4) + s(depletion.sex) + 
#     s(traps.open.within.buffer) + tern.col + boat.access + othersex.scent + 
#     samesex.scent + s(coast.dist) + s(salat.sum.means) + s(combo.heather_all.mk2) + 
#     s(builtup.mk2) + s(rough.grass.mk2) + s(pasture.mk2) + s(woodland.mk2) + 
#     s(beach.mk2) + s(bog.mk2) + s(two.month.mega, bs = "re") + 
#     s(fx.trapper, bs = "re")
# 
# Parametric coefficients:
#     Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -5.59013    0.13296 -42.043  < 2e-16 ***
#     female          0.09509    0.05256   1.809   0.0704 .  
# sub.ses.num    -0.58369    0.08580  -6.803 1.03e-11 ***
#     tern.col        0.65254    0.32425   2.012   0.0442 *  
#     boat.access     0.07277    0.08966   0.812   0.4170    
# othersex.scent  1.04876    0.10589   9.904  < 2e-16 ***
#     samesex.scent   0.62991    0.11832   5.324 1.02e-07 ***
#     ---
#     Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Approximate significance of smooth terms:
#     edf  Ref.df  Chi.sq  p-value    
# s(julian)                    3.630   8.000 576.800 4.15e-08 ***
#     s(log.dens)                  1.017   1.028   3.593  0.06045 .  
# s(day.since.start)           1.018   1.024  34.976 3.35e-09 ***
#     s(sub.sesh.night.num)        1.001   1.001  76.687  < 2e-16 ***
#     s(depletion.sex)             4.760   5.681  58.238 1.80e-10 ***
#     s(traps.open.within.buffer)  3.482   4.293  74.263 8.73e-15 ***

#     s(coast.dist)                5.982   7.135  35.730 1.57e-05 ***
#     s(salat.sum.means)           2.877   3.678  40.802 1.34e-07 ***
#     s(combo.heather_all.mk2)     1.013   1.025  29.645 6.42e-08 ***
#     s(builtup.mk2)               5.427   6.525  40.447 5.86e-07 ***
#     s(rough.grass.mk2)           1.014   1.028  16.064 6.26e-05 ***
#     s(pasture.mk2)               1.013   1.026   8.170  0.00450 ** 
#     s(woodland.mk2)              1.017   1.034   6.401  0.01234 *  
#     s(beach.mk2)                 1.012   1.024   8.575  0.00361 ** 
#     s(bog.mk2)                   1.022   1.043   5.036  0.02844 *  

#     s(two.month.mega)           61.901 119.000 219.805  < 2e-16 ***
#     s(fx.trapper)               12.383  22.000  66.999 1.74e-08 ***
#     ---
#     Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# R-sq.(adj) =  0.014   Deviance explained = 52.1%
# -ML = 4.816e+05  Scale est. = 1         n = 350922

vcov(bam.s.hab10.0)

par(mfrow = c(3,2))

interest.covar = "salat.sum.means"
interested.factor = "boat.access"
my.data = make.dummy.dat.for.model(model.name = bam.s.hab10.0 , original.df = unsex.trap, interested.covar = interest.covar , 
                                   interested.length = 100, by.fac = F, interesting.factor.name = interested.factor)
pred.bam.s.hab10.0 = predict(bam.s.hab10.0, newdata= my.data , se=TRUE, type = "link")
plot.gam.func3( model.used = bam.s.hab10.0 , expl.name = interest.covar , 
                prediction.object = pred.bam.s.hab10.0, dummy.prediciton.data = my.data, 
                orignal.data = unsex.trap , num.means = 2, mean.range = 0, main = "seaweed" , ylimer = c(0,0.01) )

interest.covar = "log.combo.heather_all.mk2"
interested.factor = "boat.access"
my.data = make.dummy.dat.for.model(model.name = bam.s.hab10.0 , original.df = unsex.trap, interested.covar = interest.covar , 
                                   interested.length = 100, by.fac = F, interesting.factor.name = interested.factor)
pred.bam.s.hab10.0 = predict(bam.s.hab10.0, newdata= my.data , se=TRUE, type = "link")
plot.gam.func3( model.used = bam.s.hab10.0 , expl.name = interest.covar , 
                prediction.object = pred.bam.s.hab10.0, dummy.prediciton.data = my.data, 
                orignal.data = unsex.trap , num.means = 2, mean.range = 0, main = "heather" , ylimer = c(0,0.01) )

interest.covar = "log.coast.dist"
interested.factor = "boat.access"
my.data = make.dummy.dat.for.model(model.name = bam.s.hab10.0 , original.df = unsex.trap, interested.covar = interest.covar , 
                                   interested.length = 100, by.fac = F, interesting.factor.name = interested.factor)
pred.bam.s.hab10.0 = predict(bam.s.hab10.0, newdata= my.data , se=TRUE, type = "link")
plot.gam.func3( model.used = bam.s.hab10.0 , expl.name = interest.covar , 
                prediction.object = pred.bam.s.hab10.0, dummy.prediciton.data = my.data, 
                orignal.data = unsex.trap , num.means = 2, mean.range = 0, main = "Coast dist" , ylimer = c(0,0.01) )

interest.covar = "log.ALT3"
interested.factor = "boat.access"
my.data = make.dummy.dat.for.model(model.name = bam.s.hab10.0 , original.df = unsex.trap, interested.covar = interest.covar , 
                                   interested.length = 100, by.fac = F, interesting.factor.name = interested.factor)
pred.bam.s.hab10.0 = predict(bam.s.hab10.0, newdata= my.data , se=TRUE, type = "link")
plot.gam.func3( model.used = bam.s.hab10.0 , expl.name = interest.covar , 
                prediction.object = pred.bam.s.hab10.0, dummy.prediciton.data = my.data, 
                orignal.data = unsex.trap , num.means = 2, mean.range = 0 , main = "ALT", ylimer = c(0,0.01) )

interest.covar = "log.rough.grass.mk2"
interested.factor = "boat.access"
my.data = make.dummy.dat.for.model(model.name = bam.s.hab10.0 , original.df = unsex.trap, interested.covar = interest.covar , 
                                   interested.length = 100, by.fac = F, interesting.factor.name = interested.factor)
pred.bam.s.hab10.0 = predict(bam.s.hab10.0, newdata= my.data , se=TRUE, type = "link")
my.data$log.combo.heather_all.mk2 = my.data$log.combo.heather_all.mk2/1000000
plot.gam.func3( model.used = bam.s.hab10.0 , expl.name = interest.covar , 
                prediction.object = pred.bam.s.hab10.0, dummy.prediciton.data = my.data, 
                orignal.data = unsex.trap , num.means = 2, mean.range = 0, main = "Rough Grass" , ylimer = c(0,0.01) )

interest.covar = "log.builtup.mk2"
interested.factor = "boat.access"
my.data = make.dummy.dat.for.model(model.name = bam.s.hab10.0 , original.df = unsex.trap, interested.covar = interest.covar , 
                                   interested.length = 100, by.fac = F, interesting.factor.name = interested.factor)
pred.bam.s.hab10.0 = predict(bam.s.hab10.0, newdata= my.data , se=TRUE, type = "link")
my.data$log.rough.grass.mk2 = my.data$log.rough.grass.mk2/1000000
plot.gam.func3( model.used = bam.s.hab10.0 , expl.name = interest.covar , 
                prediction.object = pred.bam.s.hab10.0, dummy.prediciton.data = my.data, 
                orignal.data = unsex.trap , num.means = 2, mean.range = 0 , main = "built up", ylimer = c(0,0.01) )

interest.covar = "log.combo.bog_all.mk2"
interested.factor = "boat.access"
my.data = make.dummy.dat.for.model(model.name = bam.s.hab10.0 , original.df = unsex.trap, interested.covar = interest.covar , 
                                   interested.length = 100, by.fac = F, interesting.factor.name = interested.factor)
pred.bam.s.hab10.0 = predict(bam.s.hab10.0, newdata= my.data , se=TRUE, type = "link")
my.data$log.pasture.mk2 = my.data$log.pasture.mk2/1000000
plot.gam.func3( model.used = bam.s.hab10.0 , expl.name = interest.covar , 
                prediction.object = pred.bam.s.hab10.0, dummy.prediciton.data = my.data, 
                orignal.data = unsex.trap , num.means = 2, mean.range = 0 , main = "bog", ylimer = c(0,0.01) )

interest.covar = "log.pasture.mk2"
interested.factor = "boat.access"
my.data = make.dummy.dat.for.model(model.name = bam.s.hab10.0 , original.df = unsex.trap, interested.covar = interest.covar , 
                                   interested.length = 100, by.fac = F, interesting.factor.name = interested.factor)
pred.bam.s.hab10.0 = predict(bam.s.hab10.0, newdata= my.data , se=TRUE, type = "link")
my.data$log.pasture.mk2 = my.data$log.pasture.mk2/1000000
plot.gam.func3( model.used = bam.s.hab10.0 , expl.name = interest.covar , 
                prediction.object = pred.bam.s.hab10.0, dummy.prediciton.data = my.data, 
                orignal.data = unsex.trap , num.means = 2, mean.range = 0 , main = "Pasture", ylimer = c(0,0.01) )

interest.covar = "log.woodland.mk2"
interested.factor = "boat.access"
my.data = make.dummy.dat.for.model(model.name = bam.s.hab10.0 , original.df = unsex.trap, interested.covar = interest.covar , 
                                   interested.length = 100, by.fac = F, interesting.factor.name = interested.factor)
pred.bam.s.hab10.0 = predict(bam.s.hab10.0, newdata= my.data , se=TRUE, type = "link")
my.data$log.pasture.mk2 = my.data$log.pasture.mk2/1000000
plot.gam.func3( model.used = bam.s.hab10.0 , expl.name = interest.covar , 
                prediction.object = pred.bam.s.hab10.0, dummy.prediciton.data = my.data, 
                orignal.data = unsex.trap , num.means = 2, mean.range = 0, , main = "Woodland", ylimer = c(0,0.01) )

