##################################################################
#
# pev.frailty
# ===========
#
# R-function for computation of Proportion of Explained Variation
# (PEV) for models with log-normal frailty term
#
# Gleiss, A. & Schemper, M. Explained variation in shared frailty 
# models,
# Statistics in Medicine 2018, 37:1482 - 1490
#
# Author:  Andreas Gleiss
# Version: 1.5 (new output)
# Date:    15 May 2024
#
# Arguments:
# ==========
# 
# ff       formula for fixed effects
#          (e.g. formula(Surv(y,uncens) ~ trt))
# frailty  variable to be used as random effect (e.g. eortc$center)
#          (assumed as consecutive integers starting with 1)
# rterms   terms for random effects 
#          (e.g. "1+trt"; default: "1" for random intercept only)
# dat      name of data frame containg model terms (e.g. eortc)
#
##################################################################

#setwd("E:/Andreas/Forschung/PROGRAMS/PEV_frailty")
source("./fsurev3.R") # install from github.com/agleiss

library(coxme)
library(rms)
library(survival)


pev.frailty <- function(ff,frailty,dat,rterms="1") {
 
  frailty_<-frailty
  tt<-terms(ff)
  ff.all<-reformulate(c(attr(tt,"term.labels"),paste("(",rterms,"|factor(frailty_))")),ff[[2]])
  ff.all2<-reformulate(c(attr(tt,"term.labels"),"offset(frailty.est2_)"),ff[[2]])
  ff.rmarg <- reformulate("offset(frailty.est2_)",ff[[2]]) # offset only
  
  randfit2<-coxme(ff.all,data=dat)
  
  frailty.est2<-(randfit2$frail$factor.frailty_) # aus coxme
  frailty.est2_<-apply(data.frame(frailty.est2),1,sum)[frailty]
  coxfit.fixed<-cph(ff, y=TRUE, surv=TRUE, 
                    method="breslow", type="kaplan-meier",
                    data=dat)
  coxfit.offs2<-cph(ff.all2, y=TRUE, surv=TRUE, 
                    method="breslow", type="kaplan-meier",
                    data=dat)
  Model_excl_frailty <- f.surev(coxfit.fixed) 
  Model_incl_frailty <- f.surev(coxfit.offs2) 

  results <- t(cbind(as.numeric(Model_excl_frailty[2:5]),
                     as.numeric(Model_incl_frailty[2:5])))
  colnames(results) <- c("D", "Dx", "V", "Vw")
  rownames(results) <- c("Fixed effects",
                         "Full model")
  # partial(!) EV for random effect as difference of the two output marginal EV's
  
  return(round(results,4))
}

#pev.frailty(formula(Surv(SURT,SURS)~LEV+INF),
#            frailty=osdata_ct$ctnum_, dat=osdata_ct)
#pev.frailty(formula(Surv(SURT,SURS)~LEV+INF),
#            frailty=osdata_ct$ctnum_, rterms="1+LEV", dat=osdata_ct)

#data(eortc) 
#pev.frailty(formula(Surv(y,uncens)~trt), frailty=eortc$center,
#            dat=eortc)
#pev.frailty(formula(Surv(y,uncens)~trt), frailty=eortc$center,
#            rterms="1+trt", dat=eortc)






