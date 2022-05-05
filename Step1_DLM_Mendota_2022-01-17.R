# Trial of general ET programs - step 1 is DLM
# SRC 2020-01-27

rm(list = ls())
graphics.off()

#library(Langevin)
#library(plotrix)
library(parallel)
#library(svMisc)  # used only if eigenvalues are bootstrapped

#source('DriftDiffJumpFunction.r')
source('ODLMAR_NoBoot_2018-10-20.R')

options(mc.cores = parallel::detectCores())

# Load data
load(file='Me_BGA+cov_1min_2021.Rdata')
title = c('Mendota 2021 lBGA')
Fname = c('DLM_result_2021.Rdata')

print(MeBGA.nona[1,])

# Subset to days of summer stratification, 15 May (135) or 1 June (152) to 15 Sept (258)
MeBGAs = subset(MeBGA.nona,subset=(DOY >= 152 & DOY <= 258))
# Trim sensor interruptions: lBGA < 2
#MeBGAs$lBGA = ifelse(MeBGAs$lBGA < 2,2,MeBGAs$lBGA)
# Save trimmed data to .csv
#write.csv(MeBGAs,file='Me_1min_2011trim.csv')

windows(width=12,height=6)
par(mfrow=c(1,1),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(MeBGAs$DOY,MeBGAs$lBGA,type='l',
     xlab='Day of Year',ylab='log10 BGA',ylim=c(1.5,4.5),main=title)
grid()

# Select variates
Tscore = MeBGAs$DOY
pig = MeBGAs$lBGA

# Start DLM
# z score function
zscore=function(x) {
  xbar = mean(x,na.rm=T)
  xsd = sd(x,na.rm=T)
  z = (x - xbar)/xsd
  return(z)
}

# optionally convert to logs and take zscore
Tstep = Tscore
X.raw = pig
X.z = zscore(X.raw)
X.rawmean = mean(X.raw,na.rm=T)
X.rawsd = sd(X.raw,na.rm=T)
print(c('mean & s.d. of X.raw ',X.rawmean,X.rawsd),quote=F)

# Select input to dlm
X.dlm = X.z

windows(width=12,height=6)
plot(Tstep,X.dlm,type='l',col='forestgreen',xlab='DoY index',ylab='x.dlm',
     main='time series for DLM')
grid()

# Set up DLM
nobs = length(X.dlm)
nl = 1 # number of lags
delta = 0.95 # 0<delta<1; see advice in functions

# Run DLM
ODL.out = ODLMAR(nl,delta,X.dlm,Tstep,title)

# Output matrices are stored sideways, like MARSS
Yyhat = ODL.out[[1]]
EigenVals = ODL.out[[2]]
B.ests = ODL.out[[3]]  
B.sd = ODL.out[[4]]
errvar = ODL.out[[5]] # updated error variance

# Post process DLM -----------------------------------------------

# Calculate moving equilibrium
X.eq = B.ests[1,]/(1 - B.ests[2,])
# Calculate its variance
SDterm1 = X.eq*X.eq
SDterm2 = (B.sd[1,]*B.sd[1,] + errvar)/(B.ests[1,]*B.ests[1,])
SDterm3 = (B.sd[2,]*B.sd[2,])/((1 - B.ests[2,])*(1 - B.ests[2,]))
SD.eq = sqrt(SDterm1*(SDterm2 + SDterm3))
# Z score
Z.eq = X.eq/SD.eq

# Time steps start at 2
Nstep = length(Tstep)

# Plot components of steady-state estimate
windows(width=6,height=12)
par(mfrow=c(3,1),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(Tstep[2:Nstep],X.eq,type='l',col='blue',ylim=c(-10,10),
     ylab='Steady State',xlab='DoY',
     main='Local Steady-State estimate, sd, and ratio')
grid()
plot(Tstep[2:Nstep],SD.eq,type='l',col='red',ylim=c(0,10),
     ylab='S.D.',xlab='DoY')
grid()
plot(Tstep[2:Nstep],Z.eq,type='l',col='purple',
     ylab='Z score',xlab='DoY')
grid()

# Calculate level estimates
level = B.ests[1,]
stdlevel = B.ests[1,]/B.sd[1,]

# Plot components of level estimate
windows(width=12,height=9)
par(mfrow=c(2,1),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(Tstep[2:Nstep],level,type='l',col='blue',#ylim=c(-10,10),
     ylab='level',xlab='DoY',
     main='Level and Std level estimate')
grid()
plot(Tstep[2:Nstep],stdlevel,type='l',col='red',#ylim=c(0,10),
     ylab='Std Level',xlab='DoY')
grid()

# analyze subset
xrange = c(170,180)
windows(width=12,height=12)
par(mfrow=c(5,1),mar=c(3, 4.2, 2, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(Tstep[2:Nstep],stdlevel,xlim=xrange,type='l',col='red',#ylim=c(0,10),
     ylab='Std Level',xlab='DoY')
grid()
plot(Tstep,X.z,xlim=xrange,type='l',col='ForestGreen',#ylim=c(0,10),
     ylab='X.z',xlab='DoY')
grid()
plot(Tstep,MeBGAs$v.wind,xlim=xrange,type='l',col='blue',#ylim=c(0,10),
     ylab='V.wind',xlab='DoY')
grid()
plot(Tstep,MeBGAs$d.wind,xlim=xrange,type='l',col='blue',#ylim=c(0,10),
     ylab='d.wind',xlab='DoY')
grid()
plot(Tstep,MeBGAs$wtemp,xlim=xrange,type='l',col='magenta',#ylim=c(0,10),
     ylab='wtemp',xlab='DoY')
grid()

# Save post-processed DLM data
# X.raw is original series, X.z is z score, X.rawmean and X.rawsd are mean & sd for z score
# nl is DLM lags, delta is DLM delta, X.dlm is input to DLM, Tstep is time,
# local equilibria are X.eq and SD.eq, with ratio z.eq
# level (intercept) measures are level or stdlevel (for level/sd)
# Yyhat, B.ests, B.sd, errvar are DLM outputs
save(MeBGAs,nl,delta,X.raw,X.z,X.rawmean,X.rawsd,X.dlm,Tstep,X.eq,SD.eq,
     level,stdlevel,Z.eq,Yyhat,B.ests,B.sd,errvar,file=Fname)
print(Fname,quote=F)

