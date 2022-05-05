# Trial of general ET programs - step 1 is DLM
# SRC 2020-01-27

rm(list = ls())
graphics.off()

library('moments')
library('parallel')

source('ODLMAR_for_Bootstrap_2020-11-23.R')

options(mc.cores = parallel::detectCores())

# Input result of nominal DLM and set up bootstrap

# save results
# Save post-processed DLM data
# X.raw is original series, X.z is z score, X.rawmean and X.rawsd are mean & sd for z score
# nl is DLM lags, delta is DLM delta, X.dlm is input to DLM, Tstep is time,
# local equilibria are X.eq and SD.eq, with ratio z.eq
# level (intercept) measures are level or stdlevel (for level/sd)
# Yyhat, B.ests, B.sd, errvar are DLM outputs
#save(nl,delta,X.raw,X.z,X.rawmean,X.rawsd,X.dlm,Tstep,X.eq,SD.eq,level,stdlevel,Z.eq,
#     Yyhat,B.ests,B.sd,errvar,file=Fname)

load(file='DLM_Result_2021.Rdata')

epsilon = Yyhat[,3] - Yyhat[,2]
Yoriginal = Yyhat[,2] # save the original Y with another name
Yhat.nominal = Yyhat[,3] # save Yhat from nominal model

windows()
plot(Yyhat[,1],epsilon,type='l',lwd=1,col='blue',xlab='Time step',ylab='epsilon')

print('Descriptive stats of errors',quote=F)
print(summary(epsilon))
print('N, S.D., skewness, kurtosis of epsilon',quote=F)
print(c(length(epsilon),sd(epsilon),skewness(epsilon),kurtosis(epsilon)),quote=F)

windows(width=8,height=5)
par(mfrow=c(1,2),mar=c(4, 4.3, 4, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
eps.acf=acf(epsilon,lag.max=10)
eps.pacf=pacf(epsilon,lag.max=10)

print('',quote=F)
print(eps.acf)

Nboot = 100  # number of bootstrap cycles
Tstep = Yyhat[,1]
Ndlm = length(Yyhat[,1]) - 1
Tstep = Yyhat[1:Ndlm,1]
Bootlevel = matrix(0,nr=Ndlm,nc=(1+Nboot))
Bootlevel[,1] = Tstep
Bootsdlevel = Bootlevel

# DLM run details
nobs = length(X.dlm)
title = c('Bootstrap')

# BOOTSTRAP
tstart = Sys.time()

for(i in 1:Nboot) {
  eps.rand = sample(epsilon,replace=T) # randomize eps
  Ypsuedo = Yhat.nominal + eps.rand
  print(c('boot cycle ',i),quote=F)
  ODL.out = ODLMAR(nl,delta,Ypsuedo,Tstep,title)
  # Output matrices are stored sideways, like MARSS
  #Yyhat = ODL.out[[1]]
  #EigenVals = ODL.out[[2]]  # only if eigenvalues computed
  B.ests = ODL.out[[2]]    # 3 if there are eigenvalues
  B.sd = ODL.out[[3]]      # 4 if eigenvalues
  #errvar = ODL.out[[5]] # updated error variance; 5 if eigenvalues
  level = B.ests[1,]
  stdlevel = B.ests[1,]/B.sd[1,]
  Bootlevel[,(i+1)] = level
  Bootsdlevel[,(i+1)] = stdlevel
}

tstop = Sys.time()
runtime = tstop-tstart
print(c('Bootstrap run time = ',runtime),quote=F)

# Plot results
windows(width=12,height=6)
yrange=range(Bootlevel[,2:(Nboot+1)],na.rm=T)
plot(Bootlevel[,1],Bootlevel[,2],type='l',lwd=1,col='blue',ylim=yrange,
     xlab='timestep',ylab='level')
points(Bootlevel[,1],Bootlevel[,3],type='l',lwd=1,col='blue')
for(iplot in 4:(Nboot+1)) {
  points(Bootlevel[,1],Bootlevel[,iplot],type='l',lwd=1,col='blue')
}

windows(width=12,height=6)
#xx1=Bootsdlevel[is.finite(rowSums(Bootsdlevel)),]
#xx2=subset(xx1,subset=(rowSums(xx1)<=1.e4))
yrange=range(Bootsdlevel[,2:(Nboot+1)],na.rm=T)
plot(Bootsdlevel[,1],Bootsdlevel[,2],type='l',lwd=1,col='blue',ylim=yrange,
     xlab='timestep',ylab='std level')
points(Bootsdlevel[,1],Bootsdlevel[,3],type='l',lwd=1,col='blue')
for(iplot in 4:(Nboot+1)) {
  points(Bootsdlevel[,1],Bootsdlevel[,iplot],type='l',lwd=1,col='blue')
}


# Level statistics vs pigment
levelmean = apply(Bootlevel[,2:(Nboot+1)],1,mean,na.rm=T)
levelsd = apply(Bootlevel[,2:(Nboot+1)],1,sd,na.rm=T)

windows()
plot(Bootlevel[,1],levelmean,type='l',lwd=1,col='blue',xlab='timestep',ylab='mean level')

windows()
plot(Bootlevel[,1],levelsd,type='l',lwd=1,col='red',xlab='timestep',ylab='sd of level')

levelratio = levelmean/levelsd

windows()
plot(Bootlevel[,1],levelratio,type='l',lwd=1,col='purple',xlab='timestep',ylab='mean/sd for level')

save(Nboot,Bootlevel,Bootsdlevel,levelsd,levelmean,file='DLM_boot_2021.Rdata')

