# Test DDbintau
# SRC 2021-02-16

rm(list = ls())
graphics.off()

library(stats)
library(parallel)

source('DDbintau+D4.r')

options(mc.cores = parallel::detectCores())

# Results of DLM
#save(Nboot,Bootlevel,Bootsdlevel,file='DLM_boot_output.Rdata')
# first column of Bootsdlevel is doy, next 100 columns are bootstrap DLM series

load(file='DLM_boot_2021.Rdata')
fname=c('DDbintau_boot_2021.Rdata')

#load(file='DLM_boot_Tuesday1315.Rdata')
#fname=c('DDbintau_boot_Tuesday1315.Rdata')

# set constants for the run
title = c('StdLevel_LogPhycocyanin')
ntau = 9
nbin = 100

# define output arrays
D1ar = array(0,dim=c(Nboot,nbin,2)) # Nboot*nbin*2 matrices
D2ar = D1ar
sigar = D1ar
xeqmat = matrix(0,nr=Nboot,3)

tstart = Sys.time()
ngoodboot=0
for(iboot in 1:Nboot) {
  Xvar0 = Bootsdlevel[ ,(iboot+1)]
  #Xvar = na.omit(Xvar0) # usual filter- removes NA
  Xvar1 = Xvar0[is.finite(Xvar0)] # removes infinite values and NA
  Xvar = subset(Xvar1,subset=(abs(Xvar1) < 20))
  #bw <- 0.3*sd(Xvar)  # mesh method: try between 0.1*sd and 0.5*sd 
  bw = 0.1*diff(range(Xvar)) # tie bw to range of time series for bins
  # run function
  DDout = DDbins(Xvar,bw,ntau,nbin)
  # extract smoothed output
  D1s = DDout[[1]]
  D2s = DDout[[2]]
  sigmas = DDout[[3]]
  # screen for quality - reject runs with missing bins
  nnona = length(na.omit(D1s$y))
  if(nnona < nbin) { next }
  nnona = length(na.omit(D2s$y))
  if(nnona < nbin) { next }
  nnona = length(na.omit(sigmas$y))
  if(nnona < nbin) { next }
  # check for 3 roots
  sdrift = sign(D1s$y)
  dsdrift = c(0,-diff(sdrift))
  xeq = D1s$x[which(!dsdrift == 0)]
  if(length(xeq) != 3) { next }
  xeqmat[iboot,] = xeq
  # fill output arrays
  D1ar[iboot, ,1] = D1s$x
  D1ar[iboot, ,2] = D1s$y
  D2ar[iboot, ,1] = D2s$x
  D2ar[iboot, ,2] = D2s$y
  sigar[iboot, ,1] = sigmas$x
  sigar[iboot, ,2] = sigmas$y
  ngoodboot = ngoodboot+1
}

# Find xeq rows with less than 3 equilibria
test1 = rep(Nboot,1)
for(ib in 1:Nboot) {
  test1[ib] = sum(xeqmat[ib,])
}
OKboot = which(test1 != 0)
NOK = length(OKboot)
print(c('Number of bootstraps with 3 roots',NOK),quote=F)

# extract xeq, D1ar, and D2ar arrays with only 3-root cases
xeqmat.OK = xeqmat[OKboot, ]
D1ar.OK = D1ar[OKboot, , ]
D2ar.OK = D2ar[OKboot, , ]
sigar.OK = sigar[OKboot, , ]

print('confirm dimensions of arrays with only 3-root cases',quote=F)
print(c(dim(xeqmat.OK),'xeqmat.OK'),quote=F)
print(c(dim(D1ar.OK),'D1ar.OK'),quote=F)
print(c(dim(D2ar.OK),'D2ar.OK'),quote=F)
print(c(dim(sigar.OK),'sigar.OK'),quote=F)

windows()
plot(D1ar.OK[1, ,1],D1ar.OK[1, ,2],type='l',lwd=1,col='blue',
     xlab='pigment',ylab='D1')
for(i in 2:NOK) {
  points(D1ar.OK[i, ,1],D1ar.OK[i, ,2],type='l',lwd=1,col='blue')
}
abline(h=0,lty=2,lwd=2,col='darkred')
abline(v=-3,lty=3,lwd=2)
abline(v=3,lty=3,lwd=2)
grid()

windows()
plot(D2ar.OK[1, ,1],D2ar.OK[1, ,2],type='l',lwd=1,col='blue',
     xlab='pigment',ylab='D2')
for(i in 2:NOK) {
  points(D2ar.OK[i, ,1],D2ar.OK[i, ,2],type='l',lwd=1,col='blue')
}

windows()
plot(sigar.OK[1, ,1],sigar.OK[1, ,2],type='l',lwd=1,col='blue',
     xlab='pigment',ylab='sigma')
for(i in 2:NOK) {
  points(sigar.OK[i, ,1],sigar.OK[i, ,2],type='l',lwd=1,col='blue')
}

tstop = Sys.time()
print('',quote=F)
print('----------------------------------------------------------',quote=F)
runtime = difftime(tstop,tstart,units='mins')
print(c('runtime, minutes ',runtime),quote=F)
print(c('number of bootstraps that were OK ',ngoodboot),quote=F)
save(Nboot,nbin,D1ar,D2ar,sigar,xeqmat,OKboot,NOK,xeqmat.OK,
     D1ar.OK,D2ar.OK,sigar.OK,file=fname)

