# Test DDbintau
# SRC 2021-02-16

rm(list = ls())
graphics.off()

library(stats)
library(moments)

options(mc.cores = parallel::detectCores())

source('DDbintau+D4.r')

# save results
# Save post-processed DLM data
# X.raw is original series, X.z is z score, X.rawmean and X.rawsd are mean & sd for z score
# nl is DLM lags, delta is DLM delta, X.dlm is input to DLM, Tstep is time,
# local equilibria are X.eq and SD.eq, with ratio z.eq
# level (intercept) measures are level or stdlevel (for level/sd)
# Yyhat, B.ests, B.sd, errvar are DLM outputs
#save(nl,delta,X.raw,X.z,X.rawmean,X.rawsd,X.dlm,Tstep,X.eq,SD.eq,level,stdlevel,Z.eq,
#     Yyhat,B.ests,B.sd,errvar,file=Fname)
#
load(file='DLM_result_2021.Rdata')
fname=c('Mendota2021_DDbintau.Rdata')
plot.title=c('Mendota_2021')

Xvar = na.omit(stdlevel[10:length(stdlevel)])
Tstep = Tstep[1:length(Xvar)]  # if Xvar = DLM output then trim Tstep
SDXvar = sd(Xvar)
nxvar = length(Xvar)
dXvar = Xvar[2:nxvar] - Xvar[1:(nxvar-1)]
M2xvar = mean(dXvar^2)
M4xvar = mean(dXvar^4)
K.xvar = M4xvar/(M2xvar^2)
R.K.xvar = kurtosis(dXvar)

print('stats for Xvar',quote=F)
print(c('SD',SDXvar),quote=F)
print(c('M2',M2xvar),quote=F)
print(c('M4',M4xvar),quote=F)
print(c('K = M4/M2',K.xvar),quote=F)
print(c('R kurtosis formula ',R.K.xvar),quote=F)

windows(width=10,height=4)
par(mfrow=c(1,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5)
plot(Tstep,Xvar,type='l',lwd=2,xlab='day',
     ylab='std level')
grid()

windows(width=10,height=4)
par(mfrow=c(1,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5)
plot(Tstep[1:(nxvar-1)],diff(Xvar),type='l',lwd=2,xlab='day',
     ylab='diff(Xvar')
grid()

dxsort = sort(dXvar)
dxrank = c(1:length(dXvar))/nxvar
#windows(width=10,height=8)
#par(mfrow=c(1,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5)
#plot(dxrank,dxsort,type='l',lwd=2,xlab='rank',
#     ylab='sorted Xvar')
#grid()

# QQ plots
windows()
qqnorm(Xvar, pch = 19, cex=0.7, col='gray', main='NPP of X', frame = FALSE)
qqline(Xvar, col = "steelblue", lwd = 2)

# Set up for binning method
# DDbins = function(Xvar,bw,ntau,nbin)
title = c('Pigment variate')
ntau = 9
nbin = 100
bw <- 0.3*sd(Xvar)  # try between 0.1*sd and 0.5*sd 
# run function
DDout = DDbins(Xvar,bw,ntau,nbin)
# extract smoothed output
#outlist = list(D1s,D2s,sigmas,D4s,bin.mid,D1,D2,sigma,D4,SD) 
D1s = DDout[[1]]
D2s = DDout[[2]]
sigmas = DDout[[3]]
D4s = DDout[[4]]
bin.mid= DDout[[5]]
D1 = DDout[[6]]
D2 = DDout[[7]]
D4 = DDout[[9]]
SD = DDout[[10]]
SDall = DDout[[11]]

# Find equilibria
sdrift = sign(D1s$y)
dsdrift = c(0,-diff(sdrift))
xeq = D1s$x[which(!dsdrift == 0)]
ixeq = which(!dsdrift == 0)  # indices of the equilibria

print('equilibria and their indices',quote=F)
print(xeq,quote=F)
print(ixeq,quote=F)

# find crossings of xeq[2], the threshold
x.xeq2 = Xvar - xeq[2]
sx.xeq2 = sign(x.xeq2)
dsx.xeq2 = diff(sx.xeq2)
ixcross = which(!dsx.xeq2 == 0) # indices of crossings
dx.nocross = dXvar[-ixcross]

M2nocross = mean(dx.nocross^2)
M4nocross = mean(dx.nocross^4)
K.nocross = M4nocross/(M2nocross^2)

print('-------------------------------------------------------',quote=F)
print('stats for dXvar with no threshold crossings',quote=F)
print(c('M2nocross',M2nocross),quote=F)
print(c('M4nocross',M4nocross),quote=F)
print(c('K = M4/M2',K.nocross),quote=F)

# look at results of DDbintau
windows()
par(mfrow=c(3,1),mar=c(4,4.5,2,2)+0.1,cex.axis=2,cex.lab=2)
plot(D1s$x,D1s$y,type='l',lwd=2,col='blue',xlab='state',ylab='D1')
abline(h=0,lty=2)
plot(sigmas$x,sigmas$y,type='l',lwd=2,col='red',xlab='state',ylab='sigma')
plot(D2s$x,D2s$y,type='l',lwd=2,col='red',xlab='state',ylab='D2')

windows()
par(mfrow=c(1,1),mar=c(4,4.5,2,2)+0.1,cex.axis=2,cex.lab=2)
yrange = range(c(D2s$y^2,D4s$y),na.rm=T)
plot(D2s$x,D2s$y^2,type='l',lwd=2,col='red',ylim=yrange,log='y',
     xlab='state',ylab='D2^2 & D4',main='smoothed')
points(D4s$x,D4s$y,type='l',lwd=2,col='blue',xlab='state',ylab='D2')
legend('topright',legend=c('D2^2','D4'),lwd=c(2,2),col=c('red','blue'))
grid()

windows()
par(mfrow=c(1,1),mar=c(4,4.5,2,2)+0.1,cex.axis=2,cex.lab=2)
yrange = range(c(D2^2,D4),na.rm=T)
plot(bin.mid,D2^2,type='b',lwd=1,pch=19,cex=0.7,col='red',ylim=yrange,log='y',
     xlab='state',ylab='D2^2 & D4 by Bin for dt=1',main=plot.title)
points(bin.mid,D4,type='b',lwd=1,pch=19,cex=0.7,col='blue',xlab='state',ylab='D2')
legend('bottomleft',legend=c('D2^2','D4'),lwd=c(2,2),col=c('red','blue'))
grid()

windows()
par(mfrow=c(1,1),mar=c(4,4.5,2,2)+0.1,cex.axis=2,cex.lab=2)
yrange = range(c(D2,D4),na.rm=T)
plot(bin.mid,D2,type='b',lwd=1,pch=19,cex=0.7,col='red',ylim=yrange,log='y',
     xlab='state',ylab='D2 & D4 by Bin for dt=1',main=plot.title)
points(bin.mid,D4,type='b',lwd=1,pch=19,cex=0.7,col='blue',xlab='state',ylab='D2')
legend('bottomleft',legend=c('D2','D4'),lwd=c(2,2),col=c('red','blue'))
grid()


# scale D2 and D4 by SD
D2scale = D2/(SD^2)
D4scale = D4/(SD^4)

windows()
par(mfrow=c(1,1),mar=c(4,4.5,2,2)+0.1,cex.axis=2,cex.lab=2)
yrange = range(c(D2scale,D4scale),na.rm=T)
plot(bin.mid,D2scale,type='p',pch=19,cex=0.7,log='y',ylim=yrange,col='red',
     xlab='state',ylab='D2scale red, D4scale blue', main='scale by SD')
points(bin.mid,D4scale,type='p',pch=19,cex=0.7,col='blue')

save(Tstep,Xvar,bin.mid,D1s,D2s,sigmas,xeq,file=fname)
