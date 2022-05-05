# Solve backward Fokker-Planck for May model with Babak's pseudocode Dec 28 2019 email
# Exit times and probability weighting from Appendix of Science manuscript
# see comparison of integration methods 2020-01-10; they get the same answer

rm(list = ls())
graphics.off()

library(bvpSolve)
library(cubature)
library(parallel)

options(mc.cores = parallel::detectCores())

# load information for boundary value problem; this version needs only avec, Drift.vec, TotSig.vec
# Save drift and diffusion curves
#save(Tstep,Xvar,avec,sigma,D1,D2,xeq,file=fname)
#
#
load('Mendota2021_DDbintau.Rdata')
Fname.out = c('Results_Step4ET_Mendota2021.Rdata')
piglab=c('Phycocyanin')

windows()
par(mfrow=c(1,1),mar=c(4,4.5,2,2)+0.1,cex.axis=1.6,cex.lab=1.6)
xrange=range(c(D1s$x,sigmas$x))
yrange = range(c(D1s$y,sigmas$y),na.rm=T)
plot(D1s$x,D1s$y,type='l',lwd=2,col='black',xlim=xrange,ylim=yrange,
     xlab=piglab,
     ylab='Drift and Diffusion')
points(sigmas$x,sigmas$y,type='l',lty=2,lwd=2,col='black')
grid()
abline(h=0,col='gray',lwd=2,lty=3)
legend('topleft',legend=c('Drift','Diffusion as s.d.'),lty=c(1,2),lwd=c(2,2),
       col=c('black','black'),bty='n',cex=1.6)

# Set up for bvp

# Interpolate with approxfun() for D1 and D2
# NOTE: This works if D1 and D2 do not have a sharp change in slope
D1fct = function(x) {
  yhat = approx(x=D1s$x, y=D1s$y,xout=x,method='linear',rule=2)$y
  return(yhat)
}

D2fct = function(x)  {
  yhat = approx(x=D2s$x, y=D2s$y,xout=x,method='linear',rule=2)$y
  return(yhat)
}

# Write D1 and D2 as functions; this is necessary if D1 or D2 have a sharp change in slope
# first trim D1 and D1 to eliminate missing D2
D1x = D1s$x
D1y = D1s$y
D2x = D2s$x
D2y = D2s$y
D1spline = smooth.spline(x=D1x,y=D1y)
D1fun = function(x) {
  #yhat = approx(x=avec, y=Drift.vec,xout=x,method='linear',rule=2)$y
  yhat=predict(D1spline,x)$y
  return(yhat)
}
D2spline = smooth.spline(x=D2x,y=D2y)
D2fun = function(x)  {
  #yhat = approx(x=avec, y=D2.from.Sig,xout=x,method='linear',rule=2)$y
  yhat=predict(D2spline,x)$y
  return(yhat)
}

# Find equilibria 
sdrift = sign(D1y)
dsdrift = c(0,-diff(sdrift))
xeq = D1x[which(!dsdrift == 0)]
ixeq = which(!dsdrift == 0)  # indices of the equilibria

print('equilibria & their indices',quote=F)
print(xeq)
print(ixeq)

Dratio = D1s$y/D2s$y
D1.over.D2 = function(x) {   # D2 = 0.5*sigma^2(x) for standard Langevin
  yhat = approx(x=D1s$x, y=Dratio, xout=x, method='linear',rule=2)$y
  return(yhat)
}

# boundaries for  attractors
#xlow = c(0.5,xeq[2])  # lower attractor
#xlow = c(xeq[2],xeq[3]+0.2)  # upper attractor

# plot drift and diffusion
#xvec = seq(avec[1]-0.1,xeq[3]+0.2,length.out=100) # original range
xvec = seq(-3.1,3.1,length.out=100)
drift = rep(0,100)
diff = rep(0,100)
for(i in 1:100) {
  drift[i] = D1fun(xvec[i])
  diff[i] = D2fun(xvec[i])
}

# Plot potential function
PF = rep(0,100) 
for(i in 2:100) {
  x0 = xvec[1]
  x1 = xvec[i]
  PF[i] = integrate(D1fun,lower=x0,upper=x1)$value
}
negPF = -1*PF

windows()
par(mfrow=c(3,1),mar=c(4,4.5,2,2)+0.1,cex.axis=2,cex.lab=2)
plot(xvec,drift,type='l',lwd=2,col='black',xlab='Chl Level',ylab='Drift')
abline(h=0,lty=2)
grid()
plot(xvec,negPF,type='l',lwd=2,col='black',xlab='Chl Level',ylab='-Potential')
grid()
plot(xvec,diff,type='l',lwd=2,col='black',xlab='Chl Level',ylab='Diffusion')
grid()

# Calculate and plot effective potential  ------------------------------------------
eprange = range(xvec)
xvec.ep = seq(eprange[1],eprange[2],length.out=100)

# Check D1/D2 
D1D2 = rep(0,100)
D1D2adj = rep(0,100)
for(i in 1:100) {
  D1D2[i] = D1.over.D2(xvec.ep[i])
  D1D2adj[i] = log(0.5*D2fun(xvec.ep[i])^2) - D1D2[i]
}

windows()
par(mfrow=c(1,1),mar=c(4,4,2,2)+0.1,cex.lab=1.5,cex.axis=1.5)
plot(xvec.ep,D1D2,type='l',lwd=2,col='darkred',xlab='x',ylab='D1/D2')
abline(h=0,lty=3,lwd=2)
grid()

windows()
par(mfrow=c(1,1),mar=c(4,4,2,2)+0.1,cex.lab=1.5,cex.axis=1.5)
plot(xvec.ep,D1D2adj,type='l',lwd=2,col='darkred',xlab='x',ylab='-(D1/D2)+log(D2)')
abline(h=0,lty=3,lwd=2)
grid()

# Plot effective potential function
EPF = rep(0,99)
for(i in 1:99) {
  x0 = xvec.ep[1]
  x1 = xvec.ep[i+1]
  xhalf = (x0 + x1)/2
  # These integrate commands return 'bad behavior of integrand'
  #integral = integrate(D1.over.D2,lower=x0,upper=x1)$value
  #integral = integrate(Vectorize(D1.over.D2),lower=x0,upper=x1)$value
  # This function from the cubature package seems to work
  integral = hcubature(f=D1.over.D2,lowerLimit=x0,upperLimit=x1)$integral
  logdiff = log(0.5*D2fun(xhalf)^2)  # D2 = 0.5*sigma^2 for standard Langevin
  EPF[i] = -1*integral + logdiff
}

windows(width=7,height=10)
par(mfrow=c(2,1),mar=c(4,4,2,2)+0.1,cex.lab=1.5,cex.axis=1.5)
plot(xvec,negPF,type='l',lwd=2,col='black',xlab='Chl Level',ylab='Potential')
plot(xvec.ep[2:100],EPF,type='l',lwd=2,col='black',xlab='Chl Level',ylab='Effective Potential')

windows(width=8,height=12)
par(mfrow=c(3,1),mar=c(4,4.5,2,2)+0.1,cex.axis=2,cex.lab=2)
plot(xvec,drift,type='l',lwd=2,col='black',xlab='Chl Level',ylab='Drift')
abline(h=0,lty=3)
abline(v=xeq,lty=2)
plot(xvec,diff,type='l',lwd=2,col='black',xlab='Chl Level',ylab='Diffusion')
abline(v=xeq,lty=2)
plot(xvec.ep[2:100],EPF,type='l',lwd=2,col='black',xlab='Chl Level',ylab='Effective Potential')
abline(v=xeq,lty=2)

# Calculate Exit Times =======================================================

# function for solving the boundary value problem as two differential equations 
#  for T (col 1) and dT/dx (col 2)
feval2 = function(x,y,plist) {
  out1 = y[2]
  out2 = -(D1fun(x)*y[2]+1)/D2fun(x)
  return( list(c(out1,out2)) )
}

feval2new = function(x,y,plist) { # D2() is sigma not D2 = 0.5*sigma^2
  out1 = y[2]
  out2 = -(D1fun(x)*y[2]+1)/(0.5*D2fun(x)^2)
  return( list(c(out1,out2)) )
}

# Left basin

# set up for bvpSolve
yini = c(NA,0)
yend = c(0,NA)

# solve the left basin from x = 0 (reflecting) to x=xeq[2] (absorbing)
#x = seq(xeq[1]-1,xeq[2],length.out=30)  # x vector, original
x = seq(-3.1,xeq[2],length.out=30)  # custom interval

# solve with bvpcol or bvptwp
trycol <- bvptwp(yini = yini, x = x, func = feval2, yend = yend, parm=plist)

windows()
par(mfrow=c(1,1),mar=c(3,4,3,2)+0.1,cex.lab=1.5,cex.axis=1.5)
plot(trycol[,1],trycol[,2],type='l',lty=1,lwd=2,col='blue',xlab='x',
     ylab='Exit Time',main='left basin')
abline(v=xeq[1],lty=1,col='magenta')
abline(v=xeq[2], lty=2, col='red')

# save solution for further analysis
ETL = trycol # left exit time

# end of solution for basin 1

# Right basin 

  # set up for bvpSolve
  yini = c(0, NA)
  yend = c(NA, 0)
  
  # right basin from x=xeq[2] (absorbing) to x > xeq[3] (reflecting)
  #x = seq(xeq[2],xeq[3]+1,length.out=30)  # x vector original
  x = seq(xeq[2],3.1,length.out=30)  # custom interval
  # solve with bvpcol or bvptwp
  trycol <- bvptwp(yini = yini, x = x, func = feval2, yend = yend, parm=plist)
  
  # save solution for further analysis
  ETR = trycol # right exit time
  
  windows()
  par(mfrow=c(1,1),mar=c(3,4,3,2)+0.1,cex.lab=1.5,cex.axis=1.5)
  plot(trycol[,1],trycol[,2],type='l',lty=1,lwd=2,col='blue',xlab='x',
       ylab='Exit Time',main = 'right basin')
  abline(v=xeq[2], lty=2, col='red')
  abline(v=xeq[3],lty=1,col='magenta')
  
# end of solution for basin 2
  
# CALCULATE WEIGHTS FOR AVERAGE EXIT TIME -----------------------------------------------------  
# Weight of x is p(x) from stationary distribution of Langevin equation  
# Based on appendix of Carpenter & Brock, Ecology Letters 2006 based on the book 
# 'Noise-Induced Transitions' by Horsthemke and Lefever 1984
  
# function for inner integral
finner = function (z) { 
  fx = D1fun(z)/(D2fun(z)) 
  return(fx)
}

# function for g^2 weights  
gg = function(z) {
  fx = 1/(D2fun(z))  
  return(fx)
}

# sigma for stationary distrib is not estimated so set it to something
sigma.sd = 1

# Calculate weighted average ET for both basins ===================================
# ETL[N,1] is the same as ETR[1,1]; they connect at xeq[2]
nL = length(ETL[,1])
nR = length(ETR[,1])
x = c(ETL[1,1]-0.01,ETL[2:nL,1],ETR[2:nR,1],ETR[nR,1]+0.01) # x has gaps that match ETL+ETR
dx = diff(x)
ETboth = c(ETL[1,2]-0.01,ETL[2:nL,2],ETR[2:nR,2],ETR[nR,2]+0.01)
nx = length(x)

# Check behavior of finner over range of x 
yfinner = finner(x)
windows()
plot(x,yfinner,type='l',lwd=2,col='darkred',xlab='x',ylab='finner')

# Weights by method in Science paper
# See 2020-08-20 version of this code for the H&L version
wtraw = rep(0,nx)
for(i in 2:(nx)) {
  intpart = hcubature(f=finner,lowerLimit=x[1],upperLimit=x[i])$integral
  epart = exp(intpart)
  wtraw[i] = epart*gg(x[i])
}

# normalize weights
wts = wtraw/sum(wtraw)

# artistic convenience on first weight
wts[1] = 0.98*wts[2]

windows()
par(mfrow=c(1,1),mar=c(4,4.5,2,2)+0.1,cex.axis=2,cex.lab=2)
plot(x,wts,type='l',lwd=2,col='black',xlab='Chl level',ylab='weight')
#     main='  left basin            right basin')
abline(v=xeq[1],lty=1,lwd=2,col='gray')
abline(v=xeq[2], lty=2,lwd=3,col='gray')
abline(v=xeq[3],lty=1,lwd=2,col='gray')
#text(1.5,0.04,'left basin',cex=1.5,font=2)
#text(2.2,0.04,'right basin',cex=1.5,font=2)

print('',quote=F)
print('weights',quote=F)
print(wts)

# Calculate mean exit time for left basin ---------------------------------------------

#
meanETl = sum(ETL[,2]*wts[1:nL])/sum(wts[1:nL])
print('',quote=F)
print('Mean ET for left basin, hours',quote=F)
print(meanETl/60)
print('-----------------------------------------------',quote=F)

# save axes
xL = x[1:nL]
wtsL = wts[1:nL]

# Calculate weighted average ET for right basin ===========================================================

#
meanETr = sum(ETR[,2]*wts[(nL+1):nx])/sum(wts[(nL+1):nx])
print('',quote=F)
print('Mean ET for right basin, hours',quote=F)
print(meanETr/60)
print('-----------------------------------------------',quote=F)

xR = x[(nL+1):nx]
wtsR = wts[(nL+1):nx]

# 2 panels, ET only
windows(width=10,height=6)
par(mfrow=c(1,2),mar=c(4,4.5,2,2)+0.1,cex.axis=2,cex.lab=2)
# top left = left basin
plot(ETL[,1],ETL[,2],type='l',lty=1,lwd=2,col='black',xlab='Chl Level',
     ylab='Exit Time, min',main='left basin')
abline(v=xeq[1],lty=1,lwd=2,col='gray')
abline(v=xeq[2], lty=2,lwd=2, col='gray')
# top right = right basin
plot(ETR[,1],ETR[,2],type='l',lty=1,lwd=2,col='black',xlab='Chl Level',
     ylab='Exit Time, min',main = 'right basin')
abline(v=xeq[2], lty=2,lwd=2,col='gray')
abline(v=xeq[3],lty=1,lwd=2,col='gray')

# 4 panel graph, including weights
windows()
par(mfrow=c(2,2),mar=c(4,4.5,2,2)+0.1,cex.axis=2,cex.lab=2)
# top left = left basin
plot(ETL[,1],ETL[,2],type='l',lty=1,lwd=2,col='black',xlab='Chl Level',
     ylab='Exit Time, min',main='left basin')
abline(v=xeq[1],lty=1,lwd=2,col='gray')
abline(v=xeq[2], lty=2,lwd=2, col='gray')
# top right = right basin
plot(ETR[,1],ETR[,2],type='l',lty=1,lwd=2,col='black',xlab='Chl Level',
     ylab='Exit Time, min',main = 'right basin')
abline(v=xeq[2], lty=2,lwd=2,col='gray')
abline(v=xeq[3],lty=1,lwd=2,col='gray')
# lower left = left weights
xrange=range(ETL[,1])
plot(xL,wtsL,type='l',lwd=2,xlim=xrange,xlab='Chl Level',ylab='weight')
abline(v=xeq[1],lty=1,lwd=2,col='gray')
abline(v=xeq[2], lty=2,lwd=2, col='gray')
# lower right = right weights
plot(xR,wtsR,type='l',lwd=2,xlab='Chl Level',ylab='weight')
abline(v=xeq[2], lty=2,lwd=2,col='gray')
abline(v=xeq[3],lty=1,lwd=2,col='gray')

# Save results for further plotting
# avec, drift.vec, TotSig.vec, xeq are results of DDJ 
# xvec, drift, diff, negPF are results of interpolation functions for next steps
# xvec.ep, EPF are effective potential
# x, wts are probability weights for the full range
# ETL and ETR are exit times for left and right basin
# xL, wtsL are weights for left basin
# xR, wtsR are weights for right basin
save(Tstep,Xvar,D1s,D2s,sigmas,xeq,xvec,drift,diff,negPF,
     xvec.ep,EPF,x,wts,ETL,ETR,xL,wtsL,xR,wtsR,
     file=Fname.out)
