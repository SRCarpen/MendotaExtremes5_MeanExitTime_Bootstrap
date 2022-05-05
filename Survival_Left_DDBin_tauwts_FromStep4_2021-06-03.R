# Attempt Survival Function with DeSolve
# SRC using Babak notes 23 Sept 2020

rm(list = ls())
graphics.off()

library('deSolve')
library('grDevices')
library('numDeriv')

options(mc.cores = parallel::detectCores())

# Load results of Step 4

# Save results for further plotting
# avec, drift.vec, TotSig.vec, xeq are results of DDJ 
# xvec, drift, diff, negPF are results of interpolation functions for next steps
# xvec.ep, EPF are effective potential
# x, wts are probability weights for the full range
# ETL and ETR are exit times for left and right basin
# xL, wtsL are weights for left basin
# xR, wtsR are weights for right basin
#save(Tstep,Xvar,D1s,D2s,sigmas,xeq,xvec,drift,diff,negPF,
#     xvec.ep,EPF,x,wts,ETL,ETR,xL,wtsL,xR,wtsR,
#     file=Fname.out)

load(file='Results_Step4ET_Mendota2021.Rdata')
fname = c('Survival_Left_Mendota2021.Rdata')

windows()
par(mfrow=c(1,1),mar=c(4,4.5,2,2)+0.1,cex.axis=1.6,cex.lab=1.6)
xrange=range(c(D1s$x,sigmas$x),na.rm=T)
yrange = range(c(D1s$y,sigmas$y),na.rm=T)
plot(D1s$x,D1s$y,type='l',lwd=2,col='black',xlim=xrange,ylim=yrange,
     xlab='Chlorophyll Level',
     ylab='Drift and Diffusion')
points(sigmas$x,sigmas$y,type='l',lty=2,lwd=2,col='black')
abline(h=0,col='gray',lwd=2,lty=3)
legend('topright',legend=c('Drift','sigma'),lty=c(1,2),lwd=c(2,2),
       col=c('black','black'),bty='n',cex=1.6)

# calculate mean exit times
nL = length(ETL[,1])
nR = length(ETR[,1])
nx = length(x)
#
meanETl = sum(ETL[,2]*wts[1:nL])/sum(wts[1:nL])
print('',quote=F)
print('Mean ET for left basin, hours',quote=F)
print(meanETl/60)
print('-----------------------------------------------',quote=F)

#
meanETr = sum(ETR[,2]*wts[(nL+1):nx])/sum(wts[(nL+1):nx])
print('',quote=F)
print('Mean ET for right basin, hours',quote=F)
print(meanETr/60)
print('-----------------------------------------------',quote=F)


# Write D1 and D2 as functions 
D1dat = na.omit( as.data.frame(cbind(D1s$x,D1s$y)) )
D1spline = smooth.spline(x=D1dat[,1],y=D1dat[,2])
D1fct = function(x) {
  #yhat = approx(x=avec, y=Drift.vec,xout=x,method='linear',rule=2)$y
  yhat=predict(D1spline,x)$y
  return(yhat)
}
D2dat = na.omit( as.data.frame(cbind(D2s$x,D2s$y)) )
D2spline = smooth.spline(x=D2dat[,1],y=D2dat[,2])
D2fct = function(x)  {
  #yhat = approx(x=avec, y=D2.from.Sig,xout=x,method='linear',rule=2)$y
  yhat=predict(D2spline,x)$y
  return(yhat)
}

# Survival curve for left basin

# Make grid
N = 1000
# derived from equilibria
#xs = seq(from=xeq[1]-0.5,to=(xeq[2]-0.01),length.out=N)
#dx = diff(range(xs))/N
# derived from probability weights
xrange = range(ETL[,1],na.rm=T)
xs = seq(from=xrange[1],to=xrange[2],length.out=N)
dx = diff(range(xs))/N

# Langevin components on mesh
drif = D1fct(xs)
dif = D2fct(xs)

# Derivatives on mesh
#dD1 = grad(func=D1,x=xs,method='Richardson') # first derivative D1 on mesh
#dD2 = grad(func=D2,x=xs,method='Richardson') # first derivative D2 on mesh
#d2D2 = rep(0,N)  # vector for second derivative D2 on mesh
#for(i in 1:N) {
#  d2D2[i] = hessian(func=D2,x=xs[i],method='Richardson')  
#}

# test plots
windows() #windows(height=12,width=6)
par(mfrow=c(2,1),mar=c(4,4.5,2,2)+0.1,cex.axis=1.6,cex.lab=1.6)
plot(xs,drif,type='l')
plot(xs,dif,type='l')
#plot(xs,dD1,type='l')
#plot(xs,dD2,type='l')
#plot(xs,d2D2,type='l')

Pars = c(drif=drif,
          dif=dif,
          #dD1=dD1,
          #dD2=dD2,
          #d2D2=d2D2,
          dx=dx)

Bonehead = function(time,state,parms,N) {
  with ( as.list(parms), {
    S=state
    #S1 = ifelse(state>1,1,state) # prevent S>1
    #S = ifelse(S1<0,1.e-5,S1) # prevent S<0
    S[N] = 0 # impose right boundary condition; length remains N
    dS = c(0,diff(S)) # include left boundary condition; length remains N
    dSdx = dS/dx
    #d2S = c(dS[2]-dS[1],diff(dS)) # duplicate left-side value, length remains N
    d2S = c(diff(dS),(dS[N]-dS[N-1])) # duplicate right-side value, length remains N
    d2Sdx2 = d2S/(dx*dx)
    # calculate rate using the survival formula directly
    dS = drif*dSdx + dif*d2Sdx2  # drif is D1, dif is D2
    dSout = list(as.vector(dS))
    return(dSout)
  } )
}

# set the duration and intervals
TMAX = 1000
times = seq(from=0,to=TMAX,by=0.1)

# set the initial conditions
Sini = rep(1,N)  # intial S=1
#Sini = runif(N,min=0.99,max=1) # random initial near 1 

# Run the model
tstart = Sys.time()
out = ode.1D(y=Sini,times=times,func=Bonehead,parms=Pars,nspec=1, N=N)
tstop = Sys.time()
runtime = tstop-tstart
print(c('Run time = ',runtime),quote=F)

# set up for plots
# rows are times, columns are pigment values
#Smat = out[,2:(N+1)]  # full Smat
# subset for large number of time steps
tsub = times #seq(1,length(times),1000) # time step values for subset
tsub = tsub/60 # convert time axis to hours
Smat = out[,2:(N+1)] # out[times,2:(N+1)]
print(c('range of S = ',round(range(Smat),2)),quote=F)

#windows()
#filled.contour(x=times,y=xs,z=Smat,color=gray.colors,xlab='Time',ylab='Pigment')

# Choose pigment values to plot
ivalues = c(2,round(N/4),round(N/2),round(0.8*N),(N-1))
windows(height=12,width=5)
par(mfrow=c(5,1),mar=c(4,4.5,2,2)+0.1,cex.axis=1.6,cex.lab=1.6)
for(i in 1:5) {
  inow = ivalues[i]
  pigment = xs[inow]
  plot(tsub,Smat[,inow],type='l',xlab='time',ylab='S',log='y',#ylim=c(-0.1,1),
                  main=paste("pigment= ",round(pigment,2)) )
  abline(h=0.5,lty=2,lwd=2,col='red')
}

# surface plot functions (slow) ======================================================================================================
leftsurf = function(tsub,xs,Smat) {
  windows()
  par(mfrow=c(1,1),mar=c(4,4.5,2,2)+0.1,cex.axis=1.6,cex.lab=1.6)
  zvals = c(0.1,0.3,0.5,0.7,0.9)
  opt = options(max.contour.segments=200000)
  # Left basin
  contour(x=tsub,y=xs,z=Smat,xlab='Time, hours',
          # Right basin
          #zrevcol = apply(out[,(2:length(times)+1)],2,rev)  # flip z matrix on y axis for right basin only
          #contour(x=out[,1],y=rev(xgrid$x.mid),z=zrevcol,xlab='Time, minutes',
          ylab='Phycocyanin',main='Survival Function',
          levels=zvals,labcex=1,opt)
}

rightsurf = function(tsub,xs,Smat) {
  windows()
  par(mfrow=c(1,1),mar=c(4,4.5,2,2)+0.1,cex.axis=1.6,cex.lab=1.6)
  zvals = c(0.1,0.3,0.5,0.7,0.9)
  opt = options(max.contour.segments=200000)
  # Left basin
  #contour(x=tsub,y=xs,z=Smat,xlab='Time, minutes',
          #Right basin
          zrevcol = apply(Smat,2,rev)  # flip z matrix on y axis for right basin only
          contour(x=out[,1],y=rev(xgrid$x.mid),z=zrevcol,xlab='Time, minutes',
          ylab='Phycocyanin',main='Survival Function',
          levels=zvals,labcex=1,opt)
}
# End surface functions ================================================================================================================

# Select surface plot if needed (slow)
leftsurf(tsub,xs,Smat)
#rightsurf(tsub,xs,Smat)

# Half-life
# Left basin
options("max.contour.segments"=200000)
zhalf0 = contourLines(x=tsub,y=xs,z=Smat,nlevels=1,levels=0.5)
# Right basin
#zrevcol = apply(Smat,2,rev)  # flip z matrix on y axis for right basin only
#zhalf0 = contourLines(x=out[,1],y=rev(xgrid$x.mid),z=zrevcol,nlevels=1,levels=zvals)
zhalf1 = as.data.frame(zhalf0)

# extract a subset for plotting
nhalf = length(zhalf1$x)
halfkeep = seq(1,nhalf,length.out=100)
zhalf = zhalf1[halfkeep,]

windows()
par(mfrow=c(1,1),mar=c(4,4.5,2,2)+0.1,cex.axis=1.6,cex.lab=1.6)
plot(zhalf$y,zhalf$x,type='l',lwd=2,xlab='phycocyanin',ylab='S(x,t) = 0.5 (hours)')
abline(v=xeq,lty=c(1,2,1),col=c('blue','red','blue'))

# Two panels - contour() cannot share a window
windows()
par(mfrow=c(1,1),mar=c(4,4.5,2,2)+0.1,cex.axis=1.6,cex.lab=1.6)
zvals = c(0.1,0.3,0.5,0.7,0.9)
opt = options(max.contour.segments=200000)
contour(x=xs,y=tsub,z=t(Smat),xlab='Phycocyanin',col='blue',lwd=2,
        # Right basin
        #zrevcol = apply(out[,(2:length(times)+1)],2,rev)  # flip z matrix on y axis for right basin only
        #contour(x=out[,1],y=rev(xgrid$x.mid),z=zrevcol,xlab='Time, minutes',
        ylab='Time (hours)',main='Clear Water Survival Probability',
        levels=zvals,labcex=1,opt)
abline(v=xeq,lty=c(2,2,2),lwd=c(2,2,2),col=c('blue','red','blue'))

windows()
par(mfrow=c(1,1),mar=c(4,4.5,2,2)+0.1,cex.axis=1.6,cex.lab=1.6)
plot(zhalf$y,zhalf$x,type='l',lwd=2,col='blue',xlab='phycocyanin',ylab='Half-Life, hours',
     main='Clear-Water Half Life')
abline(v=xeq,lty=c(2,2,2),lwd=c(2,2,2),col=c('blue','red','blue'))

# compute weighted mean survival time
# Write stationary density (weights) as a function 
Wspline = smooth.spline(x=x,y=wts)
Wfct = function(x) {
  #yhat = approx(x=avec, y=Drift.vec,xout=x,method='linear',rule=2)$y
  yhat=predict(Wspline,x)$y
  return(yhat)
}

# weighted mean; remember that x and y are reversed for zhalf
NS = length(zhalf$y)
SW = rep(0,NS)
SxW = rep(0,NS)
for(i in 1:NS) {
  SW[i] = Wfct(zhalf$y[i])  # weight
  SxW[i] = zhalf$x[i]*SW[i] # piece of weighted avg
}

Shalf.wts = sum(SxW)/sum(SW)
print(c('Weighted Half-Life ',Shalf.wts),quote=F)

save(tsub,xs,Smat,xeq,zhalf,D1s,D2s,sigmas,x,wts,Shalf.wts,
     file=fname)