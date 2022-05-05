# Solve BVP for ET and integrate the density for each bootstrapped Langevin eq
# SRC 2021-02-16

rm(list = ls())
graphics.off()

library('stats')
library('parallel')
library('bvpSolve')
library('cubature')

library('deSolve')
library('grDevices')
library('numDeriv')

options(mc.cores = parallel::detectCores())

source('DDbintau+D4.r')

# Functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

Bonehead.left = function(time,state,parms,N) {
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
} # end Bonehead (function for ode1D)

Bonehead.right = function(time,state,parms,N) {
  with ( as.list(parms), {
    S=state
    #S1 = ifelse(state>1,1,state) # prevent S>1
    #S = ifelse(S1<0,1.e-5,S1) # prevent S<0
    S[1] = 0 # impose right boundary condition; length remains N
    dS = c(diff(S),0 ) # include left boundary condition; length remains N
    dSdx = dS/dx
    d2S = c(dS[2]-dS[1],diff(dS)) # duplicate left-side value, length remains N
    #d2S = c(diff(dS),(dS[N]-dS[N-1])) # duplicate right-side value, length remains N
    d2Sdx2 = d2S/(dx*dx)
    # calculate rate using the survival formula directly
    dS = drif*dSdx + dif*d2Sdx2  # drif is D1, dif is D2
    dSout = list(as.vector(dS))
    return(dSout)
  } )
} # end Bonehead (function for ode1D)

# END FUNCTIONS  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

# Enter Data - Results of bootstrap DDbintau ----------------------------------------------------
#save(Nboot,nbin,D1ar,D2ar,sigar,xeqmat,OKboot,NOK,xeqmat.OK,
#    D1ar.OK,D2ar.OK,sigar.OK,file=fname)
load(file='DDbintau_boot_2021.Rdata')
fname=c('ET+Sboot_2021.Rdata')
# or
#load('DDbintau_Tuesday1315.Rdata')
#fname=c('ET+Sboot_Tuesday1315.Rdata')

# Simulation control parameters
tstart = Sys.time() # start clock

ETLR = matrix(0,nr=NOK,nc=2)
ShalfLR = matrix(0,nr=NOK,nc=2)

ETfun = for(iOK in 1:NOK)  {  # usually 1:NOK
  ib = OKboot[iOK]
  xeq = xeqmat.OK[iOK,]
  D1 = D1ar.OK[iOK, , ]
  D2 = D2ar.OK[iOK, , ]
  
  D1x = D1ar.OK[iOK, ,1 ]
  D1y = D1ar.OK[iOK, ,2 ]
  D2x = D2ar.OK[iOK, ,1 ]
  D2y =D2ar.OK[iOK, ,2 ]
  # Spline functions for D1 and D2
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
  x = seq( -3.1,xeq[2],length.out=10) # mesh of 20 (vs 30) does not blow up bvptwp
  
  # solve with bvpcol or bvptwp
  trycol <- tryCatch({
    bvptwp(yini = yini, x = x, func = feval2, yend = yend, parm=plist)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")} )
  
  # save solution for further analysis
  ETL = trycol # left exit time
  
  # end of solution for basin 1
  
  # Right basin 
  
  # set up for bvpSolve
  yini = c(0, NA)
  yend = c(NA, 0)
  
  # right basin from x=xeq[2] (absorbing) to x > xeq[3] (reflecting)
  x = seq(xeq[2],3.1,length.out=10) # mesh of 20 (vs 30) does not blow up bvptwp
  # solve with bvpcol or bvptwp
  trycol <- tryCatch({
    bvptwp(yini = yini, x = x, func = feval2, yend = yend, parm=plist)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")} )
  
  # save solution for further analysis
  ETR = trycol # right exit time
  
  # end of solution for basin 2
  
  # CALCULATE WEIGHTS FOR AVERAGE EXIT TIME -----------------------------------------------------  
  # Weight of x is p(x) from stationary distribution of Langevin equation  
  
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
  
  # sigma for stationary distrib is not estimated so set it to 1
  sigma.sd = 1
  
  # Calculate weighted average ET for both basins ===================================
  # ETL[N,1] is the same as ETR[1,1]; they connect at xeq[2]
  nL = length(ETL[,1])
  nR = length(ETR[,1])
  x = c(ETL[1,1]-0.01,ETL[2:nL,1],ETR[2:nR,1],ETR[nR,1]+0.01) # x has gaps that match ETL+ETR
  dx = diff(x)
  ETboth = c(ETL[1,2]-0.01,ETL[2:nL,2],ETR[2:nR,2],ETR[nR,2]+0.01)
  nx = length(x)
  
  # Weights by method in Science paper
  # See 2020-08-20 version of this code for the H&L version
  wtraw = rep(0,nx)
  for(i in 2:(nx)) {
    intpart = hcubature(f=finner,lowerLimit=x[1],upperLimit=x[i])$integral
    epart = exp(intpart)
    # normalization in Science paper
    wtraw[i] = gg(x[i])*epart 
  }
  
  # normalize weights
  wts = wtraw/sum(wtraw)
  
  # artistic convenience on first weight
  wts[1] = 0.98*wts[2]
  
  # Calculate mean exit time for left basin ---------------------------------------------
  
  #
  meanETl = sum(ETL[,2]*wts[1:nL])/sum(wts[1:nL])
  
  # save axes
  xL = x[1:nL]
  wtsL = wts[1:nL]
  
  # Calculate weighted average ET for right basin ===========================================================
  
  #
  meanETr = sum(ETR[,2]*wts[(nL+1):nx])/sum(wts[(nL+1):nx])
  
  ETLR[iOK,] = c(meanETl,meanETr) # output result
  tnow = Sys.time()
  etime = difftime(tnow,tstart,units='mins')
  print(c(iOK,ETLR[iOK,],etime),quote=F) 
  
  # Calculate left and right median survival times ********************************************************
  
  # Left survival
  # Make grid for ode.1D
  N = 1000
  xrange = range(ETL[,1],na.rm=T)
  xs = seq(from=xrange[1],to=xrange[2],length.out=N) # left basin
  dx = diff(range(xs))/N
  # Langevin components on mesh
  drif = D1fun(xs)
  dif = D2fun(xs)
  # Set up ODE
  Pars = c(drif=drif,
           dif=dif,
           dx=dx)
  
  # set the duration and intervals
  TMAX = 1000 # up to 5000 may be needed
  times = seq(from=0,to=TMAX,by=0.1)
  
  # set the initial conditions
  Sini = rep(1,N)  # intial S=1
  
  # Run the survival model for left side
  out = ode.1D(y=Sini,times=times,func=Bonehead.left,parms=Pars,nspec=1, N=N)
  # subset for large number of time steps
  tsub = times #seq(1,length(times),1000) # t values for subset
  tsub = tsub/12 # convert time axis to hours
  Smat = out[,2:(N+1)] # out[times,2:(N+1)]
  # Half-life
  # Left basin
  options("max.contour.segments"=200000)
  zhalf0 = contourLines(x=tsub,y=xs,z=Smat,nlevels=1,levels=0.5)
  # Right basin
  #zrevcol = apply(Smat,2,rev)  # flip z matrix on y axis for right basin only
  #zhalf0 = contourLines(x=out[,1],y=rev(xgrid$x.mid),z=zrevcol,nlevels=1,levels=zvals)
  zhalf1 = as.data.frame(zhalf0)
  # extract a subset
  nhalf = length(zhalf1$x)
  halfkeep = seq(1,nhalf,length.out=100)
  zhalf = zhalf1[halfkeep,]
  
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
  
  Shalf.left = sum(SxW)/sum(SW)
  ShalfLR[iOK,1] = Shalf.left
  print(c('Weighted Half-Life, Left ',Shalf.left),quote=F)
  
  # Right survival time ---------------------------------------------------------------
  
  # Make grid for ode.1D
  N = 1000
  xrange = range(ETR[,1],na.rm=T)
  xs = seq(from=xrange[1],to=xrange[2],length.out=N)
  dx = diff(range(xs))/N
  # Langevin components on mesh
  drif = D1fun(xs)
  dif = D2fun(xs)
  # Set up ODE
  Pars = c(drif=drif,
           dif=dif,
           dx=dx)
  
  # set the duration and intervals
  TMAX = 1000 # up to 5000 may be needed
  times = seq(from=0,to=TMAX,by=0.1)
  
  # set the initial conditions
  Sini = rep(1,N)  # intial S=1
  
  # Run the survival model for the right side 
  out = ode.1D(y=Sini,times=times,func=Bonehead.right,parms=Pars,nspec=1, N=N)
  # subset for large number of time steps
  tsub = times #seq(1,length(times),1000) # t values for subset
  tsub = tsub/12 # convert time axis to hours
  Smat = out[,2:(N+1)] # out[times,2:(N+1)]
  # Half-life
  # Left basin
  options("max.contour.segments"=200000)
  zhalf0 = contourLines(x=tsub,y=xs,z=Smat,nlevels=1,levels=0.5)
  # Right basin
  #zrevcol = apply(Smat,2,rev)  # flip z matrix on y axis for right basin only
  #zhalf0 = contourLines(x=out[,1],y=rev(xgrid$x.mid),z=zrevcol,nlevels=1,levels=zvals)
  zhalf1 = as.data.frame(zhalf0)
  # extract a subset
  nhalf = length(zhalf1$x)
  halfkeep = seq(1,nhalf,length.out=100)
  zhalf = zhalf1[halfkeep,]
  
  # compute weighted mean survival time
  # Write stationary density (weights) as a function 
  Wspline = smooth.spline(x=x,y=wts)
  Wfct = function(x) {
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
  
  Shalf.right = sum(SxW)/sum(SW)
  ShalfLR[iOK,2] = Shalf.right
  print(c('Weighted Half-Life, Right ',Shalf.right),quote=F)
  print('',quote=F)
}

tstop = Sys.time()
print('',quote=F)
print('----------------------------------------------------------',quote=F)
runtime = difftime(tstop,tstart,units='mins')
print(c('runtime, minutes ',runtime),quote=F)

# unpack result
save(ETLR,ShalfLR,file=fname)
