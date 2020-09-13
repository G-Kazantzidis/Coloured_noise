#' Creates a trajectory using coloured noise
#' @export
#' @param nu numeric variable (type 0 for white, 1 for pink or 2 for brown)
#' @param tmax numeric variable (length of the trajectory)

noise.f<-function(nu,tmax)
{
  tmax<-tmax+1
  # Parameters
  beta<-4 # Base for temporal spacing
  cutoffM<-8 # Upper cutoff
  cutoffm<-0 # Lower cutoff

  # Initialisation
  k<- cutoffm:cutoffM
  lags<-length(k)
  t<-1:tmax
  Dphi <- log(beta)

  tau<-rep(NA,lags)
  ro<-rep(NA,lags)
  w<-rep(NA,lags)
  W<-rep(NA,lags)
  a<-matrix(NA,lags,tmax)
  aw<-matrix(NA,lags,tmax)
  n<-matrix(NA,lags,tmax)

  for(ki in 1:lags)
  {
    tau[ki]<-beta^k[ki]
    ro[ki]<-exp(-1/tau[ki])
    w[ki]<- exp((nu-1)*ki*Dphi)
  }


  for(ki in 1:lags)
  {
    a[ki,1]<-rnorm(1,0,0.01) # the variance set to 0.01(?)
    W[ki]<-w[ki]/sum(w)

    for(t in 1:(tmax-1))
    {
      n[ki,t]<-rnorm(1,0,0.01) # the variance set to 0.01(?)
      a[ki,t+1]<-ro[ki]*a[ki,t]+sqrt(1-ro[ki]^2)*n[ki,t]
      aw[ki,t]<-a[ki,t]*W[ki]
    }

  }

  noise<-rep(NA,tmax)

  for (i in 1:tmax)
  {noise[i]<-sum(aw[1:lags,i])}
  return(noise[1:(tmax-1)]) #The noise we want
}


#' Calculates the p value of a timeseries according to the given colour
#' @export
#' @param nu numeric variable (0 for white, 1 for pink or 2 for brown)
#' @param reps numeric variable (number of replications)
#' @param p double (data, including NAs if any)

Col_sign<-function(p,nu,reps)

{
  li<-data.frame(p,c(1:length(p)))
  li<-na.omit(li)
  reg<-summary(lm(li[,1]~li[,2]))
  tmax<-tail(li[,2],n=1)
  slop<-reg$coefficients[2,1]
  slope_for_test<-abs(reg$coefficients[2,1])
  pv<-rep(NA,(reps-1))

  for(k in 1:(reps-1))
  {
    tran<- as.numeric(noise.f(nu = nu ,tmax = tmax))
    tran<-tran[1:(tmax-1)] #  the last value is always NA..
    tran<-tran-mean(tran)#substruct the mean
    tran<-tran/sd(tran) # divide with StDev
    tran<-tran*sd(li[,1])#multiply with stdev of the population
    tran<-tran+mean(li[,1])# add the mean of the population
    form<-summary(lm(tran~c(1:(tmax-1)))) #find the slope of noise
    pv[k]<-abs(form$coefficients[2,1]) # the absolute slope
  }
  Coloured.m<- (sum(pv>slope_for_test)+1)/(2*reps)
  return(Coloured.m)
}
