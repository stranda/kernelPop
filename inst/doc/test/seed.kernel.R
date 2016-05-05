landscape.new.kernel.test<-function (sdkernel = 1, sdscale=c(90,600), sdshape=c(1,200), sdmix=1,
                                 polkernel=1, polscale=c(90,600), polshape=c(1,200), polmix=1,
                                 offspring=1000)
{
    rland <- NULL
    rland <- landscape.new.empty()
    rland <- landscape.new.intparam(rland, h = 1, s = 2, np = 0, 
        totgen = 20000)
    rland <- landscape.new.switchparam(rland, mp = 0)
    rland <- landscape.new.floatparam(rland, s = 0, seedscale = sdscale, seedshape = sdshape,
                                      seedmix=sdmix, pollenscale = polscale,
                                      pollenshape = polshape, pollenmix=polmix,
                                      asp = 1)

    S <- matrix(c(0, 0, 1.03*1/offspring, 0), byrow = T, nrow = 2)
    R <- matrix(c(0, offspring, 0, 0), byrow = T, nrow = 2)
    M <- matrix(c(0, 0, 0, 1), byrow = T, nrow = 2)
    rland <- landscape.new.local.demo(rland, S, R, M)
    S <- matrix(rep(0, (2 * rland$intparam$habitat)^2), nrow = 2 * 
        rland$intparam$habitat)
    R <- matrix(rep(0, (2 * rland$intparam$habitat)^2), nrow = 2 * 
        rland$intparam$habitat)
    M <- matrix(rep(0, (2 * rland$intparam$habitat)^2), nrow = 2 * 
        rland$intparam$habitat)

    rland <- landscape.new.epoch(rland, S = S, R = R, M = M, 
                                 carry = 20000,
                                 extinct = rep(0.05, rland$intparam$habitat),
                                 leftx = 0, rightx = 4000, boty = 0, topy = 4000,
                                 maxland = c(0,0,4000,4000))

    rland$demography$epochs[[1]]$seedkern[,1] <- rep(sdkernel,rland$intparam$habitats*rland$intparam$stages)

    rland <- landscape.new.locus(rland,type=0, ploidy=1, mutationrate=0.005, numalleles=3, frequencies = c(0.2, 0.2, 0.6))
    rland <- landscape.new.locus(rland,type=1,ploidy=2,mutationrate=0.001,transmission=0,numalleles=5)
#print(unlist(rland$floatparam))
    rland <- landscape.new.individuals(rland, c(0,1))
    rland$individuals[1,c(4:5)] <- c(1500,1500)
    rland
}

###tries out the different kernels.  Plots offspring and mothers
###also compares the dispersal distances to what we intended
landscape.kernel.demo <- function()
{
##  variances for all distributions held at 10,000
  ssc1 <- c(100,120,50) #exp mean, weibull scale, scale of weibull component of mixture
  ssh1 <- c(0,1.2,1) # exp scale doesn't matter, weibull scale adjusted to give variance of 100
  ssc2 <- c(0,0,400)
  ssh2 <- c(0,0,250)
  smix <- c(0,0,0.94)

  knames <- c("Exponential","Weibull","Mixed")
  
  psc1 <- c(50,50,50)
  psh1 <- c(0,2,0)
  psc2 <- c(0,0,500)
  psh2 <- c(0,0,50)
  pmix <- rep(0.5,3)

  par(mfrow=c(3,1))
  for (k in 1:3)
    {
      l <- landscape.new.kernel.test(sdkernel=k, sdscale=c(ssc1[k],ssc2[k]), sdshape=c(ssh1[k],ssh2[k]),
                                     sdmix=smix[k])
      l <- landscape.reproduce(l)
      dists <- sqrt((l$ind[l$ind[,6]>0,4]-l$ind[l$ind[,6]>0,6])^2+(l$ind[l$ind[,6]>0,5]-l$ind[l$ind[,6]>0,7])^2)
      hist(dists,freq=F,breaks=25,xlab="Seed dispersal distance from mother",main=paste(knames[k],"kernel"),
           xlim=c(0,900))
      
      pts <- seq(floor(min(dists)),ceiling(max(dists)),length=100)

      if (k==1)
        {
          print(mean(dists))
          print(var(dists))
          points(pts,dweibull(pts,scale=ssc1[k],shape=1),type="l",lwd=0.5)
        }
      if (k==2)
        {
          print(mean(dists))
          print(var(dists))
          points(pts,dweibull(pts,scale=ssc1[k],shape=ssh1[k]),type="l",lwd=0.5)
        }
      if (k==3)
        {
          print(mean(dists))
          print(var(dists))
          points(pts,(smix[k]*dweibull(pts,scale=ssc1[k],shape=ssh1[k])+
                              (1-smix[k])*dnorm(pts,mean=ssc2[k],sd=ssh2[k])/(1-pnorm(0,mean=ssc2[k],sd=ssh2[k]))),type="l",lwd=0.5)
        }
    }
  par(mfrow=c(1,1))
}


#
#
# here the functions in this file are invoked
#
library(kernelPop)

