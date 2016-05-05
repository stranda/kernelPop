##the landscape produced by this function is inteded to be an example of how different stages within populations can have different
##seed dispersal parameters (can also differ among populations too)
## this is set up as a demostration, using landscape.simulate will not work, but landscape.reproduce will produce the seeds and
##seed dispersal shadow
landscape.stage.kernel.test<-function (sdkernel = 1, polkernel=1, polscale=c(90,600), polshape=c(1,200), polmix=1,
                                 offspring=200)
{
    rland <- NULL
    rland <- landscape.new.empty()
    rland <- landscape.new.intparam(rland, h = 1, s = 6, np = 0, 
        totgen = 20000)
    rland <- landscape.new.switchparam(rland, mp = 0)
    rland <- landscape.new.floatparam(rland, s = 0, pollenscale = polscale,
                                      pollenshape = polshape, pollenmix=polmix,
                                      asp = 1)

    ##all non-seed stages survive forever.
    ##seeds (stage 0 die every year)
    S <- matrix(c(0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0,
                  0, 0, 0, 1, 0, 0,
                  0, 0, 0, 0, 1, 0,
                  0, 0, 0, 0, 0, 1
                  ), byrow = T, nrow = rland$intparam$stages)
    
    ##lower number stages produce less offspring
    ##
    R <- matrix(c(0, 0, 0, offspring, 0, 0,
                  0, 0, 0, 0, offspring, 0,
                  0, 0, 0, 0, 0, offspring,
                  0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0
                  ), byrow = T, nrow = rland$intparam$stages)
    
    M <- matrix(c(0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0,
                  0, 0, 0, 1, 1, 1,
                  0, 0, 0, 1, 1, 1,
                  0, 0, 0, 1, 1, 1
                  ), byrow = T, nrow = rland$intparam$stages)
    
    rland <- landscape.new.local.demo(rland, S, R, M)
    S <- matrix(rep(0, (rland$intparam$stages * rland$intparam$habitat)^2), nrow = rland$intparam$stages * 
        rland$intparam$habitat)
    R <- matrix(rep(0, (rland$intparam$stages * rland$intparam$habitat)^2), nrow = rland$intparam$stages * 
        rland$intparam$habitat)
    M <- matrix(rep(0, (rland$intparam$stages * rland$intparam$habitat)^2), nrow = rland$intparam$stages * 
        rland$intparam$habitat)

   
    if ((dim(sdkernel)[1]==dim(S)[1])&&(dim(sdkernel)[2]==6))
      {
        rland <- landscape.new.epoch(rland, S = S, R = R, M = M, 
                                     carry = 20000,
                                     extinct = rep(0.05, rland$intparam$habitat),
                                     seed.kernels = sdkernel,
                                     leftx = 0, rightx = 4000, boty = 0, topy = 4000,
                                     maxland = c(0,0,4000,4000))
      }
    else
      {
        stop("dimensions of the kernel object does not correspond to the landscape")
      }

    rland <- landscape.new.locus(rland,type=0, ploidy=1, mutationrate=0.005, numalleles=3, frequencies = c(0.2, 0.2, 0.6))
    rland <- landscape.new.locus(rland,type=1,ploidy=2,mutationrate=0.001,transmission=0,numalleles=5)
    ##print(unlist(rland$floatparam))
    rland <- landscape.new.individuals(rland, c(0,0,0,1,1,1))
    rland$individuals[1,c(4:5)] <- c(1000,2000)
    rland$individuals[2,c(4:5)] <- c(2000,2000)
    rland$individuals[3,c(4:5)] <- c(3000,2000)
    rland
}

###tries out the different kernels.  Plots offspring and mothers
###also compares the dispersal distances to what we intended
landscape.stage.kernel.demo <- function()
{
  ##this is the object that describe seed kernels for different stages in a landscape. This expects 4 stages
  ##with the first stage non-reproducing.  I'm using an exponential, weibull and mixed distribution for the four stages
  ##column 1 of this matrix.  The other columns are scale and shape parameters (two sets because some distributions are mixtures)
  ##and a mixting parameter in the last col.
  
  sk <- matrix(c(
                 0,   0, 0,   0,   0,    0,
                 0,   0, 0,   0,   0,    0,
                 0,   0, 0,   0,   0,    0,
                 3, 100, 1,   0,   0,    1,    #exponential with mean 100
                 3, 215.87, 2,   0,   0,    1,    #weibull with scale 200 and shape 2   
                 3, 10    , 1, 130, 125, 0.5     #mixture with weibull scale 50, shape 2,
                                                   #normal mean 300, sd 200, mixture 75% weibull
                 ), ncol=6,byrow=T)
  l <- landscape.stage.kernel.test(sdkernel=sk)
  l <- landscape.reproduce(l)


  # the rest of this function is just fancy plotting
  par(bg = "white")
  split.screen(c(2,1))        # split display into two screens
  split.screen(c(1,3), screen = 1) # now split the bottom half into 3
  screen(2)
  
  xmi <- min(c(l$ind[l$ind[,1]<3,4]))
  xma <- max(c(l$ind[l$ind[,1]<3,4]))
    
  ymi <- min(c(l$ind[l$ind[,1]<3,5]))
  yma <- max(c(l$ind[l$ind[,1]<3,5]))
  par(mar=c(5,4,0,1)+0.1)
  plot(l$ind[l$ind[,1]<3,5]~l$ind[l$ind[,1]<3,4],ylab="Y-coordinate",xlab="X-coordinate",pch=l$ind[,1]+1,ylim=c(ymi,yma),xlim=c(300,3700))
  text(c(1300,2350,3350),c(2400,2400,2400),c("Exponential","Weibull","Mixed Weibull-Gaussian"))
  par(mar=c(3,4,0,1)+0.1)

  screen(3)
  l1 <- l
  l1$individuals <- l1$individuals[l1$individuals[,1]==0,]
  hist(seed.dist(l1),xlab="",main="",freq=F,ylim=c(0,0.012),xlim=c(0,600))
  points(0:600,dweibull(0:600,scale=100,shape=1),type="l")
  text(0.011~400,"Exponential")
  text(0.009~400,paste("Actual variance:",round(var(seed.dist(l1)),0)))
  screen(4)
  l1 <- l
  l1$individuals <- l1$individuals[l1$individuals[,1]==1,]
  hist(seed.dist(l1),xlab="dispersal distance from mother",main="",freq=F,ylim=c(0,0.012),xlim=c(0,600))
  points(0:600,dweibull(0:600,scale=215.87,shape=2),type="l")
  text(0.011~400,"Weibull")
  text(0.009~400,paste("Actual variance:",round(var(seed.dist(l1)),0)))  
  screen(5)
  l1 <- l
  l1$individuals <- l1$individuals[l1$individuals[,1]==2,]
  hist(seed.dist(l1),xlab="",main="",freq=F,ylim=c(0,0.012),xlim=c(0,600))
  points(0:600,dmixed(0:600,10,1,130,125,0.5),type="l")
  text(0.011~400,"Mixed Weibull-Gaussian")
  text(0.009~400,paste("Actual variance:",round(var(seed.dist(l1)),0)))  
  close.screen(all=T)
  par(mar=c(5,4,4,1)+0.1)

#  ind.kde <- kde2d(l$ind[l$ind[,1]==0,4],l$ind[l$ind[,1]==0,5],lims=c(xmi,xma,xmi,xma),h=rep(800,2),n=50)
#  contour(ind.kde,nlevels=20)
#  persp(ind.kde,phi = 30, theta = 20, d = 5)
#  par(mfrow=c(1,1))
}


#
#
# here the functions in this file are invoked
#

library(kernelPop)
landscape.stage.kernel.demo()
