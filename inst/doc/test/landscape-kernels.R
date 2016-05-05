source("mixed-dist-prob-functions.R")

##
#just gives the densities in up to the first four kernels in the local demography
pollen.density <- function(l)
  {
    pk <- l$demography$epochs[[1]]$pollenkern
    pollenprod <- which(apply(l$demography$localdem[[1]]$LocalM,1,sum)>0)
    pkplot=pk[pollenprod,]
    pkplot=cbind(pkplot,pollenprod)
    pkplot=pkplot[1:min(dim(pkplot)[1],4),]
    par(mfrow=c(dim(pkplot)[1],1))
    apply(pkplot,1,function(x,mindens){
      if (x[1]==1)
        {
          xmax=max(which((dweibull(0:10000,scale=x[2],shape=1)-mindens)>0))
          plot(0:xmax,dweibull(0:xmax,scale=x[2],shape=1)-mindens,type="n",lwd=1.5,
               main=paste("Stage:",x[7]),xlab="Pollen dispersal distance",
               ylab="Density",
               ylim=c(0,max(dweibull(0:xmax,scale=x[2],shape=1)-mindens)))
          points(0:xmax,ifelse(dweibull(0:xmax,scale=x[2],shape=1)-mindens>0,dweibull(0:xmax,scale=x[2],shape=1)-mindens,NA),type="l")
        }
      if (x[1]==2)
        {
          xmax=max(which((dweibull(0:10000,scale=x[2],shape=x[3])-mindens)>0))
          plot(0:xmax,dweibull(0:xmax,scale=x[2],shape=x[3])-mindens,type="l",lwd=1.5,main=paste("Stage:",x[7]),xlab="Pollen dispersal distance",
               ylab="Density",ylim=c(0,max(dweibull(0:xmax,scale=x[2],shape=x[3])-mindens)))
          points(0:xmax,ifelse(dweibull(0:xmax,scale=x[2],shape=x[3])-mindens>0,dweibull(0:xmax,scale=x[2],shape=x[3])-mindens,NA),type="l")
        }
      if (x[1]==3)
        {
          xmax=max(which((dmixed(0:10000,x[2],x[3],x[4],x[5],x[6])-mindens)>0))
          plot(0:xmax,dmixed(0:xmax,x[2],x[3],x[4],x[5],x[6])-mindens,type="n",lwd=1.5,main=paste("Stage:",x[7]),xlab="Pollen dispersal distance",
               ylab="Density",ylim=c(0,max(dmixed(0:xmax,x[2],x[3],x[4],x[5],x[6])-mindens)))
          points(0:xmax,ifelse(dmixed(0:xmax,x[2],x[3],x[4],x[5],x[6])-mindens>0,dmixed(0:xmax,x[2],x[3],x[4],x[5],x[6])-mindens,NA),type="l")
        }
    },mindens=l$floatparam$mindens)
    
    par(mfrow=c(1,1))
  }


seed.density <- function(l)
  {
    pk <- l$demography$epochs[[1]]$seedkern
    seedprod <- which(apply(l$demography$localdem[[1]]$LocalR,2,sum)>0)
    pkplot=pk[seedprod,]
    pkplot=cbind(pkplot,seedprod)
    pkplot=pkplot[1:min(dim(pkplot)[1],4),]
    par(mfrow=c(dim(pkplot)[1],1))
    apply(pkplot,1,function(x){
      xmax=1000
      if (x[1]==1)
        {
          plot(0:xmax,dweibull(0:xmax,scale=x[2],shape=1),type="l",lwd=1.5,
               main=paste("Stage:",x[7]),xlab="Seed dispersal distance",
               ylab="Density",
               ylim=c(0,max(dweibull(0:xmax,scale=x[2],shape=1))))
        }
      if (x[1]==2)
        {
          plot(0:xmax,dweibull(0:xmax,scale=x[2],shape=x[3]),type="l",lwd=1.5,main=paste("Stage:",x[7]),xlab="Seed dispersal distance",
               ylab="Density",ylim=c(0,max(dweibull(0:xmax,scale=x[2],shape=x[3]))))
        }
      if (x[1]==3)
        {
          plot(0:xmax,dmixed(0:xmax,x[2],x[3],x[4],x[5],x[6]),type="l",lwd=1.5,main=paste("Stage:",x[7]),xlab="Pollen dispersal distance",
               ylab="Density",ylim=c(0,max(dmixed(0:xmax,x[2],x[3],x[4],x[5],x[6]))))
        }
    })
    
    par(mfrow=c(1,1))
  }
