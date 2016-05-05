###################################################
### chunk number 1: load-library
###################################################
library(kernelPop)
set.seed(123) ###change this if you use this code in production!


###################################################
### chunk number 2: create-skeleton
###################################################
land <-  landscape.new.empty()
names(land)


###################################################
### chunk number 3: add-intparams
###################################################
args(landscape.new.intparam)
land <- landscape.new.intparam(land,h=4,s=6)


###################################################
### chunk number 4: add-floatparams
###################################################
args(landscape.new.floatparam)
land <- landscape.new.floatparam(land)


###################################################
### chunk number 5: add-boolean
###################################################
args(landscape.new.switchparam)
land <- landscape.new.switchparam(land)


###################################################
### chunk number 6: setup-local
###################################################
S <- matrix(c( 
              0,    0,  0,     0,    0,  0,
              0,    0,  0,     0,    0,  0,
              0.18, 0,  0,     0,    0,  0,
              0,  0.14,  0,  0.26,    0,  0,
              0,    0, 0.7,    0, 0.09,  0,
              0,    0,  0,   0.2,    0,  0.18), byrow=T,nrow=6)
              
R <- matrix(c(
              0,  0,  0,  0,  14,  0,
              0,  0,  0,  0,  8.5,  0,
              0,  0,  0,  0,  0,  0,
              0,  0,  0,  0,  0,  0,
              0,  0,  0,  0,  0,  0,
              0,  0,  0,  0,  0,  0), byrow=T,nrow=6)
              
M <- matrix(c(
              0,  0,  0,  0,  0,  0,
              0,  0,  0,  0,  0,  0,
              0,  0,  0,  0,  0,  0,
              0,  0,  0,  0,  0,  0,
              0,  0,  0,0.25, 0,  0.75,
              0,  0,  0,  0,  0,  0), byrow=T,nrow=6)
args(landscape.new.local.demo)
land <- landscape.new.local.demo(land, S, R, M)
print ("add a new local demography with different Reproduction")

R2 <- matrix(c(
              0,  0,  0,  0,  8,  0,
              0,  0,  0,  0,  5.5,  0,
              0,  0,  0,  0,  0,  0,
              0,  0,  0,  0,  0,  0,
              0,  0,  0,  0,  0,  0,
              0,  0,  0,  0,  0,  0), byrow=T,nrow=6)

land <- landscape.new.local.demo(land, S, R2, M)



###################################################
### chunk number 7: 
###################################################
zeromat <- matrix(0,nrow=4*6,ncol=4*6)


###################################################
### chunk number 8: setup-carry-exinct
###################################################
extnct <- c(0.0,0.1,0.0,0.1)
k <- c(1000,600,600,1000)


###################################################
### chunk number 9: setup-localdemography
###################################################
ldem <- c(0.5,0.5)


###################################################
### chunk number 10: seedkernel
###################################################
sk <- matrix(0,nrow=4*6,ncol=6)
sk[,1] <- rep(3,4*6)
sk[,2] <- rep(10,4*6)
sk[,3] <- rep(1.1,4*6)
sk[,4] <- rep(100,4*6)
sk[,5] <- rep(50,4*6)
sk[,6] <- rep(0.5,4*6)
sk


###################################################
### chunk number 11: pollenkernel
###################################################
pk <- matrix(0,nrow=4*6,ncol=6)
pk[,1] <- rep(3,4*6)
pk[,2] <- rep(5,4*6)
pk[,3] <- rep(2,4*6)
pk[,4] <- rep(100,4*6)
pk[,5] <- rep(50,4*6)
pk[,6] <- rep(1,4*6)
pk


###################################################
### chunk number 12: locations
###################################################
lx <- c(0,0,800,800)
rx <- c(600,600,1400,1400)
bty <- c(0,800,0,800)
ty <- c(600,1400,600,1400)


###################################################
### chunk number 13: make-epoch
###################################################
args(landscape.new.epoch)
land <- landscape.new.epoch(land,S=zeromat,R=zeromat,M=zeromat,
                            extinct=extnct,carry=k,localprob=ldem,
                            pollen.kernels=pk,seed.kernels=sk,
                            leftx=lx,rightx=rx,boty=bty,topy=ty)


###################################################
### chunk number 14: make-loci
###################################################
args(landscape.new.locus)
land <- landscape.new.locus(land,type=0,ploidy=1,transmission=1,numalleles=5)
land <- landscape.new.locus(land,type=1,ploidy=2,transmission=0,numalleles=3)
land <- landscape.new.locus(land,type=2,ploidy=2,transmission=0,numalleles=3,
                            allelesize=125)
length(land$loci)


###################################################
### chunk number 15: demo-some-functions
###################################################
landscape.ploidy(land)
landscape.democol()


###################################################
### chunk number 16: add-individuals
###################################################
vlen <- land$intparam$habitats*land$intparam$stages
vec <- rep(100,vlen)
land <- landscape.new.individuals(land,vec)


###################################################
### chunk number 17: Simulation
###################################################

l1 <- landscape.simulate(land,10)
l2 <- landscape.simulate(land,10)
l3 <- landscape.simulate(land,10)
l4 <- landscape.simulate(land,10)



###################################################
### chunk number 18: plot-landscape
###################################################
par(mfrow=c(2,2))
landscape.plot.locations(l1)
landscape.plot.locations(l2)
landscape.plot.locations(l3)
landscape.plot.locations(l4)
par(mfrow=c(1,1))


###################################################
### chunk number 19:  eval=FALSE
###################################################
## save(file="l1.rda",l1)
## rm(l1)
## load(file="l1.rda",l1)


###################################################
### chunk number 20: change-extinct
###################################################
names(l1$demography$epochs[[1]])

sk <- matrix(0,nrow=4*6,ncol=6)
sk[,1] <- rep(3,4*6)
sk[,2] <- rep(10,4*6)
sk[,3] <- rep(1.1,4*6)
sk[,4] <- rep(400,4*6)
sk[,5] <- rep(100,4*6)
sk[,6] <- rep(0.5,4*6)
sk

l1$demography$epochs[[1]]$seedkern <- sk



###################################################
### chunk number 21: simulate-again
###################################################
l1 <- landscape.simulate(l1,10)
l2 <- landscape.simulate(l2,10)
l3 <- landscape.simulate(l3,10)
l4 <- landscape.simulate(l4,10)


###################################################
### chunk number 22: plot-again
###################################################
par(mfrow=c(2,2))
landscape.plot.locations(l1)
landscape.plot.locations(l2)
landscape.plot.locations(l3)
landscape.plot.locations(l4)
par(mfrow=c(1,1))


###################################################
### chunk number 23: seed-kernel
###################################################
source("../test/distance-functions.R")
par(mfrow=c(2,2))
hist(seed.dist(l1),breaks=30,xlab="seed dispersal distance")
hist(seed.dist(l2),breaks=30,xlab="seed dispersal distance")
hist(seed.dist(l3),breaks=30,xlab="seed dispersal distance")
hist(seed.dist(l4),breaks=30,xlab="seed dispersal distance")
par(mfrow=c(1,1))


###################################################
### chunk number 24: pollen-kernel
###################################################
par(mfrow=c(2,2))
hist(pollination.dist(l1),breaks=30,xlab="pollination dispersal distance")
hist(pollination.dist(l2),breaks=30,xlab="pollination dispersal distance")
hist(pollination.dist(l3),breaks=30,xlab="pollination dispersal distance")
hist(pollination.dist(l4),breaks=30,xlab="pollination dispersal distance")
par(mfrow=c(1,1))


###################################################
### chunk number 25: popsize
###################################################
table(landscape.populations(l1))


###################################################
### chunk number 26: simulation3
###################################################
gland <- land
landlist <- vector("list",11)
landlist[[1]] <- gland
for(i in 2:11)
{
  print(table(landscape.populations(gland)))
  gland <- landscape.simulate(gland,10)
  landlist[[i]] <- gland
  if (i==6)
    {
      sk <- gland$demography$epochs[[1]]$seedkern
      sk[,4] <- rep(500,dim(sk)[1])
      sk[,6] <- rep(0.8,dim(sk)[1])
      gland$demography$epochs[[1]]$seedkern <- sk
    }
}


###################################################
### chunk number 27: summary-stat
###################################################
plot.ob <- do.call(rbind,lapply(landlist,function(l)
                                {
                                  c(l$intparam$currentgen,
                                    mean(landscape.amova(l,ns=24)))
                                  }))
print(xyplot(plot.ob[,2]~plot.ob[,1],type=c("b","smooth"),xlab="Time in years",ylab="Mean Phi-ST"))


###################################################
### chunk number 28: genepopout
###################################################
landscape.genepop.output(gland)


###################################################
### chunk number 29: arlequin-out
###################################################
l <- landscape.write.foreign(gland,fn="diploid-arlequin.arb",fmt="arlequin")
#l <- landscape.write.foreign(gland,fn="haploid-arlequin.arb",fmt="arlequinhap")


###################################################
### chunk number 30: migrate
###################################################
l <- landscape.write.foreign(gland,fn="migrate.infile",fmt="migrate")


###################################################
### chunk number 31: migrate
###################################################
l <- landscape.write.foreign(gland,fn="biosys.txt",fmt="biosys")


###################################################
### chunk number 32: fdist
###################################################
landscape.write.fdist(gland)


