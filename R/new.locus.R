"landscape.new.locus" <-
function(rland,type=0,ploidy=1,mutationrate=0,transmission=1,numalleles=2,frequencies=NULL,states=NULL)
{
  if (!(is.list(rland$loci)))
    {
      rland$loci <- list(list(type=0,ploidy=0,trans=0,rate=0,alleles=0))
      locusnum <- 1
    }
  else
    {
      locusnum <- length(rland$loci) + 1
      rland$loci[[locusnum]] <- list(type=0,ploidy=0,rate=0,trans=0,alleles=0)
    }

  rland$intparam$locusnum <- locusnum
  
  if(type >= 0 && type <= 1)
    {
      rland$loci[[locusnum]]$type <- typelookup(type)
    }
  else
    {
      stop("Invalid type of locus")
    }

  if(ploidy == 1 || ploidy == 2)
    {
      rland$loci[[locusnum]]$ploidy <- ploidy
    }
  else
    {
      stop("Invalid ploidy count")
    }

  rland$loci[[locusnum]]$rate <- mutationrate

  if(transmission == 0 || transmission == 1)
    {
      rland$loci[[locusnum]]$trans <- as.integer(transmission)
    }
  else
    {
      stop("Invalid transmission number")
    }

  if(numalleles >= 0)
    {
      rland$loci[[locusnum]]$alleles <- makealleles(type,numalleles,allelesize,frequencies,states)
    }
  else
    {
      stop("Need non-negative numbers of alleles")
    }

  rland
}

