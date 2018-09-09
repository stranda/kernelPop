#probably needs a good bit of error checking written in
                                        #

###aindex is a vector of the aindices for each locus that confer an additive effect
###hsq is a vector of heritabilities (now always 1)
landscape.new.expression <-
function(rland,expmat,addstates=NULL,hsq=NULL)  
{
    if (is.null(hsq)) hsq <- rep(1,dim(expmat)[2])  #change if we change the heritabilities for each trait
                                        #if diff than 0, used to add noise
    
    if (is.null(addstates)) addstates <- rep(1,
                                           length=dim(expmat)[1])
    
    if (dim(expmat)[1]!=length(rland$loci)) {stop("expression mat should have nloc rows and nphen columns")}
    rland$expression <- list(expmat=expmat,addstates=addstates,hsq=hsq)
    ##set the nphen in intparam
    rland$intparam$nphen <- dim(expmat)[2]
    
  rland
}



landscape.new.gpmap <- function(rland,gpdisp=NULL,gpdemo=NULL)
{
    if ((rland$intparam$nphen<1)||(is.null(rland$expression)))
    {
        rland$gpmap= list(gpdisp=cbind(rep(-1,5),rep(0,5),rep(1,5)),
                          gpdemo=cbind(rep(-1,3),rep(0,3),rep(1,3)))
    } else {

        if ((is.null(gpdisp))||(!is.matrix(gpdisp)))
        {
            gpdisp=cbind(rep(-1,5),rep(0,5),rep(1,5))
        }

        if ((is.null(gpdemo))||(!is.matrix(gpdemo)))
        {
            gpdemo=gpdemo=cbind(rep(-1,3),rep(0,3),rep(1,3))
        }
        
        rland$gpmap <- list(gpdisp=gpdisp,gpdemo=gpdemo)
        
    }
    rland
}
