##
## here is a function that returns the coordinates of each population in a landscape
## we now depend on dplyr
##
landscape.popcoords <- function(rland)
{
    as.data.frame(data.frame(pop=landscape.populations(rland),
               x=rland$individuals[,4],
               y=rland$individuals[,5]) %>%
        group_by(pop) %>% summarize(x=mean(x),y=mean(y)))
    
}

##makes a data frame of phenotypic means and variances for each population phenotype comparison
landscape.phenosummary  <- function(l)
{
    if ((is.landscape(l))&&(l$intparam$nphen>0))
    {
        phens <- as.data.frame(cbind(pop=landscape.populations(l),landscape.phenotypes.c(l)))
        names(phens)[-1] <- paste0("phen",1:l$intparam$nphen)
        phensum <- left_join(group_by(phens,pop) %>% summarise_all(.funs=c("mean","sd")),
                             group_by(phens,pop) %>% summarise(n=n()))
        data.frame(full_join(data.frame(pop=1:l$intparam$habitats),phensum)%>%arrange(pop))
    }  else {NULL}
}


###this function takes a list of items that describe calcuation of summary statistics (essntially stats functions)
landscape.gensummary <- function(l,analyses=list(He=list(format=1, #1 is landscape 2 is genetics::genotype
                                                         name="ExpHet",
                                                         func=function(l){
                                                             as.data.frame(cbind(pop=unique(landscape.populations(l)),landscape.exp.het(l)))
                                                         }
                                                         ),
                                                 LD=list(format=1, #1 is landscape 2 is genetics::genotype
                                                         name="LD",
                                                         func=function(l){
                                                             as.data.frame(landscape.LD(l))
                                                         }
                                                         ))
                                 )
{
    if (is.landscape(l))
    {
        lst=lapply(analyses, function(x)
        {
            if (x$format==1) res = x$func(l)
            else if (x$format==2) res = x$func(g)
            else stop("specify a correct input format for the statistic")
            res$stat=x$name
            res
        })
        
    } else {NULL}
}

###the genetics package genotypes
### returns a list of genotypes
landscape2genotype <- function(l)
{
    lapply(which(landscape.ploidy(l)==2),function(loc)
    {
        tloc <- matrix(landscape.locus(loc,l)[,-1:-9],ncol=2)
        tloc2=cbind(ifelse((tloc[,1]>tloc[,2]),tloc[,2],tloc[,1]),
                    ifelse((tloc[,1]>tloc[,2]),tloc[,1],tloc[,2]))
        tloc2f=factor(paste(tloc2[,1],tloc2[,2],sep="/"))
        attr(tloc2f,"allele.names")=as.character(unique(c(tloc)))
        attr(tloc2f,"allele.map")=do.call(rbind,strsplit(levels(tloc2f),"/"))
        attr(tloc2f,"genotypeOrder")=as.character(c("1/1","1/2","2/1","2/2"))
        class(tloc2f) <- "genotype"
        class(tloc2f) <- append(class(tloc2f),"factor")
        tloc2f
    })
}


landscape.LD <- function(l)
{

    ret <- NULL
    ##rsq is standardized with the overall allele freqs
    afreqs <- matrix(NA,nrow=sum(landscape.ploidy(l)==2),ncol=2)
    locs <- list()
    pops=landscape.populations(l)
    for (i in which(landscape.ploidy(l)==2))  ###be careful because could have haploid low-numbered loci (probably should just toss the haploids
    {
        locs[[i]] <- landscape.locus(i,l)[,-1:-9]
        tbl= table(locs[[i]][,1])
        afreqs[i,] <- tbl/sum(tbl)
    }

    ret <- do.call(rbind,lapply(which(landscape.ploidy(l)==2),function(i)
    {
        rdf=NULL
        for (j in 1:i)
            if (j<i)
                if (((afreqs[i,1]>0)&(afreqs[i,1]<1))&((afreqs[j,1]>0)&(afreqs[j,1]<1))) #polymorphism
                {
                    df <- data.frame(pop=pops,
                                     av=locs[[i]][,1],
                                     ev=locs[[i]][,1]==locs[[j]][,1],
                                     gv=locs[[i]][,1]>locs[[j]][,1]
                                     )
                    mns <- df%>%mutate(p11=(av==1)&(ev&(!gv)),
                                       p22=(av==2)&(ev&(!gv)),
                                       p12=((!ev)&(!gv)),
                                       p21=((!ev)&(gv))) %>%
                        group_by(pop)%>%summarise(p11=mean(p11),p12=mean(p12),p21=mean(p21),p22=mean(p22))
                    
                    rdf <- rbind(rdf,mutate(mns,
                                            D=(p11*p22)-(p12*p21),
                                            rsq= (((p11*p22)-(p12*p21))/sqrt(afreqs[i,1]*afreqs[j,1]*afreqs[i,2]*afreqs[j,2]))^2,
                                            loc1=i,
                                            loc2=j)
                                 )
                }
        rdf
    }))
    ret[apply(ret[,2:5],1,function(v)max(v)<1),]
}




landscape.LD.old <- function(l)
{

    afreqs <- matrix(NA,nrow=sum(landscape.ploidy(l)==2),ncol=2)
    ret <- NULL
    for (i in which(landscape.ploidy(l)==2))
    {
        loc1 = landscape.locus(i,l)[,-1:-9]
        if (length(unique(c(landscape.locus(i,l)[,-1:-9])))>1)
        {
            if (is.na(afreqs[i,1]))
            {
                tbl= table(loc1[,1])
                afreqs[i,] <- tbl/sum(tbl)
            }
            
            for (j in which(landscape.ploidy(l)==2))
                if (j<i)
                {
                    loc2 = landscape.locus(j,l)[,-1:-9]
                    if (length(unique(c(loc2)))>1)
                    {
                 
                        if (is.na(afreqs[j,1]))
                        {
                            tbl= table(loc2[,1])
                            afreqs[j,] <- tbl/sum(tbl)
                        }

#                        genos=unique(cbind(loc1[,1],loc2[,1]))
                        tbl=table(paste(loc1[,1],loc2[,1],sep="/"))

                         if (sum(names(tbl)=="1/1")>0) p11=tbl[names(tbl)=="1/1"]/sum(tbl) else p11 = 0
                        if (sum(names(tbl)=="1/2")>0) p12=tbl[names(tbl)=="1/2"]/sum(tbl) else p12=0
                        if (sum(names(tbl)=="2/1")>0) p21=tbl[names(tbl)=="2/1"]/sum(tbl) else p21=0
                         if (sum(names(tbl)=="2/2")>0) p22=tbl[names(tbl)=="2/2"]/sum(tbl) else p22=0
                        print(c(p11,p12,p21,p22))
                        atbl1 = table(loc1[,1])
                        atbl2 = table(loc2[,1])
                        p1=afreqs[i,1]
                        p2=afreqs[i,2]
                        q1=afreqs[j,1]
                        q2=afreqs[j,2]
                        D=p11*p22 - p12*p21
                        rsq = (D/sqrt(p1*p2*q1*q2))^2
                        print(c(i,j,D,rsq))
                        ret <- rbind(ret,cbind(loc1=i,loc2=j,D=D,rsq=rsq))
                    }
                }
            }
    }
    ret
    
          
}


###
###
###
landscape.dist2origin <- function(l,origin=NULL)
{
    if (is.null(origin)) {origin <- sample(unique(landscape.populations(l)),1)}
    crd=data.frame(landscape.popcoords(l))
    rownames(crd) <- crd$pop
    d <- as.matrix(dist(crd[,c("x","y")],diag=T,upper=T))
    d[rownames(d)==as.character(origin),]
}

landscape.neighborhood <- function(rland)
{
    mo  <- landscape.geodist(rland,cmp=c("offspring","father"))
    fo  <- landscape.geodist(rland,cmp=c("offspring","father"))
    mp2 <- ((mo+fo)/2)^2
    df  <- data.frame(pop=landscape.populations(rland),mp2=mp2) %>%
        group_by(pop)%>%summarise(sigma2 = mean(mp2,na.rm=T))
    area=((rland$demography$epochs[[1]]$rightx-rland$demography$epochs[[1]]$leftx) *
          (rland$demography$epochs[[1]]$topy-rland$demography$epochs[[1]]$boty))[as.numeric(as.character(df$pop))]
    df <- as.data.frame(df)
    df$D =  rland$demography$epochs[[1]]$Carry[as.numeric(as.character(df$pop))]/area
    
    df$Nn = c(4*pi*df[,"sigma2"]*df$D)
    df[,c("pop","Nn")]
}
