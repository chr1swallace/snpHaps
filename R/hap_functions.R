library(parallel)
library(snpStatsWriter)
library(mice)
library(magrittr)

fix.alleles <- function(long,short) {
  strsplit(long,"") %>%
    mapply(setdiff, . , short) %>%
      lapply(.,paste0,collape="") %>%
        unlist()
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @return string with location of snphap executable
##' @author Chris Wallaced
snphap.exec <- function() {
    system.file("bin/snphap",package="snpHaps")
}

##' Generate haplotypes 
##'
##' Write a snpStats object to file, and use snphap to phase with multiple imputation
##' @param data snpStats object
##' @param slist snp list
##' @param snps snps dataframe with columns "allele.A", "allele.B" denoting alleles
##' @param samples samples dataframe
##' @param A1 colname of allele 1 in snps
##' @param A2 colname of allele 2 in snps
##' @param cols.copy columns to be copied to resulting phased object
##' @param f.in filename for snphap input (will be overwritten)
##' @param redo set to TRUE to regenerate snphap output regardless of whether it exists
##' @return data.frame containing the best phased haplotypes
##' @export
##' @importFrom snpStatsWriter write.snphap
##' @author Chris Wallace
genhaps <- function(data, slist=colnames(data), snps, samples,A1="allele.A",A2="allele.B",cols.copy=character(0),
                    f.in="snphap/haps.in", redo=FALSE) {
  if(length(cols.copy))
    cols.copy <- intersect(cols.copy,colnames(samples))
  if(!file.exists(dirname(f.in)))
    dir.create(dirname(f.in),recursive=TRUE)
  ## any simple corrections to multi-base SNPs needed?
  a1=snps[slist,A1]
  a2=snps[slist,A2]
  if(is.factor(a1))
      a1 <- as.character(a1)
  if(is.factor(a2))
      a2 <- as.character(a2)
  if(any(nchar(a1)!=1) || any(nchar(a2)!=1)) { # can we fix it?
    w <- which(nchar(a1)>1 & nchar(a2)==1)
    if(length(w))
        for(i in w) 
            a1[i] <- fix.alleles(a1[i],a2[i])[1]
    w <- which(nchar(a2)>1 & nchar(a1)==1)
    if(length(w))
        for(i in w) 
            a2[i] <- fix.alleles(a2[i],a1[i])[1]
  }
  message("writing snphap input file ",f.in)
  write.snphap(data[,slist],
               a1=a1,
               a2=a2,
               file=f.in)
  f.out1 <- paste0(sub(".in$","",f.in),".out1")
  f.out2 <- paste0(sub(".in$","",f.in),".out2")
  if(redo || !file.exists(f.out2) || !file.exists(f.out1)) {
    message("running snphap (this may take a while)")
    command <- paste(snphap.exec(), "-q -nh -mi 10 -ss",f.in,f.out1,f.out2," > kk")
    message(command)
    system(command)
  } else {
    message("Using existing snphap output in ",dirname(f.out1),". Run with redo=TRUE to regenerate.")
  }
  message("reading output")
  d <- hap.read(f.out2)
  if("cc" %in% colnames(samples))
    d$cc <- samples[d$id,"cc"]
  if(length(cols.copy))
    d <- cbind(d, samples[d$id,cols.copy,drop=FALSE])
  return(d)
}

hap.read <- function(f) {
  d <- read.table(f,sep="\t",as.is=TRUE,header=TRUE)
  d <- d[,-c(2,3)]
  d$hap <- apply(d[,-c(1,ncol(d)),drop=FALSE],1,paste,collapse="")
  d
}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param d
##' @param q
##' @param thr
##' @param base
##' @return
##' @export
##' @author Chris Wallace
model <- function(d,q=NULL,thr=0.001,base=NULL) {
  tt <- tapply(d$Prob,d$hap,sum)
  tt <- 100*tt/sum(tt)
##  print(100*tt/sum(tt))
  use <- names(tt)[tt/sum(tt)>thr]
  if(!is.null(base)) {
      if(!(base %in% names(tt))) {
          message("haplotypes found:")
          cat(names(tt),sep="\n")
          stop("requested base haplotype not found: ", base)
      }
  } else {
      base <- names(tt)[which.max(tt)]
  }
  d$hap <- relevel(factor(d$hap), ref=base)
  if("cc" %in% colnames(d)) {
    f <- "cc ~ hap"
    fam <- "binomial"
  } else {
    f <- "y ~ hap"
    fam <- "gaussian"
  }
  if(!is.null(q)) {
    d$q <- q
    f <- paste(f,"+ q")
  }
  m <- glm(as.formula(f), data=d, family=fam,subset=hap %in% use)
  null.value=0
  names(m$coefficients)[1] <- base
  ss <- data.frame(row.names=names(m$coefficients),
                   beta=c(null.value,m$coefficients[-1]),
                   se=c(NA,sqrt(diag(vcov(m)))[-1]),
                   p=c(NA,2*pnorm(abs(m$coefficients)/sqrt(diag(vcov(m))),lower=FALSE)[-1]))
  ss <- ss[,c("beta","se","p")]
  ss$Fq <- tt[sub("hap","",rownames(ss))]
  o <- order(ss$Fq[-1],decreasing=TRUE) + 1
  ss <- ss[c(1,o),]
  invisible(list(result=ss,BIC=BIC(m)))                   
}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param dir
##' @param df
##' @param thr
##' @param maxi
##' @param haps.pattern
##' @param covars
##' @param phenotype
##' @return 
##' @export
##' @author Chris Wallace
hapfreq <- function(dir,df,thr=0.001,maxi=NULL,haps.pattern="all-haps.out",
                     covars="q", phenotype="y") {
  haps <- hap.read(file.path(dir,haps.pattern))
  tt <- tapply(haps$Prob,haps$hap,sum)
  tt <- 100*tt/sum(tt)
  use <- names(tt)[tt/sum(tt)>thr]
  base <- names(tt)[which.max(tt)]
  lev <- names(sort(tt[use],decreasing=TRUE))

  ## read MI and average freqs
  files <- list.files(dir,pattern=paste0(haps.pattern,"."),full=TRUE)
  if(!is.null(maxi) && length(files)>maxi)
    files <- files[1:maxi]
  message("reading first multiply imputed dataset from ",dir)
 # fits <- lapply(seq_along(files), function(i) {
    cat(".")
    haps <- read.table(files[[1]],header=TRUE,sep="\t",as.is=TRUE)
    haps <- haps[,-c(2,3)]
    haps$hap <- apply(haps[,-c(1,ncol(haps)),drop=FALSE],1,paste,collapse="")
    haps <- haps[ haps$id %in% rownames(df), ,drop=FALSE]
    m <- match(haps$id, rownames(df))
    if(length(covars))
      haps <- cbind(haps, df[m,covars,drop=FALSE])

    haps$y <- df[ m, phenotype ]

  if(!length(covars)) {
    haps$covar <- "ALL"
    covars <- "ALL"
  }
  haps$cv <- apply(haps[,covars],1,paste,collapse="/")
  
    if(length(kk<-setdiff(use,haps$hap)))
      stop(length(kk)," haps not found")
    haps$hap[ !(haps$hap %in% use) ] <- "OTHER"
    lev <- c(lev,"OTHER")
##    haps <- subset(haps, hap %in% use)
    haps$hap <- factor(haps$hap, levels=lev)
  table(haps$cv,haps$hap)
  
}
##' Model haplotypes using multiple imputation
##'
##' .. content for \details{} ..
##' @title 
##' @param dir 
##' @param df 
##' @param family 
##' @param thr 
##' @param maxi 
##' @param haps.pattern 
##' @param covars 
##' @param phenotype 
##' @importFrom mice as.mira pool
##' @return 
##' @export
##' @author Chris Wallace
model.mi <- function(dir,df,family="binomial",thr=0.001,maxi=NULL,haps.pattern="haps.out2",
                     covars="q", phenotype="y",base=NULL) {

    ## fix haplotype base and set of haps to use
    haps <- hap.read(file.path(dir,haps.pattern))
    tt <- tapply(haps$Prob,haps$hap,sum)
    tt <- 100*tt/sum(tt)
    use <- names(tt)[tt/sum(tt)>thr]
    if(!is.null(base)) {
        if(!(base %in% use)) {
            message("haplotypes found:")
            cat(use,sep="\n")
            stop("requested base haplotype not found: ", base)
        }
    } else {
        base <- names(tt)[which.max(tt)]
    }
    lev <- names(sort(tt[use],decreasing=TRUE))
    message(length(lev)," haplotypes found meeting frequency threshold ",thr)
    lev <- c(lev,"OTHER")
    wh <- which(lev==base)
    if(wh!=1)
        lev <- c(lev[wh],lev[-wh])

    ## model formula
    f <- "y ~ hap"
    if(length(covars)) {
        covars <- intersect(covars,colnames(df))
        f <- paste(f, paste(covars, collapse="+"), sep=" + ")
    }
    
    ## read MI and fit models
    files <- list.files(dir,pattern=paste0(haps.pattern,"."),full=TRUE)
    if(!is.null(maxi) && length(files)>maxi)
        files <- files[1:maxi]
    message("reading ",length(files)," multiply imputed datasets from ",dir)
    fits <- lapply(seq_along(files), function(i) {
        cat(".")
        haps <- hap.read(files[[i]])
        haps <- haps[ haps$id %in% rownames(df), ]
        m <- match(haps$id, rownames(df))
        if(length(covars))
            haps <- cbind(haps, df[m,covars,drop=FALSE])
        
        haps$y <- df[ m, phenotype ]
        ## if(length(kk<-setdiff(use,haps$hap)))
        ##   stop(length(kk)," haps not found")
        haps$hap[ !(haps$hap %in% lev) ] <- "OTHER"
        ##    haps <- subset(haps, hap %in% use)
        haps$hap <- factor(haps$hap, levels=lev)
        glm(as.formula(f), data=haps, family=family)
    })
    
    ## make mira
    ## print(table(sapply(lapply(fits,coefficients),length)))
    message("averaging results over the MI fits")
    mi <- as.mira(fits)
    ss <- as.data.frame(summary(pool(mi)))
    null.value=0
    rownames(ss)[1] <- paste0("hap",base)
    ss$Fq <- unname(tt[sub("hap","",rownames(ss))])
    ss$est[1] <- 0
    ss$se[1] <- NA
    ss$Pr[1] <- NA
    ss$p <- ss$Pr
    ss$beta <- ss$est
    ss <- ss[,c("beta","se","p","Fq")]
    invisible(ss)                   
}
expand.haps <- function(m,snps,df) {
    h <- structure(vector("list",length(snps)), names=paste(snps,df[snps,"dbSNP"],sep="/"))
    msnps <- sub("hap","",rownames(m))
    for(i in seq_along(snps)) {
        h[[i]] <- substr(msnps,i,i)
    }
    cbind(as.data.frame(h),m)
}


getdips <- function(f) {
    haps <- read.table(f,header=TRUE,sep="\t",as.is=TRUE)
                                        #  haps <- haps[,-c(2,3)]
  haps$hap <- apply(haps[,-c(1:3,ncol(haps))],1,paste,collapse="")
  haps <- haps[,c(1:3,ncol(haps)-1,ncol(haps))]
  haps <- haps[order(haps$id,haps$ass,haps$chr,haps$hap),]
  dips <- merge(subset(haps,chr==1,select=c("id","ass","Probability","hap")),
                subset(haps,chr==2,select=c("id","ass","hap")),
                by=c("id","ass"),suffixes=c(".1",".2"))
  dips$dip <- paste(dips$hap.1,dips$hap.2,sep="/")
   return(dips)
}

diplotypes.mi <- function(dir,df,family="binomial",thr=0.001,maxi=NULL,haps.pattern="all-haps.out",
                     covars="q", phenotype="y") {

  ## fix haplotype base and set of haps to use
  dips <- getdips(paste0(dir,"/",haps.pattern))
  tt <- tapply(dips$Prob,dips$dip,sum)
  tt <- 100*tt/sum(tt)
  use <- names(tt)[tt/sum(tt)>thr]
  base <- names(tt)[which.max(tt)]
  lev <- names(sort(tt[use],decreasing=TRUE))

  ## model formula
  f <- "y ~ dip"
  if(length(covars)) {
    covars <- intersect(covars,colnames(df))
    f <- paste(f, paste(covars, collapse="+"), sep=" + ")
  }

  ## read MI and fit models
  files <- list.files(dir,pattern=paste0(haps.pattern,"."),full=TRUE)
  if(!is.null(maxi) && length(files)>maxi)
    files <- files[1:maxi]
  message("reading ",length(files)," multiply imputed datasets from ",dir)
  fits <- lapply(seq_along(files), function(i) {
    cat(".")
    dips <- getdips(files[[i]])
    dips <- dips[ dips$id %in% rownames(df), ]
    m <- match(dips$id, rownames(df))
    if(length(covars))
      dips <- cbind(dips, df[m,covars,drop=FALSE])

    dips$y <- df[ m, phenotype ]
    if(length(kk<-setdiff(use,dips$dip)))
      stop(length(kk)," dips not found")
    dips$dip[ !(dips$dip %in% use) ] <- "OTHER"
    lev <- c(lev,"OTHER")
##    dips <- subset(dips, dip %in% use)
    dips$dip <- factor(dips$dip, levels=lev)
    glm(as.formula(f), data=dips, family=family)
  })

  ## make mira
  print(table(sapply(lapply(fits,coefficients),length)))
  message("averaging results over the MI fits")
  mi <- as.mira(fits)
  ss <- as.data.frame(summary(pool(mi)))
  null.value=0
  rownames(ss)[1] <- paste0("dip",base)
  ss$Fq <- tt[sub("dip","",rownames(ss))]
  ss$est[1] <- 0
  ss$se[1] <- NA
  ss$Pr[1] <- NA
  ss$p <- ss$Pr
  ss$beta <- ss$est
  ss <- ss[,c("beta","se","p","Fq")]
  invisible(ss)                   
}

library(reshape)
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param results
##' @param hsnps
##' @param hfreq
##' @param thr
##' @param y.max
##' @param y.min
##' @return 
##' @export
##' @author Chris Wallace
plotter <- function(results,hsnps,hfreq,thr=0.01,y.max=NA,y.min=NA) {

  if(!is.na(thr)) {
    results <- subset(results,Fq>(thr*100))
    haps.found <- levels(results$hap)[ levels(results$hap) %in% results$hap ]
    results$hap <- factor(as.character(results$hap),levels=haps.found)
    hfreq <- hfreq[,sub("hap","",levels(results$hap)),drop=FALSE]
  }
  
  if(!is.na(y.max) & !is.na(y.min)) {
    results <- within(results, {
      ymin <- pmax(beta - 1.96*se, y.min)
      ymax <- pmin(beta + 1.96*se, y.max)
    })
  } else {
    results <- within(results, {
      ymin <- beta - 1.96*se
      ymax <- beta + 1.96*se
    })    
  }
  
  diplotypes <- "dip" %in% colnames(results)
  if(diplotypes) {
    results$hap <- results$dip
    hsnps <- c(hsnps,"/",hsnps)
  }

  hlab <- results[!duplicated(results[,"hap"]),]
  hlab <- within(hlab, {
    x <- as.numeric(hlab$hap) + 0.25
    hap <- sub("hap|dip","",hap)
  }) %>%
    subset(., hap!="OTHER" & Fq>thr*100, select=c("hap","x"))
  hlab <- cbind(hlab, strsplit(hlab$hap,"") %>% do.call("rbind",.))

  h <- melt(hlab,c("hap","x"))
  h$snp <- hsnps[as.numeric(h$variable)]
  major <- melt(hlab[1,-which(colnames(hlab)=="x")],c("hap"))
  major$major <- TRUE
  if(diplotypes) {
    dash <- which(hsnps=="/")
    major$value[(dash+1):length(hsnps)] <- major$value[1:(dash-1)]
  }
  h <- merge(h,major[,c("variable","value","major")],all=TRUE)
  h$major[is.na(h$major)] <- FALSE
  h$y <- as.numeric(h$variable)

  hplot <- ggplot(h,aes(x=x,y=y,label=value,colour=major,xmin=x-0.5,xmax=x+0.5,ymin=y-0.5,ymax=y+0.5,fill=major)) + 
    geom_rect() +
      geom_text(angle=90) +
  theme_minimal() +
    scale_y_continuous(breaks=seq_along(hsnps),labels=hsnps,limits=c(0.5,length(hsnps) + 0.5)) +
      scale_fill_manual("major allele",values=c("grey20","white")) +
        scale_colour_manual("major allele",values=c("white","grey20")) + ylab("SNP") +
         scale_x_continuous("Haplotype",breaks=as.numeric(hlab$hap)+0.25,
                            labels=sub("hap","",hlab$hap)) +
  theme(legend.position="bottom")

  colours <- scales::hue_pal()(length(levels(results$disease)))
  
  vlines <- unique(as.numeric(results$hap)) - 0.25
  oplot <- ggplot(results, aes(x=as.numeric(hap)+as.numeric(disease)*0.1,y=beta,ymin=ymin,ymax=ymax,col=disease)) + geom_pointrange() +
    geom_hline(yintercept=0) +
      geom_vline(xintercept=vlines,colour="grey") +
      scale_x_continuous("Haplotype",breaks=as.numeric(hlab$hap)+0.25,
                         labels=sub("hap","",hlab$hap)) +
                           scale_colour_manual(values=colours) +
                           theme_minimal() +
                           theme(axis.text.x=element_text(angle=90,family="monospace",size=12),
                                 legend.position="top") +
                                   ylab("Log OR")
  if(!is.na(y.min) & !is.na(y.max))
    oplot <- oplot + ylim(y.min,y.max)

  if(diplotypes) {
    dfreq <- matrix(NA,nrow(hfreq),nrow(hlab),dimnames=list(rownames(hfreq),hlab$hap))
    ss <- strsplit(colnames(dfreq),"/")
    n <- structure(1:ncol(hfreq),names=colnames(hfreq))
    homs <- sapply(ss,function(x) x[[1]]==x[[2]])
    for(i in seq_along(ss)) {
      ni <- n[ss[[i]]]
      if(ni[1]==ni[2]) {
        dfreq[,i] <- hfreq[,ni[1]]^2
      } else {
        dfreq[,i] <- 2 * hfreq[,ni[1]] * hfreq[,ni[2]]
      }
    }
    hfreq <- dfreq
  }
  
  h2 <- melt(hfreq,varnames=c("group","hap"))
  h2 <- subset(h2, hap %in% hlab$hap & !grepl("/NA",group))
  ss <- strsplit(as.character(h2$group),"/")
  h2$disease <- sapply(ss,"[[",2)
  h2$country <- sapply(ss,"[[",1)
  h2 <- subset(h2,disease %in% c(levels(results$disease)## ,"CONTROL"
                                 ) & country=="UK")
  h2$disease <- factor(h2$disease,levels=c(levels(results$disease),"CONTROL"))
  h2$hap <- factor(h2$hap,levels=hlab$hap)

  fplot <- ggplot(h2, aes(x=as.numeric(hap)+0.1*as.numeric(disease),y=100*value,ymin=0,ymax=100*value,col=disease)) + geom_pointrange() +
    geom_hline(yintercept=0) +
      geom_vline(xintercept=vlines,colour="grey") +
      scale_x_continuous("Haplotype",breaks=as.numeric(hlab$hap)+0.25,
                         labels=sub("hap","",hlab$hap)) +
                           scale_colour_manual(values=c(## "grey30",
                                                 colours)) +
                           theme_minimal() +
                           theme(axis.text.x=element_text(angle=90,family="monospace",size=12),
                                 legend.position="top") +
                             scale_y_continuous("Fq (%)")
  
  return(list(haps=hplot,odds=oplot,freq=fplot,hlab=hlab))
}
