#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(gtx))
suppressPackageStartupMessages(library(gaston))
options(warn=1)

TOTNUMPOP=5
POPPERSUP=4
P_THRESHOLD=5e-4

option_list = list(
  make_option(c("-s", "--num-samples"), type="integer", default=NULL,
              help="Total number of samples to simulate.", metavar="integer"),
  make_option(c("-a", "--anc-af-distr"), type="character", default=NULL,
              help="Comma-separated list of proportions to simulate from each ancestry group. Should be a vector of 5 probabilities and sum to 1. Alternatively, if this is only 1 probability, it is assumed to mean the proportion of EUR.", metavar="character"),
  make_option(c("-p", "--p-thresh"), type="numeric", default=NULL,
              help="P-value threshold at which to select SNPs for inclusion in the PRS.", metavar="numeric"),
  make_option(c("-m", "--model"), type="character", default=NULL,
              help="Run only the log-normal model ('ln') or the log-normal, point-normal and frequency-dependent models ('all')", metavar="character"),
  make_option(c("-i", "--infile"), type="character", default=NULL,
              help="Base file input, .target.add.traw and .target.fam expected.", metavar="character"),
  make_option(c("-o", "--out"), type="character",
              help="output file prefix.", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt[["num-samples"]]) |
    is.null(opt[["anc-af-distr"]]) |
    is.null(opt[["p-thresh"]]) |
    is.null(opt[["model"]]) |
    is.null(opt[["infile"]]) |
    is.null(opt[["out"]])
  ) {
  print_help(opt_parser)
  stop("All arguments are compulsory.\n", call.=FALSE)
}


ancdist=as.numeric(strsplit(opt[["anc-af-distr"]], ",", fixed=T)[[1]])
n.samples=opt[["num-samples"]]

P_THRESHOLD=1e-4
if(opt[["p-thresh"]]<=0 | opt[["p-thresh"]]>1){
  stop("P-value threshold must be in ]0,1].")
}else{
  P_THRESHOLD=opt[["p-thresh"]]
}
infile=opt[["infile"]]
indir="."
cor.test.plus <- function(x, y) {
    x=cor.test(x,y)
  data.table(cor=x$estimate,
       Standard.Error = unname(sqrt((1 - x$estimate^2)/x$parameter)), P=x$p.value)
}


nmtbl=read.table(text="CHB	EAS
JPT	EAS
CHS	EAS
CDX	EAS
KHV	EAS
CEU	EUR
TSI	EUR
FIN	EUR
GBR	EUR
IBS	EUR
YRI	AFR
LWK	AFR
GWD	AFR
MSL	AFR
ESN	AFR
ASW	AFR
ACB	AFR
MXL	AMR
PUR	AMR
CLM	AMR
PEL	AMR
GIH	SAS
PJL	SAS
BEB	SAS
STU	SAS
ITU	SAS", sep="\t", header=F, stringsAsFactors=F)

check_anc_distr=function(anc.distr){
    if(length(anc.distr)==1){
        if(anc.distr>1 | anc.distr<0){stop(paste("Proportion of Europeans (",anc.distr,") is out of bounds."))}
        eur=anc.distr
        anc.distr=rep((1-anc.distr)/4, 4)
        anc.distr=c(anc.distr[1:3], eur, anc.distr[4])
    }
    else{
        if(sum(anc.distr)!=1 | any(anc.distr>1) | any(anc.distr<0) |length(anc.distr)!=5){
            stop(paste("Ancestry distribution should be a vector of 5 probabilities summing to 1.\nYou provided",
                      paste(anc.distr, collapse=" + "), "=", sum(anc.distr)))
        }
    }
    return(anc.distr)
}

anc.distr=check_anc_distr(as.numeric(strsplit(opt[["anc-af-distr"]], ",", fixed=T)[[1]]))
names(anc.distr)=c("AFR", "AMR", "EAS", "EUR", "SAS")
print(anc.distr)

### MAIN
resul=NULL
models=if(opt[["model"]]=="all") c("pn", "ln", "Wu") else "ln"
for(suffix in models){
  if(suffix=="pn"){cat(paste("Pointnormal model\n"))}
  if(suffix=="ln"){cat(paste("Lognormal model\n"))}
  if(suffix=="Wu"){cat(paste("Frequency-dependent model\n"))}
  cat(paste("Reading target additive dataset\n"))
  flush.console()
  target.bm=read.bed.matrix(paste0(indir, "/", infile, ".target"))
  ttfam=fread(paste0(indir, "/",infile,".target.",suffix,".fam"))
  #tmp.names <- strsplit(sub(target.bm@ped$id, pattern="\\.V", replacement = "\\."), split = "\\.")
  #target.bm@ped$id <- paste0("SAMPLE.", unlist(lapply(tmp.names, function(z) paste0(z[2], ".", z[1]))))
  
  mvmeta=fread(paste0(indir,"/mvmeta.",suffix,".out"))
  pca=fread(paste0(indir,"/pca.txt"))

  cat(paste("Building ancestry-specific scores\n"))
  flush.console()
  tst=as.data.table(t(apply(mvmeta, 1, function(x){
      snp=x["SNP"]
      x=as.numeric(x[-1])
      c(snp, t(x[4]+pca$V1*x[6]+pca$V2*x[9]+pca$V3*x[12]))})))
  pops=c("GWD", "LWK", "MSL", "YRI", "CLM", "MXL", "PEL", "PUR", "CDX", "CHB", "JPT", "KHV", "FIN", "GBR", "IBS", "TSI", "BEB", "GIH", "PJL", "STU")
  setnames(tst, c("SNP", pops))


  cat(paste("Reading target association sumstats\n"))
  flush.console()


  ## Ancestry specific stuff first
  for(p in unique(nmtbl$V1)){
      if(file.exists(paste0(indir,"/",p,".target.",suffix,".qassoc"))){
              popdat=fread(paste0(indir,"/",p,".target.",suffix,".qassoc"),select=c("SNP","BETA" ,"SE"))
      }
 
      ### indirect
      angroup=nmtbl$V2[nmtbl$V1==p]
      tosimulate=anc.distr[angroup]*n.samples
      POPPERSUP=4
      popcol=ttfam$V1[grep(p, ttfam$V1)]
      perpop=tosimulate/POPPERSUP
      if(!length(popcol)){next}
      cat(paste("\t\t Ancestry-specific, indirect method", p,"\n"))
      flush.console()
      totest=merge(tst[,c("SNP", p),with=F], popdat, by="SNP")
      totest=merge(totest, mvmeta[,c("SNP", "V2"), with=F])
      setnames(totest, c("SNP", p, "BETA", "SE", "V2"), c("SNP","mvbeta", "beta", "se","P"))
      toconvert=c("mvbeta", "beta", "se","P")
      totest[, (toconvert) := lapply(.SD, as.numeric), .SDcols = toconvert]
      totest=totest[P<P_THRESHOLD,]
           
      add=NULL
      if(nrow(totest)>0){
        gtxs=grs.summary(totest$mvbeta, totest$beta, totest$se, perpop)
        add=data.table(model=suffix, PRS=p, target=p, nsnp=nrow(totest), method="indirect", cor=sqrt(gtxs$R2rs[1]), Standard.Error=gtxs$aSE[1], P=gtxs$pval[1], varex=gtxs$R2rs[1])
      }else{
        add=data.table(model=suffix, PRS=p, target=p, nsnp=0, method="indirect", cor=NA, Standard.Error=NA, P=NA, varex=NA)
      }
      resul=rbindlist(list(resul, add), use.names=T)

      ## direct
      cat(paste("\t\t Ancestry-specific, direct method", p,"\n"))
      flush.console()
      m=merge(tst[,c("SNP", p),with=F], mvmeta[,c("SNP", "V2"), with=F],by="SNP")
      setnames(m, "V2", "P")
      m=m[P<P_THRESHOLD,]
      #Select those snps in bed matrix
      target.bm.select <- select.snps(target.bm, id %in% m$SNP)
      target=cbind(target.bm.select@snps$id, as.data.table(abs(2-t(as.matrix(target.bm.select))))) #Have to change alleles because alleles are flipped
      m=merge(m, target, by.x="SNP", by.y="V1")

      add=NULL
      if(nrow(m)>0){
        toconvert=p
        m[, (toconvert) := lapply(.SD, as.numeric), .SDcols = toconvert]
        prs=t(t(lapply(m[,..popcol], function(x){sum(m[[p]]*x)})))
        prs=data.table(id=rownames(prs), prs=unlist(prs))
        prs=merge(prs, ttfam, by.x="id", by.y="V1")
        add=cor.test.plus(prs$prs, prs$V6)
        add$model=suffix
        add$PRS=p
        add$target=p
        add$nsnp=nrow(m)
        add$varex=add$cor*add$cor
        add$method="direct"
      }else{
        add=data.table(model=suffix, PRS=p, target=p, nsnp=0, method="direct", cor=NA, Standard.Error=NA, P=NA, varex=NA)
      }
      resul=rbindlist(list(resul, add), use.names=T)
  }


  for (angroup in c("AFR", "AMR", "EAS", "EUR", "SAS", "TE")){
      cat(paste("Processing score from ", angroup,"\n"))
      flush.console()
      ancmeta=fread(paste0(indir, "/metasoft.",angroup,".",suffix,".out"),
      select=c("RSID","PVALUE_FE" ,"BETA_FE", "PVALUE_RE2", "BETA_RE"))
      #FE and RE models
      ancmeta_FE=ancmeta[PVALUE_FE<P_THRESHOLD,c("RSID", "BETA_FE")]
      ancmeta_RE=ancmeta[PVALUE_RE2<P_THRESHOLD,c("RSID", "BETA_RE")]
      colnames(ancmeta_FE) <- colnames(ancmeta_RE) <- c("RSID", "BETA_META")
      ancmeta=list(FE = ancmeta_FE, RE = ancmeta_RE)
      rm(ancmeta_FE, ancmeta_RE)
      for(meta_model in c("FE", "RE")){
        cat(paste("Number of SNPs in PRS :",nrow(ancmeta[[meta_model]]), "\n"))

        if(nrow(ancmeta[[meta_model]])==0){
          add=data.table(target=unique(nmtbl$V1))
          add[,PRS:=paste(angroup, meta_model, sep = "_")]
          add[,model:=suffix]
          add[,nsnp:=0]
          add=rbind(add,add)
          add[,method:=rep(c("indirect", "direct"), each=nrow(add)/2)]
          add[,cor:=NA]
          add[,Standard.Error:=NA]
          add[,P:=NA]
          add[,varex:=NA]
          resul=rbindlist(list(resul, add), use.names=T)
        }
        else{
          for(p in unique(nmtbl$V1)){
            if(file.exists(paste0(indir,"/",p,".target.",suffix,".qassoc"))){
                popdat=fread(paste0(indir,"/",p,".target.",suffix,".qassoc"),select=c("SNP","BETA" ,"SE"))
            }
            ### Exact method corr with phenotype
            cat(paste("\t Using target population", p, meta_model, "model\n"))
            cat(paste("\t\t Direct method", p, meta_model, "model\n"))
            flush.console()
             #Select those snps in bed matrix
            target.bm.select <- select.snps(target.bm, id %in% ancmeta[[meta_model]]$RSID)
            target=cbind(target.bm.select@snps$id, as.data.table(abs(2-t(as.matrix(target.bm.select)))))
            
            m=merge(ancmeta[[meta_model]], target, by.x="RSID", by.y="V1")
            popcol=colnames(target)[grep(p, colnames(target))]
            if(!length(popcol)){next}
            prs=t(m[,lapply(.SD, function(x){sum(BETA_META*x)}), .SDcols=popcol]) 
            prs=data.table(id=rownames(prs), prs=prs)
            prs=merge(prs, ttfam, by.x="id", by.y="V1")
            add=cor.test.plus(prs$prs.V1, prs$V6)
            add$model=suffix
            add$PRS=paste(angroup, meta_model, sep = "_") 
            add$target=p
            add$nsnp=nrow(ancmeta_subset)
            add$method="direct"
            add$varex=add$cor*add$cor
            if(is.null(resul)){resul=add}else{resul=rbindlist(list(resul,add), use.names=T)}

            ## approx method
            cat(paste("\t\t Indirect method", p, meta_model, "model\n"))
            flush.console()
            if(angroup=="TE"){perpop=n.samples}else{
              tosimulate=anc.distr[angroup]*n.samples
              POPPERSUP=4
              perpop=tosimulate/POPPERSUP
            }

            ancmeta.pop=merge(ancmeta[[meta_model]], popdat, all.x=T, by.x="RSID", by.y="SNP", suffixes = c("", paste0(".", p)))
            gtxs=grs.summary(unlist(ancmeta.pop$BETA_META), ancmeta.pop$BETA, ancmeta.pop$SE, perpop)
            add=data.table(model=suffix, PRS=paste(angroup, meta_model, sep = "_"), target=p, nsnp=nrow(ancmeta.pop), method="indirect", cor=sqrt(gtxs$R2rs[1]), Standard.Error=gtxs$aSE[1], P=gtxs$pval[1], varex=gtxs$R2rs[1])
            resul=rbindlist(list(resul, add), use.names=T)
        }

      }
    }
  }
}
fwrite(resul, paste(opt[["out"]], "prs.fit.txt", sep="."), sep="\t", quote=F, na="NA")
