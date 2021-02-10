#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(gtx))
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
  make_option(c("-i", "--infile"), type="character", default=NULL,
              help="Base file input, .target.add.traw and .target.tfam expected.", metavar="character"),
    make_option(c("-o", "--out"), type="character",
              help="output file prefix.", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt[["num-samples"]]) |
    is.null(opt[["anc-af-distr"]]) |
    is.null(opt[["p-thresh"]]) |
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
cat(paste("Reading target additive dataset\n"))
flush.console()
target=fread(paste0(indir, "/",infile,".target.add.traw"))
ttfam=fread(paste0(indir, "/",infile,".target.tfam"))

toselect=c(2,7:ncol(target))
target=target[,..toselect]
setnames(target, "SNP", "rs")
setnames(target, colnames(target), sub("_.*", "", colnames(target)))


cat(paste("Reading ancestry-specific results\n"))
flush.console()
mvmeta=fread(paste0(indir,"/mvmeta.out"))
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
popdat=list()
for(p in unique(nmtbl$V1)){
    if(file.exists(paste0(indir,"/",p,".target.qassoc"))){
            popdat[[p]]=fread(paste0(indir,"/",p,".target.qassoc"),select=c("SNP","BETA" ,"SE"))
    }
}

resul=NULL

## Ancestry specific stuff first

for(p in unique(nmtbl$V1)){
    ### indirect
    angroup=nmtbl$V2[nmtbl$V1==p]
    tosimulate=anc.distr[angroup]*n.samples
    POPPERSUP=4
    perpop=tosimulate/POPPERSUP
    popcol=colnames(target)[grep(p, colnames(target))]
    if(!length(popcol)){next}
    cat(paste("\t\t Ancestry-specific, indirect method", p,"\n"))
    flush.console()
    totest=merge(tst[,c("SNP", p),with=F], popdat[[p]], by="SNP")
    totest=merge(totest, mvmeta[,c("SNP", "V2"), with=F])
    setnames(totest, c("SNP", p, "BETA", "SE", "V2"), c("SNP","mvbeta", "beta", "se","P"))
    toconvert=c("mvbeta", "beta", "se","P")
    totest[, (toconvert) := lapply(.SD, as.numeric), .SDcols = toconvert]
    totest=totest[P<P_THRESHOLD,]
    gtxs=grs.summary(totest$mvbeta, totest$beta, totest$se, perpop)
    add=data.table(PRS=p, target=p, method="indirect", cor=gtxs$ahat[1], Standard.Error=gtxs$aSE[1], P=gtxs$pval[1], varex=gtxs$R2m[1])
    resul=rbindlist(list(resul, add), use.names=T)
    ## direct
    cat(paste("\t\t Ancestry-specific, direct method", p,"\n"))
    flush.console()
    m=merge(tst[,c("SNP", p),with=F], target, by.x="SNP", by.y="rs")
    m=merge(m, mvmeta[,c("SNP", "V2"), with=F],by="SNP")
    setnames(m, "V2", "P")
    m=m[P<P_THRESHOLD,]
    toconvert=p
    m[, (toconvert) := lapply(.SD, as.numeric), .SDcols = toconvert]
    prs=t(t(lapply(m[,..popcol], function(x){sum(m[[p]]*x)})))
    prs=data.table(id=rownames(prs), prs=unlist(prs))
    prs=merge(prs, ttfam, by.x="id", by.y="V1")
    add=cor.test.plus(prs$prs, prs$V6)
    add$PRS=p
    add$target=p
    add$varex=add$cor*add$cor
    add$method="direct"
    resul=rbindlist(list(resul, add), use.names=T)
}

for (angroup in c("AFR", "AMR", "EAS", "EUR", "SAS", "TE")){
    cat(paste("Processing score from ", angroup,"\n"))
    flush.console()
    ancmeta=fread(paste0(indir, "/metasoft.",angroup,".out"),
    select=c("RSID","PVALUE_FE" ,"BETA_FE"))
    ancmeta=ancmeta[PVALUE_FE<P_THRESHOLD,c("RSID", "BETA_FE")]
    cat(paste("Number of SNPs in PRS:",nrow(ancmeta), "\n"))
    for(p in unique(nmtbl$V1)){
        ### Exact method corr with phenotype
        cat(paste("\t Using target population", p,"\n"))
        cat(paste("\t\t Direct method", p,"\n"))
        flush.console()
        m=merge(ancmeta, target, by.x="RSID", by.y="rs")
        popcol=colnames(target)[grep(p, colnames(target))]
        if(!length(popcol)){next}
        prs=t(m[,lapply(.SD, function(x){sum(BETA_FE*x)}), .SDcols=popcol])
        prs=data.table(id=rownames(prs), prs=prs)
        prs=merge(prs, ttfam, by.x="id", by.y="V1")
        add=cor.test.plus(prs$prs.V1, prs$V6)
        add$PRS=angroup
        add$target=p
        add$method="direct"
        add$varex=add$cor*add$cor
        if(is.null(resul)){resul=add}else{resul=rbind(resul,add)}

        ## approx method
        cat(paste("\t\t Indirect method", p,"\n"))
        flush.console()
        if(angroup=="TE"){perpop=n.samples}else{
          tosimulate=anc.distr[angroup]*n.samples
          POPPERSUP=4
          perpop=tosimulate/POPPERSUP
        }

        ancmeta.pop=merge(ancmeta, popdat[[p]], all.x=T, by.x="RSID", by.y="SNP", suffixes = c("", paste0(".", p)))
        gtxs=grs.summary(ancmeta.pop$BETA_FE, ancmeta.pop$BETA, ancmeta.pop$SE, perpop)
        add=data.table(PRS=angroup, target=p, method="indirect", cor=gtxs$ahat[1], Standard.Error=gtxs$aSE[1], P=gtxs$pval[1], varex=gtxs$R2m[1])
        resul=rbindlist(list(resul, add), use.names=T)

    }
}
fwrite(resul, paste(opt[["out"]], "prs.fit.txt", sep="."), sep="\t", quote=F)
