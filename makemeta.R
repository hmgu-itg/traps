#!/usr/bin/env Rscript
library(data.table)
args=commandArgs(T)

midfix=args[2]
dir=args[1]
models=if(args[3]=="all") c("pn", "ln", "Wu") else "ln"
for(suffix in models){
    ftoread=list.files(dir, paste0(midfix, ".", suffix,".qassoc$"))
    pops=sub(paste0(midfix, ".", suffix, ".qassoc"), "", ftoread, fixed=T)

    nmtbl=fread(text="CHB	EAS
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

    ol=NULL
    for(angroup in c(unique(nmtbl$V2), "TE")){
        if(angroup=="TE"){pops=nmtbl$V1}else{pops=unlist(nmtbl$V1[nmtbl$V2==angroup])}
        group.mtx=NULL
        for(p in pops){
            fn=paste0(dir,"/", p,".", suffix, ".qassoc")
            if(!file.exists(fn)){print("did not find file");print(fn);next}
            pdat=fread(fn,
                       select = c("SNP", "BETA", "SE"))
            setnames(pdat, c("SNP", paste(c("BETA", "SE"), p, sep=".")))
            if(is.null(group.mtx)){
                group.mtx=pdat
            }else{
                group.mtx=merge(group.mtx, pdat, by="SNP")
            }
        }
        print(head(group.mtx))
        fwrite(group.mtx,
               paste0(dir,"/metasoft.", angroup, ".", suffix,".txt")
              ,sep="\t",
              quote=F,
              col.names=T)
              if(angroup != "TE"){
                if(is.null(ol)){ol=group.mtx}else{ol=merge(ol, group.mtx, by="SNP")}
              }
    }

    betacol=c(1, seq(1, ncol(ol)-1, by=2)+1)
    secol=c(1, seq(2, ncol(ol)-1, by=2)+1)

    pops=c("GWD", "LWK", "MSL", "YRI", "CLM", "MXL", "PEL", "PUR", "CDX", "CHB", "JPT", "KHV", "FIN", "GBR", "IBS", "TSI", "BEB", "GIH", "PJL", "STU")
    ordercol=paste("BETA", pops, sep=".")
    betadf=ol[,..betacol]
    setcolorder(betadf, c("SNP", ordercol))
    ordercol=paste("SE", pops, sep=".")
    sedf=ol[,..secol]
    setcolorder(sedf, c("SNP", ordercol))


    fwrite(betadf,
          paste0(dir,"/mvmeta.",suffix,".beta"),
          sep="\t",
          quote=F,
          col.names=T)
    fwrite(sedf,
          paste0(dir, "/mvmeta.",suffix,".se"),
          sep="\t",
          quote=F,
          col.names=T)
}
