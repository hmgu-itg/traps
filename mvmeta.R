#!/usr/bin/env Rscript
library(data.table)
indir="."
for(suffix in c("pn", "ln", "Wu")){
    mvbeta=fread(paste0(indir, "/mvmeta.",suffix,".beta"))
    mvse=fread(paste0(indir, "/mvmeta.",suffix,".se"))
    pca= fread('pca.txt',header=F)

    # two lines below are because pops in PCA.txt are not named.
    pops=c("GWD", "LWK", "MSL", "YRI", "CLM", "MXL", "PEL", "PUR", "CDX", "CHB", "JPT", "KHV", "FIN", "GBR", "IBS", "TSI", "BEB", "GIH", "PJL", "STU")
    pca$pop=pops


    m=merge(mvbeta, mvse, by="SNP")
    out=t(apply(m, 1, function(x){
        betavec=paste("BETA", pops, sep=".")
        betavec=as.numeric(x[betavec])
        sevec=paste("SE", pops, sep=".")
        sevec=as.numeric(x[sevec])
        sevec=sevec*sevec
        sevec=1/sevec
        pca1=pca[match(pops, pop), V1]
        pca2=pca[match(pops, pop), V2]
        pca3=pca[match(pops, pop), V3]
        model0 = glm(betavec~0,family="gaussian",weights=sevec)
        model1 = glm(betavec~1,family="gaussian",weights=sevec)
        model2 = glm(betavec~pca1+pca2+pca3,family="gaussian",weights=sevec)
        c(x["SNP"], pchisq(anova(model0,model2)[2,4],df=4,lower.tail=F),
              pchisq(anova(model1,model2)[2,4],df=3,lower.tail=F),
              pchisq(summary(model2)$deviance,df=(length(x)-1)/2-4,lower.tail=F),
              summary(model2)$coefficients[1,1],
              summary(model2)$coefficients[1,2],
              summary(model2)$coefficients[2,1],
              summary(model2)$coefficients[2,2],
              summary(model2)$coefficients[2,4],
              summary(model2)$coefficients[3,1],
              summary(model2)$coefficients[3,2],
              summary(model2)$coefficients[3,4],
              summary(model2)$coefficients[4,1],
              summary(model2)$coefficients[4,2],
              summary(model2)$coefficients[4,4],
              anova(model0,model2)[2,4],
              (length(x)-1)/2)
    }))
    out=as.data.table(out)
    fwrite(out,paste0("mvmeta.",suffix,".out"),quote=F, sep="\t")
}
