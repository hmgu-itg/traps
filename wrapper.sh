#!/bin/bash

# simulates 1000 causal variants, 0.4 heritability, 100k samples total, 80% of them European, and half of positions having heter effects between anc groups
## TODO add the lognormal model
## TODO add the Wu model (Andrew's favourite probably)

nsnp=$1
her=$2
n=$3
propheter=$4
anc=$5
baseout=$6
pthr=$7

mkdir snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout && cd snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout && cp -r /treps/* .
baseout=test100k
./simulate -n $nsnp -r $her -s $n -p $propheter -a $anc -o $baseout


for pop in GWD LWK MSL YRI CLM MXL PEL PUR CDX CHB JPT KHV FIN GBR IBS TSI BEB GIH PJL STU; do 
  awk '{if(!($1~/'$pop'/)){$NF=-9}}1' 'OFS=\t' $baseout.tfam > $pop.tfam; 
  plink --tped $baseout.tped --tfam $pop.tfam --allow-no-sex --assoc --out $pop
  awk '{if(!($1~/'$pop'/)){$NF=-9}}1' 'OFS=\t' $baseout.target.tfam > $pop.target.tfam
  plink --tped $baseout.target.tped --tfam $pop.target.tfam --allow-no-sex --assoc --out $pop.target
done


./makemeta.R $(pwd) .target
for f in `ls metasoft.*.txt`; do cp $f $f.bak; awk 'NR>1' $f | sponge $f; done

for p in AFR AMR EAS EUR SAS TE; do java -jar /Metasoft.jar -input metasoft.$p.txt -pvalue_table /HanEskinPvalueTable.txt -output metasoft.$p.out; done
./mvmeta.R

plink --tfile $baseout.target --recode A-transpose --out $baseout.target.add

#sadness, the metasoft.out files have garbage columns at the end (and spelling mistakes in col names :((( )
for f in `ls metasoft.*.out`; do cut -f1-18 $f| sponge $f; done
./processmeta.R -p $pthr -a $anc -s $n -i $baseout -o snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout.out

tar -cvjf snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout.rundata.tar.bz2  *fam *.log *.qassoc *.out metasoft.* mvmeta* $baseout.*
cp snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout.*txt snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout.rundata.tar.bz2 /transfer

sleep 5

