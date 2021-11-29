#!/bin/bash

## TODO add the Wu model

nsnp=$1
her=$2
n=$3
propheter=$4
anc=$5
baseout=$6
pthr=$7

if [[ -z "$8" ]]; then
  SAVEFILE="yes"
else
  SAVEFILE=$8
fi

echo parameters $nsnp $her $n $propheter $anc $baseout $pthr

if [[ -s "/transfer/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout.out.prs.fit.txt" ]]; then
  echo "File /transfer/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout.out.prs.fit.txt exists, nothing to do."
  exit 0
fi

if [ -z "$SCRATCH" ]; then
echo Current directory $(pwd), set SCRATCH env var to change.
mkdir -p tempFiles/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout && cd tempFiles/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout && cp -r /treps/* .

else
 mkdir -p $SCRATCH/ge64cig2/TransAncPRS/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout && cd $SCRATCH/ge64cig2/TransAncPRS/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout && cp -r /treps/* .
fi

if [[ "$SAVEFILE" == "yes" ]]; then
 ./simulate -n $nsnp -r $her -s $n -p $propheter -a $anc -o $baseout
else
  ./simulate -n $nsnp -r $her -s $n -p $propheter -a $anc -o $baseout -x
fi

for pop in GWD LWK MSL YRI CLM MXL PEL PUR CDX CHB JPT KHV FIN GBR IBS TSI BEB GIH PJL STU; do
  awk '{if(!($1~/'$pop'/)){$NF=-9}}1' 'OFS=\t' $baseout.pn.tfam > $pop.pn.tfam;
  plink --tped $baseout.tped --tfam $pop.pn.tfam --allow-no-sex --assoc --out $pop.pn
  awk '{if(!($1~/'$pop'/)){$NF=-9}}1' 'OFS=\t' $baseout.target.pn.tfam > $pop.target.pn.tfam
  plink --tped $baseout.target.tped --tfam $pop.target.pn.tfam --allow-no-sex --assoc --out $pop.target.pn

  awk '{if(!($1~/'$pop'/)){$NF=-9}}1' 'OFS=\t' $baseout.ln.tfam > $pop.ln.tfam;
  plink --tped $baseout.tped --tfam $pop.ln.tfam --allow-no-sex --assoc --out $pop.ln
  awk '{if(!($1~/'$pop'/)){$NF=-9}}1' 'OFS=\t' $baseout.target.ln.tfam > $pop.target.ln.tfam
  plink --tped $baseout.target.tped --tfam $pop.target.ln.tfam --allow-no-sex --assoc --out $pop.target.ln

done


./makemeta.R $(pwd) .target
for f in `ls metasoft.*.txt`; do cp $f $f.bak; awk 'NR>1' $f | sponge $f; done

for p in AFR AMR EAS EUR SAS TE; do java -jar /Metasoft.jar -input metasoft.$p.pn.txt -pvalue_table /HanEskinPvalueTable.txt -output metasoft.$p.pn.out; done
for p in AFR AMR EAS EUR SAS TE; do java -jar /Metasoft.jar -input metasoft.$p.ln.txt -pvalue_table /HanEskinPvalueTable.txt -output metasoft.$p.ln.out; done

./mvmeta.R

plink --tfam $baseout.target.pn.tfam --tped $baseout.target.tped --recode A-transpose --out $baseout.target.add.pn
plink --tfam $baseout.target.ln.tfam --tped $baseout.target.tped --recode A-transpose --out $baseout.target.add.ln

#the metasoft.out files have garbage columns at the end (and spelling mistakes in col names)
for f in `ls metasoft.*.out`; do cut -f1-18 $f| sponge $f; done
./processmeta.R -p $pthr -a $anc -s $n -i $baseout -o snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout.out

if [[ "$SAVEFILE" == "yes" ]]; then
tar -cvjf snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout.rundata.tar.bz2  *fam *.log *.qassoc *.out metasoft.* mvmeta* $baseout.*
cp snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout.rundata.tar.bz2 /transfer
fi

gzip snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout.*txt
cp snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout.*txt.gz /transfer

if [ -z "$SCRATCH" ]; then
 cd ../..
 rm -r tempFiles/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout
else
 rm -r $SCRATCH/ge64cig2/TransAncPRS/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout
fi

sleep 5
