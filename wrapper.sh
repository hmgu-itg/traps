#!/bin/bash

## TODO add the Wu model

nsnp=$1
her=$2
n=$3
propheter=$4
anc=$5
baseout=$6
pthr=$7
GW=$8

if [[ -z "$9" ]]; then
  SAVEFILE="yes"
else
  SAVEFILE=$9
fi

echo parameters $nsnp $her $n $propheter $anc $baseout $pthr $GW $SAVEFILE

if [[ -s "/transfer/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout.out.prs.fit.txt" ]]; then
  echo "File /transfer/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout.out.prs.fit.txt exists, nothing to do."
  exit 0
fi

echo Current directory $(pwd), set SCRATCH env var to change.
mkdir -p results/all_iter
mkdir -p tempFiles/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.$baseout && cd tempFiles/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.$baseout && cp -r /treps/* .


if [[ "$SAVEFILE" == "yes" ]]; then
 ./simulate -n $nsnp -r $her -s $n -p $propheter -a $anc -w $GW -o $baseout
else
 ./simulate -n $nsnp -r $her -s $n -p $propheter -a $anc -w $GW -o $baseout -x
fi

##MERGE the files from the different populations
echo $baseout.EAS > simus_tomerge ; echo $baseout.EAS.target > simus_tomerge_target
for superpop in EUR AFR AMR; do
  echo $baseout.$superpop >> simus_tomerge
  echo $baseout.$superpop.target >> simus_tomerge_target
done

plink --bfile $baseout.SAS --merge-list simus_tomerge --merge-equal-pos --make-bed --out $baseout
plink --bfile $baseout.SAS.target --merge-list simus_tomerge_target --merge-equal-pos --make-bed --out $baseout.target


##PLINK association tests within each population
for sp in AFR AMR EAS EUR SAS; do
  pops=$(awk -v sp=$sp '{if($2 ~ sp) {print $1} }' ./1000G_populations.txt)
  for pop in $pops; do
    awk '{if(!($1~/'$pop'/)){$NF=-9}}1' 'OFS=\t' $baseout.$sp.pn.fam > $pop.pn.fam;
    plink --bed $baseout.$sp.bed --bim $baseout.$sp.bim --fam $pop.pn.fam --allow-no-sex --assoc --out $pop.pn
    awk '{if(!($1~/'$pop'/)){$NF=-9}}1' 'OFS=\t' $baseout.$sp.target.pn.fam > $pop.target.pn.fam
    plink --bed $baseout.$sp.target.bed --bim $baseout.$sp.target.bim --fam $pop.target.pn.fam --allow-no-sex --assoc --out $pop.target.pn

    awk '{if(!($1~/'$pop'/)){$NF=-9}}1' 'OFS=\t' $baseout.$sp.ln.fam > $pop.ln.fam;
    plink --bed $baseout.$sp.bed --bim $baseout.$sp.bim --fam $pop.ln.fam --allow-no-sex --assoc --out $pop.ln
    awk '{if(!($1~/'$pop'/)){$NF=-9}}1' 'OFS=\t' $baseout.$sp.target.ln.fam > $pop.target.ln.fam
    plink --bed $baseout.$sp.target.bed --bim $baseout.$sp.target.bim --fam $pop.target.ln.fam --allow-no-sex --assoc --out $pop.target.ln

    awk '{if(!($1~/'$pop'/)){$NF=-9}}1' 'OFS=\t' $baseout.$sp.Wu.fam > $pop.Wu.fam;
    plink --bed $baseout.$sp.bed --bim $baseout.$sp.bim --fam $pop.Wu.fam --allow-no-sex --assoc --out $pop.Wu
    awk '{if(!($1~/'$pop'/)){$NF=-9}}1' 'OFS=\t' $baseout.$sp.target.Wu.fam > $pop.target.Wu.fam
    plink --bed $baseout.$sp.target.bed --bim $baseout.$sp.target.bim --fam $pop.target.Wu.fam --allow-no-sex --assoc --out $pop.target.Wu
  done
done


./makemeta.R $(pwd) .target
for f in `ls metasoft.*.txt`; do cp $f $f.bak; awk 'NR>1' $f | sponge $f; done

for p in AFR AMR EAS EUR SAS TE; do java -jar /Metasoft.jar -input metasoft.$p.pn.txt -pvalue_table /HanEskinPvalueTable.txt -output metasoft.$p.pn.out; done
for p in AFR AMR EAS EUR SAS TE; do java -jar /Metasoft.jar -input metasoft.$p.ln.txt -pvalue_table /HanEskinPvalueTable.txt -output metasoft.$p.ln.out; done
for p in AFR AMR EAS EUR SAS TE; do java -jar /Metasoft.jar -input metasoft.$p.Wu.txt -pvalue_table /HanEskinPvalueTable.txt -output metasoft.$p.Wu.out; done

./mvmeta.R

#combine the files from all the populations
for m in pn ln Wu; do
  cat $baseout.AFR.target.$m.fam $baseout.AMR.target.$m.fam $baseout.EAS.target.$m.fam $baseout.EUR.target.$m.fam $baseout.SAS.target.$m.fam > $baseout.target.$m.fam
done 

#the metasoft.out files have garbage columns at the end (and spelling mistakes in col names)
for f in `ls metasoft.*.out`; do cut -f1-18 $f| sponge $f; done
./processmeta.R -p $pthr -a $anc -s $n -i $baseout -o snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.$baseout.out

if [[ "$SAVEFILE" == "yes" ]]; then
tar -cvjf snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.$baseout.rundata.tar.bz2  *fam *.log *.qassoc *.out metasoft.* mvmeta* $baseout.*
cp snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.$baseout.rundata.tar.bz2 ../../results/
fi

gzip snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.$baseout.*txt
cp snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.$baseout.*txt.gz ../../results/

cd ../..
rm -r tempFiles/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.$baseout

sleep 5
