#!/bin/bash
model='ln'
GW='FALSE'
SAVEFILE='no'
het_effectsize='FALSE'

while getopts "S:H:n:h:a:b:p:gmes" opt; do
  case $opt in
    S) nsnp="$OPTARG" ;;
    H) her="$OPTARG" ;;
    n) n="$OPTARG" ;;
    h) propheter="$OPTARG" ;;
    a) anc="$OPTARG" ;;
    b) baseout="$OPTARG" ;;
    p) pthr="$OPTARG" ;;
    g) GW='TRUE' ;;
    m) model='all' ;;
    e) het_effectsize='TRUE' ;;
    s) SAVEFILE='yes' ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done


if [ -z "$nsnp" ]; then
  echo 'Missing -S (number of SNPs)' >&2; exit 1
fi
if [ -z "$her" ]; then
  echo 'Missing -H (heritability)' >&2; exit 1
fi
if [ -z "$n" ]; then
  echo 'Missing -n (number of samples)' >&2; exit 1
fi
if [ -z "$propheter" ]; then
  echo 'Missing -h (heterogeneity proportion)' >&2; exit 1
fi
if [ -z "$anc" ]; then
  echo 'Missing -a (ancestry)' >&2; exit 1
fi
if [ -z "$baseout" ]; then
  echo 'Missing -b (baseout name)' >&2; exit 1
fi
if [ -z "$pthr" ]; then
  echo 'Missing -p (p-value threshold)' >&2; exit 1
fi

echo parameters S:$nsnp H:$her n:$n h:$propheter a:$anc b:$baseout p:$pthr g:$GW m:$model e:$het_effectsize s:$SAVEFILE

if [[ -s "/transfer/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout.out.prs.fit.txt" ]]; then
  echo "File /transfer/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.$baseout.out.prs.fit.txt exists, nothing to do."
  exit 0
fi

echo Current directory $(pwd), set SCRATCH env var to change.
mkdir -p results/all_iter
mkdir -p tempFiles/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.$baseout && cd tempFiles/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.$baseout && cp -r /treps/* .


if [[ "$SAVEFILE" == "yes" ]]; then
 if [[ "het_effectsize" == "TRUE" ]]; then
   ./simulate_het_effectsize -n $nsnp -r $her -s $n -p $propheter -a $anc -w $GW -o $baseout
 else
   ./simulate -n $nsnp -r $her -s $n -p $propheter -a $anc -w $GW -o $baseout
 fi
else
  if [[ "het_effectsize" == "TRUE" ]]; then
   ./simulate_het_effectsize -n $nsnp -r $her -s $n -p $propheter -a $anc -w $GW -o $baseout  -x
  else
   ./simulate -n $nsnp -r $her -s $n -p $propheter -a $anc -w $GW -o $baseout  -x
 fi
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
    awk '{if(!($1~/'$pop'/)){$NF=-9}}1' 'OFS=\t' $baseout.$sp.ln.fam > $pop.ln.fam;
    plink --bed $baseout.$sp.bed --bim $baseout.$sp.bim --fam $pop.ln.fam --allow-no-sex --assoc --out $pop.ln
    awk '{if(!($1~/'$pop'/)){$NF=-9}}1' 'OFS=\t' $baseout.$sp.target.ln.fam > $pop.target.ln.fam
    plink --bed $baseout.$sp.target.bed --bim $baseout.$sp.target.bim --fam $pop.target.ln.fam --allow-no-sex --assoc --out $pop.target.ln
    
    if [[ "$model" == "all" ]]; then    
      awk '{if(!($1~/'$pop'/)){$NF=-9}}1' 'OFS=\t' $baseout.$sp.pn.fam > $pop.pn.fam;
      plink --bed $baseout.$sp.bed --bim $baseout.$sp.bim --fam $pop.pn.fam --allow-no-sex --assoc --out $pop.pn
      awk '{if(!($1~/'$pop'/)){$NF=-9}}1' 'OFS=\t' $baseout.$sp.target.pn.fam > $pop.target.pn.fam
      plink --bed $baseout.$sp.target.bed --bim $baseout.$sp.target.bim --fam $pop.target.pn.fam --allow-no-sex --assoc --out $pop.target.pn

      awk '{if(!($1~/'$pop'/)){$NF=-9}}1' 'OFS=\t' $baseout.$sp.Wu.fam > $pop.Wu.fam;
      plink --bed $baseout.$sp.bed --bim $baseout.$sp.bim --fam $pop.Wu.fam --allow-no-sex --assoc --out $pop.Wu
      awk '{if(!($1~/'$pop'/)){$NF=-9}}1' 'OFS=\t' $baseout.$sp.target.Wu.fam > $pop.target.Wu.fam
      plink --bed $baseout.$sp.target.bed --bim $baseout.$sp.target.bim --fam $pop.target.Wu.fam --allow-no-sex --assoc --out $pop.target.Wu
    fi
  done
done


./makemeta.R $(pwd) .target $model
for f in `ls metasoft.*.txt`; do cp $f $f.bak; awk 'NR>1' $f | sponge $f; done

for p in AFR AMR EAS EUR SAS TE; do java -jar /Metasoft.jar -input metasoft.$p.ln.txt -pvalue_table /HanEskinPvalueTable.txt -output metasoft.$p.ln.out; done

if [[ "$model" == "all" ]]; then
  for p in AFR AMR EAS EUR SAS TE; do java -jar /Metasoft.jar -input metasoft.$p.pn.txt -pvalue_table /HanEskinPvalueTable.txt -output metasoft.$p.pn.out; done
  for p in AFR AMR EAS EUR SAS TE; do java -jar /Metasoft.jar -input metasoft.$p.Wu.txt -pvalue_table /HanEskinPvalueTable.txt -output metasoft.$p.Wu.out; done
fi

./mvmeta.R $model

#combine the files from all the populations
if [[ "$model" == "all" ]]; then
  for m in pn ln Wu; do
    cat $baseout.AFR.target.$m.fam $baseout.AMR.target.$m.fam $baseout.EAS.target.$m.fam $baseout.EUR.target.$m.fam $baseout.SAS.target.$m.fam > $baseout.target.$m.fam
  done 
else
  cat $baseout.AFR.target.ln.fam $baseout.AMR.target.ln.fam $baseout.EAS.target.ln.fam $baseout.EUR.target.ln.fam $baseout.SAS.target.ln.fam > $baseout.target.ln.fam
fi

#the metasoft.out files have garbage columns at the end (and spelling mistakes in col names)
for f in `ls metasoft.*.out`; do cut -f1-18 $f| sponge $f; done
./processmeta.R -p $pthr -a $anc -s $n -m $model -i $baseout -o snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.$baseout.out

if [[ "$SAVEFILE" == "yes" ]]; then
tar -cvjf snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.$baseout.rundata.tar.bz2  *fam *.log *.qassoc *.out metasoft.* mvmeta* $baseout.*
cp snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.$baseout.rundata.tar.bz2 ../../results/
fi

gzip snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.$baseout.*txt
cp snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.$baseout.*txt.gz ../../results/

cd ../..
rm -r tempFiles/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.$baseout

sleep 5
