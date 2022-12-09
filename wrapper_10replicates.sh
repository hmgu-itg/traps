#!/bin/bash
model='ln'
GW='FALSE'
SAVEFILE='no'
het_effectsize='FALSE'

while getopts "S:H:n:h:a:p:T:P:gmes" opt; do
  case $opt in
    S) nsnp="$OPTARG" ;;
    H) her="$OPTARG" ;;
    n) n="$OPTARG" ;;
    h) propheter="$OPTARG" ;;
    a) anc="$OPTARG" ;;
    p) pthr="$OPTARG" ;;
    T) traps_path="$OPTARG" ;; 
    P) mypath="$OPTARG" ;;
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
if [ -z "$pthr" ]; then
  echo 'Missing -p (p-value threshold)' >&2; exit 1
fi
if [ -z "$traps_path" ]; then
  echo 'Missing -T (Path to bind to container)' >&2; exit 1
fi
if [ -z "$mypath" ]; then
  echo 'Missing -P (Path to write the output files)' >&2; exit 1
fi

command_params="-S $nsnp -H $her -n $n -h $propheter -a $anc -p $pthr"

if [ "$model" == "all" ]; then
  command_params=$(echo $command_params " -m")
fi
if [ "$GW" == "TRUE" ]; then
  command_params=$(echo $command_params " -g")
fi
if [ "$het_effectsize" == "TRUE" ]; then
  command_params=$(echo $command_params " -e")
fi
if [ "$SAVEFILE" == "yes" ]; then
  command_params=$(echo $command_params " -s")
fi

for i in {1..10}; do
  command_torun=$(echo $command_params " -b iter$i")
  FILE=$mypath/results/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.model$model.hetES$het_effectsize.iter$i.out.prs.fit.txt.gz
  if [[ ! -f "$FILE" ]]; then
    eval $(singularity exec -B $mypath -B $traps_path $traps_path/traps.sif $traps_path/traps/wrapper.sh $command_torun)
  fi
done

cat <(echo -e "nsnp\ther\tsamples\tpropheter\tanc\tpthr\tGW\titer\tmodel\tPRS\ttarget\tnsnp_inPRS\tmethod\tcor\tStandard.Error\tP\tvarex") <(for i in {1..10}; do zcat $mypath/results/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.model$model.hetES$het_effectsize.iter$i.out.prs.fit.txt.gz |tail -n+2 | while read l; do echo -e "$nsnp\t$her\t$n\t$propheter\t$anc\t$pthr\t$GW\t$i\t$l";done;done) | gzip > $mypath/results/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.model$model.hetES$het_effectsize.fullrun.out.gz
tar -cvzf $mypath/results/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.model$model.hetES$het_effectsize.fullrundata.tar.gz $mypath/results/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.model$model.hetES$het_effectsize.iter* && rm -f $mypath/results/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.model$model.hetES$het_effectsize.iter*
mv -f $mypath/results/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.model$model.hetES$het_effectsize.fullrundata.tar.gz $mypath/results/all_iter/
mv -f $mypath/results/snp$nsnp.her$her.heter$propheter.anc$anc.p$pthr.GW$GW.model$model.hetES$het_effectsize.fullrun.out.gz $mypath/results/all_iter/