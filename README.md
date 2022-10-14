# traps
TRans-Ancestry PRS Simulator
traps has been developped to assess the performances of the PRS based on simulations: an ancestry-specific PRS, a PRS based on a trans-ancestry fixed-effect model, and a PRS based on the trans-ancestry meta-regression (TAMR).  
The simulations are based on the 1000Genome project allelic frequencies and on simulated effect sizes for the variants. 
20 subpopulations (4 for each population AFR, AMR, SAS, EAS and EUR) were randomly selected. Performances are computed in each subpopulation as a measure of fit between PRS and simulated phenotypes.  

# Overview of the simulation procedure

![Treps_procedure](https://user-images.githubusercontent.com/26084630/195874633-32ff26e9-791a-4dbb-96f9-bf9000e344a1.jpg)


# Options of the traps procedure
The main script to run the whole traps procedure is `wrapper.sh` that will call successively function to simulate the data, compute the meta-analysis, compute the PRS and assess their performance.  
The following arguments need to be specified by the user:
- S, the number of SNPs to simulate
- H, the heritability of the simulated trait
- n, the number of samples to simulate
- h, the proportion of genetic heterogeneity between the populations
- a, the ancetsry distribution of the sample.
    - either a single number giving the proportion of European-ancetsry samples (the four other non-European populations will be sampled in four equal groups)
    - or a comma-separated vector giving the proportions of AFR, AMR, EAS, EUR and SAS populations
- b, the baseout name for the output files
- p, the significance threshold to select variants to include in the PRS

In addition, traps takes four optional flags:
- g, if triggered, 1e6 variants will be simulated in total, leading to the simulation of 1e6-S non-causal variants, i.e. not associated with the phenotype. If not triggered, all the simulated SNPs are considered as causal. Be careful when ticking `-g` as this option can be computationnaly intensive.
- m, if triggered, three models are used to simulate the genetic effect sizes of the variants: the point-normal, log-normal and frequency-dependant models. Otherwise, only the log-normal model is used.
- e, if triggered, the same causal variants will be present in the different populations but with heterogeneity in genetic effect sizes. If this option is not chosen, heterogeneous variants have non-null genetic effect sizes in only one of the five populations.
- s, if triggered, intermediate files are stored in `fullrundata.tar.gz`

# Output of traps
Once completed, traps will provide a `fullrundata.tar.gz` object with the final performance estimates and the intermediate files if `-s` is chosen, and the results files with the following fields:
- model: the model used to simulate the genetic effect sizes
- PRS: how the PRS was computed:
    - SAS/AMR/EUR/AFR/EAS corresponding to the ancestry-specific PRS
    - TE corresponding to the fixed-effect meta-analysis
    - CHB, ..., STU corresponding to the meta-regression
- nsnp: the number of simulated SNPs
- method: wether the PRS performances were assessed under the hypothesis that individual-level genetic data are available ('direct') or not ('indirect', relying on the gtx package)
- cor: the correlation between the simulated phenotype and the PRS
- Standard.Error: the corresponding standard error
- P: the p-value associated to the correlation test
- varex: the variance of the phenotype explained by the PRS
