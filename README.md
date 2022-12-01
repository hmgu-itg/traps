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
- **S**, the number of SNPs to simulate
- **H**, the heritability of the simulated trait (value between 0 and 1)
- **n**, the number of samples to simulate
- **h**, the proportion of genetic heterogeneity between the populations (value between 0 and 1)
- **a**, the ancetsry distribution of the sample.
    - either a single number giving the proportion of European-ancetsry samples (the four other non-European populations will be sampled in four equal groups)
    - or a comma-separated vector giving the proportions of AFR, AMR, EAS, EUR and SAS populations
- **b**, the baseout name for the output files
- **p**, the significance threshold to select variants to include in the PRS

In addition, traps takes four optional flags:
- **g**, if triggered, 1e6 variants will be simulated in total, leading to the simulation of 1e6-S non-causal variants, i.e. not associated with the phenotype. If not triggered, all the simulated SNPs are considered as causal. Be careful when ticking `-g` as this option can be computationnaly intensive.
- **m**, if triggered, three models are used to simulate the genetic effect sizes of the variants: the point-normal, log-normal and frequency-dependant models. Otherwise, only the log-normal model is used.
- **e**, if triggered, the same causal variants will be present in the different populations but with heterogeneity in genetic effect sizes. If this option is not chosen, heterogeneous variants have non-null genetic effect sizes in only one of the five populations.
- **s**, if triggered, intermediate files are stored in `fullrundata.tar.gz`

# Output of traps
Once completed, traps will provide a `.fit.txt.gz` object containing the following fields:
- **model**: the model used to simulate the genetic effect sizes
- **PRS**: how the PRS was computed:
    - SAS/AMR/EUR/AFR/EAS corresponding to the ancestry-specific PRS
    - TE corresponding to the fixed-effect meta-analysis
    - CHB, ..., STU corresponding to the meta-regression
- **nsnp**: the number of simulated SNPs
- **method**: wether the PRS performances were assessed under the hypothesis that individual-level genetic data are available ('direct') or not ('indirect', relying on the gtx package)
- **cor**: the correlation between the simulated phenotype and the PRS
- **Standard.Error**: the corresponding standard error
- **P**: the p-value associated to the correlation test
- **varex**: the variance of the phenotype explained by the PRS

# Running traps
We recommand to perform at least 10 iterations for each scenario to take into account variability of the simulations. To help the user doing it, we added the wrapper `wrapper_10iterations.sh` which takes the same arguments than `wrapper.sh` and run 10 replicates.  
Two object will be obtained: a `fullrundata.tar.gz` with the final performance estimates and the intermediate files if `-s` is chosen, and an overall results file with the additional fields compared to the output of traps:
- **her**, **samples**, **propheter**, **anc**, **pthr**, **GW** corresponding to the arguments used for the simulations  
- **iter**: the iteration number of the simulation
- **target**: in which subpopulation PRS performance was assessed
- **nsnp_inPRS**: the number of SNPs kept in PRS (passing the corresponding significance threshold)  

# Example of runs
Here are some examples on how to run the simulations with 10 replicates. The `traps/` directory and the `traps.sif` singularity container should be in `mypath`. 
500 causal SNPs, 50% heritabilitiy, 150000 individuals, no heterogeneity, 5 equal groups, using the log-normal model to simulate effect sizes:  
``
wrapper_10iterations.sh -S 500 -H 0.5 -n 150000 -h 0 -a 0.2 -p 5e-8 -P mypath 
``   
300 causal SNPs, 30% heritabilitiy, 150000 individuals, 10% heterogeneity in effect sizes, a given ancestry distribution, using the log-normal model to simulate effect sizes:   
``
wrapper_10iterations.sh -S 300 -H 0.3 -n 150000 -h 0.1 -a 0.1,0.2,0.05,0.5,0.15 -p 5e-8 -P mypath -e
``  
1000 causal SNPs, 20% heritabilitiy, 20000 individuals, 20% heterogeneity in effect sizes, 80% EUR-ancestry individuals, using three models to simulate effect sizes and saving the intermediate files:   
``
wrapper_10iterations.sh -S 1000 -H 0.2 -n 20000 -h 0.2 -a 0.8 -p 5e-8 -P mypath -m -s
``  

# References
**1000Genome project**: Sudmant, P. H. et al. An integrated map of structural variation in 2,504 human genomes. Nature 526, 75-81, doi:10.1038/nature15394 (2015).  
**Log-normal model**: O'Connor, L. J. The distribution of common-variant effect sizes. Nat Genet 53, 1243-1249, doi:10.1038/s41588-021-00901-3 (2021).  
**Frequency-dependent model**: Wu, M. C. et al. Rare-variant association testing for sequencing data with the sequence kernel association test. Am J Hum Genet 89, 82-93, doi:10.1016/j.ajhg.2011.05.029 (2011).  
**TAMR method**: Magi, R. et al. Trans-ethnic meta-regression of genome-wide association studies accounting for ancestry increases power for discovery and improves fine-mapping resolution. Hum Mol Genet 26, 3639-3650, doi:10.1093/hmg/ddx280 (2017).  
**Package gaston**: Dandine-Roulland C, Perdry H (2015) The Use of the Linear Mixed Model in Human Genetics. Hum Hered 80:196–206. https://doi.org/10.1159/000447634.  
**Package gtx**: gtx (R CRAN, 2019).
