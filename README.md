# mitoimpute

A nextflow pipeline for generating genetic scores for mtDNA-encoded transcript abundance using FUSION.

## Pre-requististes:

Nextflow

Docker/Singularity

## Usage:

```nextflow run main.nf --name "OUTNAME" --bfile "PLINKSTUB" --eqtl "EQTLFILE" --pvalue "PVALUE" --pheno "PHENONAME" --phenofile "PHENOFILE" --matchBim "MATCHEDBIM" --useMito "0/1" --useCovars "0/1" --liftmatch "0/hg19tohg20/hg20tohg19" --liftcompare "0/hg19tohg20/hg20tohg19" --model "lasso/enet/blup/top1" --targetbfile "TARGETPLINKSTUB" --targetpheno "TARGETPHENONAME" --targetphenofile "TARGETPHENOFILE" --targetcovars "TARGETCOVARS" --useRef "hg19/hg20/NA" -profile singularity/docker```

## Options:

```--name "OUTNAME"``` Name for output

```--bfile "PLINKSTUB"``` Full path to PLINK genetic file in binary format. The filename should be the name without the .bed/.bim/.fam.

```--eqtl "EQTLFILE"``` Full path to eQTL results file in PLINK format

```--pvalue "PVALUE"``` P value threshold for selecting genetic variants from the eQTL file

```--phenofile "PHENOFILE"``` Full path to a PLINK phenotype file, containing variable and covariates to be used in model building

```--pheno "PHENONAME"``` The name of the phenotype to be used in the PLINK phenotype file

```--matchBim "MATCHEDBIM"``` Full path to PLINK .bim file that you want to select SNPs from.  If you don't want to match to specific genetic variants, then supply the path to the .bim file you use for PLINKSTUB above

```--useMito "0/1"``` Include genetic variation in mtDNA when building model (1), or not (0).

```--useCovars "0/1"``` Regress out covariates against the phenotype for model building (1) or not (0).  Coavriates used will be those in the PLINK phenotype file between 'combined_s20' and 'RNA' columns (automatically generated with mitogwas pipeline).

```--liftmatch "0/hg19tohg20/hg20tohg19"``` If selecting genetic variants to include in the model (--matchBim above) state whether liftover is needed (hg19tohg20 where original data is hg20 and matched data is hg19, hg19tohg20 for the reverse, and 0 if liftover not needed).

```--liftcompare "0/hg19tohg20/hg20tohg19"``` For comparison of the genetic score model to target data (below) state whether liftover needed (hg19tohg20 where original data is hg20 and target data is hg19, hg19tohg20 for the reverse, and 0 if liftover not needed).

```--model "lasso/enet/blup/top1"``` ML method used for genetic score model

To compare predicted and actual transcript abundance values, supply target files with the following (if not wanted, just supply paths to data used to build models as above):

```--targetbfile "TARGETPLINKSTUB"``` Full path to target PLINK genetic file in binary format. The filename should be the name without the .bed/.bim/.fam. 

```--targetpheno "TARGETPHENONAME"``` The name of the phenotype to be used in the target PLINK phenotype file

```--targetphenofile "TARGETPHENOFILE"``` Full path to a target PLINK phenotype file, containing variable and covariates to be used in model building

```--targetcovars "TARGETCOVARS"``` A list of comma seperated covariates to be regressed against the target phenotype (optional)

```--useRef "hg19/hg20/0"``` Filter genetic variants to be used in model building to those on the Hapmap reference panel (either hg20 or hg19) or not (0).

```-profile singularity/docker``` use either singularity or docker for the machine image
