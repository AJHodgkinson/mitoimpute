# mitoimpute

A nextflow pipeline for generating genetic scores for mtDNA-encoded transcript abundance using FUSION.

## Pre-requististes:

Nextflow

Docker/Singularity

## Usage:

```nextflow run main.nf --name "OUTNAME" --bfile "PLINKSTUB" --eqtl "EQTLFILE" --pvalue "PVALUE" --pheno "PHENONAME" --phenofile "PHENOFILE" --matchBim "MATCHEDBIM" --useMito "0/1" --useCovars "0/1" --liftmatch "0/1" --liftcompare "0/1" --model "lasso/enet/blup/top1" --targetbfile "TARGETPLINKSTUB" --targetpheno "TARGETPHENONAME" --targetphenofile "TARGETPHENOFILE" --targetcovars "TARGETCOVARS" --useRef "hg19/hg20/NA" -profile singularity/docker```

## Options:

```--name "OUTNAME"``` Name for output

```--bfile "PLINKSTUB"``` Full path to PLINK genetic file in binary format. The filename should be the name without the .bed/.bim/.fam.

```--eqtl "EQTLFILE"``` Full path to eQTL results file in PLINK format

```--pvalue "PVALUE"``` P value threshold for selecting genetic variants from the eQTL file

```--phenofile "PHENOFILE"``` Full path to a PLINK phenotype file, containing variable and covariates to be used in model building

```--pheno "PHENONAME"``` The name of the phenotype to be used in the PLINK phenotype file

```--matchBim "MATCHEDBIM"``` Full path to PLINK .bim file that you want to select SNPs from.  If you don't want to match to specific SNPs, then supply the path to the .bim file you use for PLINKSTUB above

```--useMito "0/1"``` Include genetic variation in mtDNA when building model (1), or not (0).

```--useCovars "0/1"``` Regress out covariates against the phenotype for model building (1) or not (0).  Coavriates used will be those in the PLINK phenotype file between 'combined_s20' and 'RNA' columns (automatically generated with mitogwas pipeline).

```--liftmatch "0/1"```

```--liftcompare "0/1"```

```--model "lasso/enet/blup/top1"```

```--targetbfile "TARGETPLINKSTUB"```

```--targetpheno "TARGETPHENONAME"```

```--targetphenofile "TARGETPHENOFILE"```

```--targetcovars "TARGETCOVARS"```

```--useRef "hg19/hg20/NA"```

```-profile singularity/docker``` use either singularity or docker for the machine image
