# mitoimpute

A nextflow pipeline for generating genetic scores for mtDNA-encoded transcript abundance using FUSION.

## Pre-requististes:

Nextflow

Docker/Singularity

## Usage:

```nextflow run main.nf --name "OUTNAME" --bfile "PLINKSTUB" --eqtl "EQTLFILE" --pvalue "PVALUE" --pheno "PHENONAME" --phenofile "PHENOFILE" --matchBim "MATCHEDBIM" --useMito "0/1" --useCovars "0/1" --liftmatch "0/1" --liftcompare "0/1" --model "lasso/enet/blup/top1" --targetbfile "TARGETPLINKSTUB" --targetpheno "TARGETPHENONAME" --targetphenofile "TARGETPHENOFILE" --targetcovars "TARGETCOVARS" --useRef "hg19/hg20/NA" -profile singularity/docker```

## Options:

--rnaDir "RNAseqDIR" : Full path to directory containing gene count files.  This directory should also contain a subdirectory called RNAseQC, which contains output of RNAseQC run on STAR alignment files.

--bed "BEDFILE" : Full path to plink format bed file, containing genetic variant data for samples to be used in the analysis.

--bim "BIMFILE" : Full path to plink format bim file, containing genetic variant data for samples to be used in the analysis.

--fam "FAMFILE" : Full path to plink format fam file, containing sample names (these should match those supplied in the inforsheet DNA column) to be used in the analysis.

--infoSheet "INFOSHEET": Full path to inforsheet file. This file is plain text with header "RNA DNA", followed by the names of any covariates to be used in the eQTL analysis. Below this there should be one line per sample to be included in the analysis, with the corresponding DNA (in fam file) and RNA (stub of STAR output) code for the sample, followed by values of any covarites to be used in the analysis.

--gtfFile "GTFFILE" : Full path to GTF file used for STAR aignment

--dataName "NAME" : Name for output

-profile singularity/docker : use either singularity or docker for the machine image
