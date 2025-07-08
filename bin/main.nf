#!/usr/bin/env nextflow

params.name = ''
params.bfile = ''
params.eqtl = ''
params.pvalue = ''
params.pheno = ''
params.phenofile = ''
params.matchBim = ''
params.useMito = ''
params.useCovars = ''
params.useRef = ''
params.liftmatch = ''
params.liftcompare = ''
params.model = ''
params.targetbfile = ''
params.targetphenofile = ''
params.targetpheno = ''
params.targetcovars = ''

params.outDir = './results'

process makeFusion {

    publishDir "$params.outDir/", pattern: '*selected.bim', mode: 'copy'

    input:
    val name
    val bfile
    path eqtl
    val pvalue
    val pheno
    path phenofile
    path matchBim
    val useMito
    val useCovars
    val useRef
    val liftmatch

    output:
    path "${name}.selected.bed"
    path "${name}.selected.bim"
    path "${name}.selected.fam"
    
    script:   
    """
    perl ${baseDir}/bin/make_fusion_matching_docker.pl $name $bfile $eqtl $pvalue $pheno $phenofile $matchBim $useMito $useCovars $useRef $liftmatch
    """
}

process runFusion {

    publishDir "$params.outDir/", pattern: '*wgt*', mode: 'copy'

    input:
    val name
    path bed
    path bim
    path fam

    output: 
    path "${name}.wgt.RDat"

    script:
    """
    ln -s ./ output
    Rscript /fusion_twas-master/FUSION.compute_weights.R --bfile ${name}.selected --hsq_p 1 --tmp tmp.${name}.txt --out $name --models enet,blup,top1,lasso --PATH_plink /bin/plink --PATH_gcta /fusion_twas-master/gcta_nr_robust --PATH_gemma /fusion_twas-master/gemma
    """
}

process imputeTest {

    publishDir "$params.outDir/", pattern: '*results.txt', mode: 'copy'

    input:
    val name
    path weights
    val model
    val target_bfile
    path phenofile
    val pheno
    val liftcompare
    path starting
    val covars

    output: 
    path "${name}.results.txt"

    script:   
    """
    echo "hello3"
    perl ${baseDir}/bin/impute_and_test_docker.pl $name $weights $model $target_bfile $phenofile $pheno $liftcompare $starting $covars
    """
}

workflow {   
    makefusion_ch = makeFusion(params.name,params.bfile,params.eqtl,params.pvalue,params.pheno,params.phenofile,params.matchBim,params.useMito,params.useCovars,params.useRef,params.liftmatch)
    runfusion_ch = runFusion(params.name, makefusion_ch)
    testfusion_ch = imputeTest(params.name,runfusion_ch,params.model,params.targetbfile,params.targetphenofile,params.targetpheno,params.liftcompare,makefusion_ch.getAt(1),params.targetcovars)
}












