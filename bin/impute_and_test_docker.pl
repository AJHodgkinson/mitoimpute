use strict; 
use Statistics::R;

my $R = Statistics::R->new();
my $make;

my $covars="";

my $name=$ARGV[0]; #Name for output files
my $rweights=$ARGV[1]; #R results file from Fusion
my $model=$ARGV[2]; #Select the model for testing
my $plink_stub=$ARGV[3]; #Plink binary stub for target population
my $pheno_file=$ARGV[4]; #Phenotype file
my $pheno=$ARGV[5]; #Phenotype to be tested against
my $lift=$ARGV[6]; #0 if not required, hg19tohg20 or hg20tohg19 (fusion data genome on the left, target data genome on the right)
my $selected=$ARGV[7]; #Bim file from selected SNPs so number can be recorded
$covars=$ARGV[8]; #Covariates if required (seperated by commas)

my $report_header="Name";
my $report="$name";

#Open BIM and count selected SNPs:
my $starting=0;
open (STARTING, "$selected") || die "Unable to open weights file to read: $!\n";
while (<STARTING>) {
  $starting++;
}
close (STARTING);

$report_header.="\tStartingVariants";
$report.="\t$starting";

#Open Fusion model output
$make="data<-load(\"$rweights\")";
$R->send($make);

##Store models tested an accompanying statistics
my @models=(); my @rsq=(); my @pval=(); 
my $performance = $R->get('cv.performance');
my $state=0;
foreach my $element (@{$performance}) {
  push @models, $element if (($state==0)&&(!($element=~/rsq|pval/)));
  push @rsq, $element if (($state==1)&&(!($element=~/rsq|pval/)));
  push @pval, $element if (($state==2)&&(!($element=~/rsq|pval/)));
  $state++ if ($element=~/rsq|pval/);
}

##Calculate position of selected model for downstream marker selection
my $position;
for (my $i=0;$i<@models;$i++) {
  $position=$i if ($models[$i] eq $model);
}

my $start=@models;

#Store weights and SNPs for selected model
$make="write.table(wgt.matrix, file = \"${name}.table.weights.txt\", sep = \"\t\")";
$R->send($make);
my @weights=();
open (WEIGHTS, "${name}.table.weights.txt") || die "Unable to open weights file to read: $!\n";
while (<WEIGHTS>) {
  chomp;
  my @array=split;
  foreach my $item (@array) {
    $item=~s/\"//g;
    push @weights, $item;
  }
}
close (WEIGHTS);

my $snp_count=0;

my %weights=();
for (my $i=$start;$i<@weights;$i=$i+($start+1)) {
  $weights{$weights[$i]}=$weights[$i+1+$position];
  $snp_count++ if ($weights[$i+1+$position] != 0);
}


#Store SNP data for selected model
$make="write.table(snps, file = \"${name}.table.snps.txt\", sep = \"\t\", row.names = FALSE)";
$R->send($make);
my @snps=();
open (SNPS, "${name}.table.snps.txt") || die "Unable to open SNPs file to read: $!\n";
while (<SNPS>) {
  chomp;
  my @array=split;
  foreach my $item (@array) {
    $item=~s/\"//g;
    push @snps, $item;
  }
}
close (SNPS);

my %coordinates=(); my %weight_allele=(); my %snp_alleles=(); my %collect=();
#print "SNPS:\n";
for (my $i=6;$i<@snps;$i=($i+6)) {
  my $chr=$snps[$i]; $chr=~s/MT/26/; $chr=~s/M/26/;
  $coordinates{$snps[$i+1]}=$chr."_".$snps[$i+3];
  $weight_allele{$snps[$i+1]}=$snps[$i+4];
  $snp_alleles{$snps[$i+1]}=$snps[$i+4]."_".$snps[$i+5];
  my $tag=$chr."_".$snps[$i+3]."_".$snps[$i+4]."_".$snps[$i+5];
  $collect{$tag}=$snps[$i+1];
}
$report_header.="\tVariants";
$report.="\t$snp_count";

#Report CV statistics
my $hsq = $R->get('hsq');
my $middle="NA";
if ($hsq==1) {
  $middle=1;
}
if (!($hsq==1)) {
  my @hsq=@{$hsq};
  my $upper=$hsq[0]; my $lower=$hsq[1]; $middle=$lower+(($upper-$lower)/2);
}

$report_header.="\tHerit";
$report.="\t$middle";

my $hsqpv = $R->get('hsq.pv');
my $tot = $R->get('N.tot');

$report_header.="\tHerit_P\tIndividuals";
$report.="\t$hsqpv\t$tot";

for (my $i=0;$i<@models;$i++) {
  $report_header.="\tModel\t5F-Rsq\t5F-P" if ($models[$i] eq $model);
  $report.="\t$models[$i]\t$rsq[$i]\t$pval[$i]" if ($models[$i] eq $model);
}

#If liftover is required, create BED file for selected SNPs, run liftover, and store SNPs for later comparison
my %liftover=();

open (LIFT, ">$name.lift.bed") || die "Unable to open lift file to write to: $!\n";
for (my $i=6;$i<@snps;$i=($i+6)) {
  my $chr=$snps[$i]; $chr=~s/MT/26/; $chr=~s/M/26/;
  my $pos=$snps[$i+3]; my $pos1=$pos+1;
  print LIFT "chr$chr\t$pos\t$pos1\t${chr}_$pos\n"; #Add 'chr' so liftover will work
}

close (LIFT);

#system `/scratch/grp/hodgkinsonlab/programs_urgent/liftOver $name.lift.bed /scratch/grp/hodgkinsonlab/shared/GTEx_Data/Models/hg38ToHg19.over.chain.gz $name.lifted.bed $name.not_lifted.bed` if ($lift=~/hg19tohg20/);
#system `/scratch/grp/hodgkinsonlab/programs_urgent/liftOver $name.lift.bed /scratch/grp/hodgkinsonlab/shared/GTEx_Data/Models/hg19ToHg38.over.chain.gz $name.lifted.bed $name.not_lifted.bed` if ($lift=~/hg20tohg19/);
system `liftOver $name.lift.bed /scratch/grp/hodgkinsonlab/shared/GTEx_Data/Models/hg38ToHg19.over.chain.gz $name.lifted.bed $name.not_lifted.bed` if ($lift=~/hg19tohg20/);
system `liftOver $name.lift.bed /scratch/grp/hodgkinsonlab/shared/GTEx_Data/Models/hg19ToHg38.over.chain.gz $name.lifted.bed $name.not_lifted.bed` if ($lift=~/hg20tohg19/);
system `grep chr26 $name.lift.bed >> $name.lifted.bed` if ($lift=~/hg/); #Add Mito chromosome, as this will not have been lifted over (but co-ordinates are the same between hg19 and hg20)

system `cp $name.lift.bed $name.lifted.bed` if (!($lift=~/hg/)); #If no liftover reuqured, just copy the file
system `touch $name.not_lifted.bed` if (!($lift=~/hg/)); #If no liftover reuqured, just copy the file

open (LIFTED, "$name.lifted.bed") || die "Unable to open $name.lifted.bed to read: $!\n";
while (<LIFTED>) {
  my @array=split;
  my $tag=$array[0]."_".$array[1]; $tag=~s/chr//;
  $liftover{$tag}=$array[3];
}
close (LIFTED);


#####Copy over target SNP file and convert SNP names to be compatible with Fusion output model

system `cp $plink_stub.bed $name.full.bed`;
system `cp $plink_stub.bim $name.full.bim`;
system `cp $plink_stub.fam $name.full.fam`;

open (BIMFULL, "$name.full.bim") || die "Unable to open copied BIM file $name.full.bim to read: $!\n";
open (BIMHOLD, ">$name.hold.bim") || die "Unable to open copied BIM file $name.full.bim to read: $!\n";

my @markers=();
while (<BIMFULL>) {
  chomp;
  my @array=split;
  my $chr=$array[0]; $chr=~s/chr//; $chr=~s/MT/26/; $chr=~s/M/26/;
  print BIMHOLD "$chr\t${chr}_$array[3]_$array[5]_$array[4]\t$array[2]\t$array[3]\t$array[4]\t$array[5]\n";
  my $marker="${chr}_$array[3]_$array[5]_$array[4]";
  push @markers, $marker;
}

close (BIMFULL);
close (BIMHOLD);

system `mv $name.hold.bim $name.full.bim`;

#Loop through BIM file for taregt data and collect SNP information for SNPs in model, and also collect weight and allele information from these sites:
open (BIM, "$name.full.bim") || die "Unable to open target BIM file to read: $!\n";
open (COLLECT, ">$name.snps.txt") || die "Unable to open target BIM file to read: $!\n";
my @feature_weights=(); my @weight_allele=();

while (<BIM>) {
  my @array=split;
  my $chrpos=$array[0]."_".$array[3];
  if (!($lift=~/hg/)) {
    my $tag1=$chrpos."_".$array[4]."_".$array[5];
    my $tag2=$chrpos."_".$array[5]."_".$array[4];
    if ((exists $collect{$tag1})||(exists $collect{$tag2})) {
      print COLLECT "$array[1]\n";
      push  @feature_weights, $weights{$collect{$tag1}} if (exists $weights{$collect{$tag1}});
      push  @feature_weights, $weights{$collect{$tag2}} if (exists $weights{$collect{$tag2}});
      push @weight_allele, $weight_allele{$collect{$tag1}} if (exists $weight_allele{$collect{$tag1}});
      push @weight_allele, $weight_allele{$collect{$tag2}} if (exists $weight_allele{$collect{$tag2}});
    }
  }
  if ($lift=~/hg/) { #Use liftover conversion if required
    if (exists $liftover{$chrpos}) {
      my $chrpos1=$liftover{$chrpos};
      my $tag1=$chrpos1."_".$array[4]."_".$array[5];
      my $tag2=$chrpos1."_".$array[5]."_".$array[4];
      if ((exists $collect{$tag1})||(exists $collect{$tag2})) {
	print COLLECT "$array[1]\n";
	push  @feature_weights, $weights{$collect{$tag1}} if (exists $weights{$collect{$tag1}});
	push  @feature_weights, $weights{$collect{$tag2}} if (exists $weights{$collect{$tag2}});
	push @weight_allele, $weight_allele{$collect{$tag1}} if (exists $weight_allele{$collect{$tag1}});
	push @weight_allele, $weight_allele{$collect{$tag2}} if (exists $weight_allele{$collect{$tag2}});
      }
    }
  }  
}
close (BIM);
close (COLLECT);

my $snps_used=@feature_weights;

$report_header.="\tSNPs_Used";
$report.="\t$snps_used";

#####Filter just required SNPs from target file

system `plink --bfile $name.full --extract $name.snps.txt --recode --out $name.run`;

#####Impute mitochondria features to target file

#print "@feature_weights\n\n@weight_allele\n\n";

my %impute=();
open (PED, "$name.run.ped") || die "Unable to open target BIM file to read: $!\n";

while (<PED>) {
  my @array=split;
  $impute{$array[0]}=0;
  my $marker=0;
  for (my $i=6;$i<@array;$i=$i+2) {
    $impute{$array[0]}+=$feature_weights[$marker] if ($array[$i] eq $weight_allele[$marker]);
    $impute{$array[0]}+=$feature_weights[$marker] if ($array[$i+1] eq $weight_allele[$marker]);
    $marker++;
  }
}
close (PED);

#####Open and store selected phenotype data for comparison

my %data=(); 

if (!($covars=~/\w/)) { #Just get phenotype if no covars included
  open (PHENO1, "$pheno_file") || die "Unable to open pheno file $pheno_file to read: $!\n";
  my @head=();
  while (<PHENO1>) {
    my @array=split;
    if ($_=~/^FID/) {
      @head=@array;
    }
    else {
      for (my $i=0;$i<@array;$i++) {
	if ($head[$i] eq $pheno) {
	  $data{$array[0]}=$array[$i];
	}
      }
    }
  }
  close (PHENO1);
}

if ($covars=~/\w/) { #If covariates included, regress these out first then store data for chosen phenotype
  my @covars=split(/\,/,$covars);
  my %covar=();
  foreach my $c (@covars) {
    $covar{$c}=0;
  }
  $covars=~s/\,/\+/g;
  
  open (PHENO2, "$pheno_file") || die "Unable to open pheno file $pheno_file to read: $!\n";
  my @inds=(); my @head=(); my $exp; my %covars=(); my %miss=();
  
  my $line=0; my @head=();
  while (<PHENO2>) {
    my @array=split;
    if ($line==0) {
      @head=@array;
    }
    if ($line>0) {
      push @inds, $array[0];
      $data{$array[0]}="NA";
      for (my $i=0;$i<@array;$i++) {
	if ($head[$i] eq $pheno) {
	  $exp.=",".$array[$i];
	  $miss{$array[0]}=0 if ($array[$i] eq "NA");
	}
	if (exists $covar{$head[$i]}) {
	  my $number = "NA";
	  $number = sprintf "%.4f", $array[$i] if ($array[$i]=~/\d/);
	  $covars{$head[$i]}.=",".$number;
	  $miss{$array[0]}=0 if ($array[$i] eq "NA");
	}
      }
    }
    $line++;
  }
  close (PHENO2);
  
  
  while (my ($a,$b) = each %covars) {
    my $dat=$b; $dat=~s/,//;
    $make="$a<-c($dat)";
    $R->run($make);
  }
  
  $exp=~s/,//;
  $make="gene<-c($exp)";
  $R->run($make);
  
  $make="res<-resid(lm(gene~$covars))";
  $R->send($make);
  my $res = $R->get('res');
  my @final=();
  foreach my $element (@{$res}) {
    if ($element=~/\./) {
      push @final, $element;
    }
  }
  
  my $bump=0; my $iii="";
  for (my $i=0;$i<@inds;$i++) {
    if (!(exists $miss{$inds[$i]})) {
      $data{$inds[$i]}=$final[$i-$bump];
    }
    if (exists $miss{$inds[$i]}) { #Check of this fails if any of the covars are NA too - if they are, add in NA hash to check for this as covars are extracted.
      $data{$inds[$i]}="NA";
      $bump++;
    }
  }
}


#####Run and report association test

my $imputed; my $real; my $inds_used=0;
while (my ($a,$b) = each %data) {
  if (exists $impute{$a}) {
    $imputed.=",".$impute{$a};
    $real.=",".$b;
    $inds_used++;
  }
}

$report_header.="\tInds_Used";
$report.="\t$inds_used";

$imputed=~s/,//; $real=~s/,//;

$make="imputed<-c($imputed)\nreal<-c($real)";
$R->send($make);

$make="est<-cor.test(real,imputed)\$estimate";
$R->send($make);
my $est = $R->get('est');

$make="pval<-cor.test(real,imputed)\$p.value";
$R->send($make);
my $pval = $R->get('pval');

$report_header.="\tCor_R\tCor_P";
$report.="\t$est\t$pval";

open (RESULTS, ">$name.results.txt") || die "Unable to open results file to write to: $!\n";
print RESULTS "$report_header\n$report\n";
close (RESULTS);

