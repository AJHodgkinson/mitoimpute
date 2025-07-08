use strict;
use Statistics::R;

my $R = Statistics::R->new();
my $make;

my $use_mito=0; #Deafult to not use mitochondrial SNPs
my $use_covars=0; #Deafult is to not regress out covariates
my $use_ref=0; #Default is to not use the reference panel

my $name=$ARGV[0]; #Name of output stub
my $bfile=$ARGV[1]; #PLINK bfile stub for data
my $eqtl=$ARGV[2]; #PLINK linear results files for gene of interest (can provide ALL or EUR within this)
my $pvalue=$ARGV[3]; #P value cut off to select SNPs
my $pheno=$ARGV[4]; #Name of phenotype to use
my $pheno_file=$ARGV[5]; #Path to plink phenotype file if $use_covars==1
my $match_bim=$ARGV[6]; #Path to a BIM file to match SNPs to (either Reference LD panel or data compare)

$use_mito=$ARGV[7]; #Include MT genome (1) or not (0)
$use_covars=$ARGV[8]; #Regress out covars (1) or not (0)
$use_ref=$ARGV[9]; #Use LD Ref hg19,hg20 or not (0) or clump (clump)

my $lift=$ARGV[10]; #Flag for converting co-ordinates, 0 if not needed, hg19tohg20 or hg20tohg19 (fusion data genome on the left, target data genome on the right)

#First, collect eQTL data for each SNP and collect SNP tags:

my %eqtlp=();

open (EQTL, "gunzip -c $eqtl |") || die "Unable to open $name.temp.eqtl.gz to read: $!\n"; #open eqtl file
while (<EQTL>) {
  chomp;
  my @array=split;
  if (!($_=~/BETA/)) {
    my $chr=$array[0]; $chr=~s/chr//; $chr=~s/MT/26/; $chr=~s/M/26/; #Remove "chr" and convert Mt to 26
    my $tag=$chr."_".$array[2]."_".$array[3]; #create tag with chr_pos_A1
    $eqtlp{$tag}=$array[11]; #Store P value
  }
}

close (EQTL);

my %flip=();  my %source=();
open (BB, "$bfile.bim") || die "Unable to open clump file $name.clumped to read: $!\n";
while (<BB>) {
  my @array=split;
  my $chr1=$array[0]; $chr1=~s/chr//;
  my $chr=$array[0]; $chr=~s/chr//; $chr=~s/MT/26/; $chr=~s/M/26/;
  my $marker="${chr}_$array[3]_$array[5]_$array[4]";
  $flip{$array[1]}=$marker;
  $source{$marker}=$array[1];
}
close (BB);

#Second, store the SNPs that are in the matching file (can use target file, or LD reference file). Liftover if required:
my %match=();

open (MATCH, "$match_bim") || die "Unable to open $match_bim to rea1: $!\n";
open (LIFT, ">$name.lift.bed") || die "Unable to open lift file to write to: $!\n";

while (<MATCH>) {
  my @array=split;
  my $chr1=$array[0]; $chr1=~s/chr//;
  my $chr=$array[0]; $chr=~s/chr//; $chr=~s/MT/26/; $chr=~s/M/26/;
  my $pos=$array[3]; my $pos1=$pos+1;
  my $marker="${chr}_$array[3]_$array[5]_$array[4]";
  print LIFT "chr$chr1\t$pos\t$pos1\t${chr}_$pos\t$array[4]\t$array[5]\t$array[1]\n"; #Create liftover file
  $match{$marker}=0 if (($chr=~/26/)&&($use_mito==1)); #Just add all MT SNPs to the match hash if MT SNPs required, as they are not in liftover file
}

close (LIFT);

system `liftOver $name.lift.bed /scratch/grp/hodgkinsonlab/shared/GTEx_Data/Models/hg38ToHg19.over.chain.gz $name.lifted.bed $name.not_lifted.bed` if ($lift=~/hg20tohg19/);
system `liftOver $name.lift.bed /scratch/grp/hodgkinsonlab/shared/GTEx_Data/Models/hg19ToHg38.over.chain.gz $name.lifted.bed $name.not_lifted.bed` if ($lift=~/hg19tohg20/);
system `cp $name.lift.bed $name.lifted.bed` if (!($lift=~/hg/)); #If no liftover reuqured, just copy the file
system `touch $name.not_lifted.bed` if (!($lift=~/hg/)); #If no liftover reuqured, just copy the file

open (LIFTED, "$name.lifted.bed") || die "Unable to open $name.lifted.bed to read: $!\n"; #Read liftover file and store matching SNPs (if reference genome of original data)
while (<LIFTED>) {
  my @array=split;
  my $tag=$array[0]."_".$array[1]."_".$array[5]."_".$array[4]; $tag=~s/chr//;
  $match{$tag}=$array[6];
}
close (LIFTED);

my $reference_file="/scratch/grp/hodgkinsonlab/shared/GTEx_Data/Models/LDREF.ALL.hg38.bim";
$reference_file="/scratch/grp/hodgkinsonlab/shared/GTEx_Data/Models/LDREF.ALL.bim" if ($use_ref eq "hg19");

if ($use_ref=~/hg/) {
  my %panel=();
  
  open (PANEL, "$reference_file") || die "Unable to open reference panel file to read $reference_file\n";
  while (<PANEL>) {
    chomp;
    my @array=split;
    my $tag1=$array[0]."_".$array[3]."_".$array[5]."_".$array[4];
    my $tag2=$array[0]."_".$array[3]."_".$array[4]."_".$array[5];
    $panel{$tag1}=$match{$tag1} if (exists $match{$tag1});
    $panel{$tag2}=$match{$tag2} if (exists $match{$tag2});
  }
  close (PANEL);
  
  if ($use_mito==1) {
    #Open source BIM file and create list of MT SNPs that are in matched data:
    open (OUTMT, ">$name.mtsnps.txt") || die "Unable to open mt tag file to write to: $!\n";
    while (my ($a,$b) = each %match) {
      if ($a=~/^26_/) {
	print OUTMT "$source{$a}\n" if (exists $source{$a});
      }
    }
    close (OUTMT);
    #Create plink file with macthing MT SNPs:
    system `plink --bfile $bfile --extract $name.mtsnps.txt --make-bed --out $name.mt`;
    #Use PLINK to clump SNPs, selecting on required P value for primary SNP, and then filtering using default values
    system `plink --bfile $name.mt --clump $eqtl --out $name --clump-p1 $pvalue --clump-p2 $pvalue --clump-r2 0.8 --chr 26`;
    if (-e "$name.clumped") {
      open (CLUMP, "$name.clumped") || die "Unable to open clump file $name.clumped to read: $!\n";
      while (<CLUMP>) {
	my @array=split;
	my $collect=$flip{$array[2]};
	$panel{$collect}=0;
      }
      close (CLUMP);
    }
  }

  %match=();
  %match=%panel;
  
}

#Copy over the PLINK files so the SNP column can be edited to chr_pos_ref_alt format (this makes it easier to use models on other datasets)
system `cp $bfile.bed $name.full.bed`;
system `cp $bfile.bim $name.full.bim`;
system `cp $bfile.fam $name.full.fam`;

#Reformat SNP column, selecting only those SNPs on the autosomes or MT if included, and also only keeping SNPs that are in matched data (above):
open (BIMFULL, "$name.full.bim") || die "Unable to open copied BIM file $name.full.bim to read: $!\n";
open (BIMHOLD, ">$name.hold.bim") || die "Unable to open BIM file to write to: $!\n";
open (MARK, ">$name.markers.txt") || die "Unable to open markers file to write to: $!\n";

while (<BIMFULL>) {
  chomp;
  my @array=split;
  my $chr=$array[0]; $chr=~s/chr//; $chr=~s/MT/26/; $chr=~s/M/26/;
  print BIMHOLD "$chr\t${chr}_$array[3]_$array[5]_$array[4]\t$array[2]\t$array[3]\t$array[4]\t$array[5]\n";
  my $marker="${chr}_$array[3]_$array[5]_$array[4]";
  if ((exists $match{$marker})&&((($chr>0)&&($chr<23))||($chr==26))) {
    my $qmatch="${chr}_$array[3]_$array[4]";
    if ($eqtlp{$qmatch}<=$pvalue) {
      print MARK "$marker\n";
    }
  }
}

close (BIMFULL);
close (BIMHOLD);
close (MARK);

#Use new BIM file, and then filter for SNPs that meet criteria and are in matched file, before renaming back to original name
system `mv $name.hold.bim $name.full.bim`;
system `plink --bfile $name.full --extract $name.markers.txt --make-bed --out $name.selected`;

#Can do clumping at this point if REF selection not yet performed:

if ($use_ref=~/clump/) { #eqtl and $name.selected markers are different:

  open (EQTL, "gunzip -c $eqtl |") || die "Unable to open $name.temp.eqtl.gz to read: $!\n"; #open eqtl file
  open (EQTL_OUT, ">$name.eqtl") || die "Unable to open $name.temp.eqtl.gz to read: $!\n"; #open eqtl file

  while (<EQTL>) {
    chomp;
    if ($_=~/BETA/) {
      print EQTL_OUT "$_\n";
    }
    if (!($_=~/BETA/)) {
      my @array=split;
      print  EQTL_OUT "$array[0]\t$flip{$array[1]}\t$array[2]\t$array[3]\t$array[4]\t$array[5]\t$array[6]\t$array[7]\t$array[8]\t$array[9]\t$array[10]\t$array[11]\n";
    }
  }
  close (EQTL);
  close (EQTL_OUT);
      
  system `plink --bfile $name.selected --clump $name.eqtl --out $name --clump-p1 $pvalue --clump-p2 $pvalue --clump-r2 0.8`;
  open (CLUMP, "$name.clumped") || die "Unable to open clump file $name.clumped to read: $!\n";
  open (MARK1, ">$name.clumped.markers.txt") || die "Unable to open markers file to write to: $!\n";
  while (<CLUMP>) {
    my @array=split;
    if (!($_=~/TOTAL/)) {
      print MARK1 "$array[2]\n";
    }
  }
  close (CLUMP);

  system `plink --bfile $name.full --extract $name.clumped.markers.txt --make-bed --out $name.selected` if ($use_mito==1);
  system `plink --bfile $name.full --extract $name.clumped.markers.txt --chr 2-22 --make-bed --out $name.selected` if ($use_mito==0);
}
  
#Collect Pheno data from plink file

my %data=();

#If covars are not to be regressed out, simply store the data from the correct phenotype column
if ($use_covars == 0) {
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

#If covars are to be regressed out, collect covarites present in the files (between combined_s20 and RNA columns - these are standard output PLINK format files from the mitogwas pipeline
if ($use_covars == 1) {
  open (PHENO2, "$pheno_file") || die "Unable to open pheno file $pheno_file to read: $!\n";
  my @inds=(); my @head=(); my $exp; my %covars=(); my @gene=(); my $open=0; my %miss=();
  while (<PHENO2>) {
    my @array=split;
    if ($_=~/^FID/) {
      @head=@array;
    }
    else {
      push @inds, $array[0];
      $data{$array[0]}="NA";
      for (my $i=0;$i<@array;$i++) {
	if ($head[$i] eq $pheno) {
	  $exp.=",".$array[$i];
	  push @gene, $array[$i];
	  $miss{$array[0]}=0 if ($array[$i] eq "NA");
	}
	if ($head[$i]=~/RNA/) {
	  $open=0;
	}
	if ($open==1) {
	  my $number = "NA";
	  $number = sprintf "%.4f", $array[$i] if ($array[$i]=~/\d/);
	  $covars{$head[$i]}.=",".$number;
	  $miss{$array[0]}=0 if ($array[$i] eq "NA");
	}
	if ($head[$i]=~/combined_s20/) {
	  $open=1;
	}
      }
    }
  }
  close (PHENO2);
   
  my $covar_string;
  while (my ($a,$b) = each %covars) {
    my $dat=$b; $dat=~s/,//;
    $make="$a<-c($dat)";
    print "$make\n";
    $R->send($make);
    $covar_string.="+".$a
      
    #eval { $R->send($make) };
    #$R->send($make) if (!$@);
    #$covar_string.="+".$a if (!$@);
  }

  $covar_string=~s/\+//;
  
  $exp=~s/,//;
  $make="gene<-c($exp)";
  $R->run($make);
  
  $make="res<-resid(lm(gene~$covar_string))";
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

while (my ($a,$b) = each %data) {
  #print "$a\t$b\n";
}

#Edit Plink fam file to insert the phenotype data for the chosen feature

open (FAM, "$name.selected.fam") || die "Unable to open fam $name.selected.fam file to read: $!\n";
open (OUTFAM, ">$name.hold.fam") || die "Unable to open fam $name.hold.fam file to write to: $!\n";
open (REM, ">$name.rem.txt") || die "Unable to open fam $name.hold.fam file to write to: $!\n";

my $rem=0;
while (<FAM>) {
  my @array=split;
  my $data="NA";
  $data=$data{$array[0]} if (exists $data{$array[0]});
  print OUTFAM "$array[0] $array[1] $array[2] $array[3] $array[4] $data\n";
  if ($data=~/NA/) {
    print REM "$array[0] $array[1]\n";
    $rem++;
  }
}
close (FAM);

system `mv $name.hold.fam $name.selected.fam`;

if ($rem>0) {
  system `plink --bfile $name.selected --remove $name.rem.txt --make-bed --out $name.rem`;
  system `mv $name.rem.fam $name.selected.fam`;
  system `mv $name.rem.bed $name.selected.bed`;
  system `mv $name.rem.bim $name.selected.bim`;
}

