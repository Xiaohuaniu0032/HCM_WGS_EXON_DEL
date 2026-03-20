use strict;
use warnings;
use File::Basename;


# FHOD3 GENE BAM. NOT FULL BAM FILE
my ($child_bam,$father_bam,$mother_bam,$outdir) = @ARGV;

# software
my $whatshap = "/home/fulongfei/miniconda3/bin/whatshap";
my $freebayes = "/ehpcdata/fulongfei/software/freebayes";
my $bgzip = "/home/fulongfei/miniconda3/bin/bgzip";
my $tabix = "/home/fulongfei/miniconda3/bin/tabix";
my $samtools = "/opt/PUB/software/Bin/samtools-1.16.1/bin/samtools";

# database file
my $ref = "/ehpcdata/fulongfei/database/ref/hg19/hg19.fa";



my $child_name = (split /\./, basename($child_bam))[0]; # FHOD3_25B09089386.bam => FHOD3_25B09089386
$child_name =~ s/^FHOD3_//;

my $father_name = (split /\./, basename($father_bam))[0];
$father_name =~ s/^FHOD3_//;

my $mother_name = (split /\./, basename($mother_bam))[0];
$mother_name =~ s/^FHOD3_//;


# make outdir
if (!-d $outdir){
	`mkdir -p $outdir`;
}

# creat a run shell
my $name = (split /\./, basename($child_bam))[0]; # FHOD3_25B09089387.bam => FHOD3_25B09089387
my $runsh = "$outdir/$name\.phasing.sh";
open SH, ">$runsh" or die;


# use freebayes
my $gene_vcf = "$outdir/$name\.freebayes.trio.vcf";
my $cmd = "$freebayes --bam $child_bam --bam $father_bam --bam $mother_bam --fasta-reference $ref --min-alternate-fraction 0.1 --min-alternate-count 3 --min-coverage 10 --ploidy 2 --vcf $gene_vcf";
print SH "$cmd\n\n";


# phasing
my $phased_vcf = "$outdir/$name\.freebayes.trio.phased.vcf";
my $pedigree_ped = "$outdir/$name\_pedigree.ped";
# https://whatshap.readthedocs.io/en/latest/guide.html#phasing-pedigrees
# Fields: family, individual_id, paternal_id, maternal_id, sex, phenotype
# FAMILY01 the_child father mother 0 1

# WhatsHap checks the column count, but uses only
# column 2: individual ID
# column 3: paternal ID
# column 4: maternal ID

open PED, ">$pedigree_ped" or die;
print PED "FAMILY01\t$child_name\t$father_name\t$mother_name\t0\t1\n";
close PED;

my $reads_used_for_phase = "$outdir/$name\.reads.used.for.phase.txt";
$cmd = "$whatshap phase -o $phased_vcf --ped $pedigree_ped --reference $ref --output-read-list $reads_used_for_phase $gene_vcf $child_bam $father_bam $mother_bam";
print SH "$cmd\n\n";

$cmd = "$bgzip $phased_vcf";
print SH "$cmd\n";

$cmd = "$tabix $phased_vcf\.gz";
print SH "$cmd\n\n";


# haplo-tag
my $haplotag_list = "$outdir/$name\.HAPLOTAG_LIST";
my $haplotag_out = "$outdir/$name\.haplotag.out";
$cmd = "$whatshap haplotag --reference $ref --output-haplotag-list $haplotag_list -o $haplotag_out $phased_vcf\.gz $child_bam";
print SH "$cmd\n\n";

# split bam into hap.bam
my $hap1 = "$outdir/$name\.hap1.bam";
my $hap2 = "$outdir/$name\.hap2.bam";

$cmd = "$whatshap split --output-h1 $hap1 --output-h2 $hap2 $child_bam $haplotag_list";
print SH "$cmd\n\n";

$cmd = "$samtools index $hap1";
print SH "$cmd\n";
$cmd = "$samtools index $hap2";
print SH "$cmd\n\n";


# call variants for each hap
# my $hap1_gvcf = "$outdir/$name\.hap1.gvcf";
# USE .vcf NOT .gvcf

my $hap1_vcf = "$outdir/$name\.hap1.vcf";
$cmd = "$freebayes --bam $hap1 --fasta-reference $ref --min-alternate-fraction 0.7 --min-alternate-count 3 --min-coverage 5 --ploidy 2 \>$hap1_vcf";
print SH "$cmd\n\n";

my $hap2_vcf = "$outdir/$name\.hap2.vcf";
$cmd = "$freebayes --bam $hap2 --fasta-reference $ref --min-alternate-fraction 0.7 --min-alternate-count 3 --min-coverage 5 --ploidy 2 \>$hap2_vcf";
print SH "$cmd\n";


close SH;
`chmod 755 $runsh`;





# select N SNV to construct HAPLO SEQUENCE


# COMPARE EACH SNV POS AMONG HAPLO





# generate consensus

# clustal align with ref

# minimap2 align hap to ref

# sam2tsv

# process tsv

# extract CDS

# blast CDS against IMGT CDS.fa



