#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

# ============================================================
# extract_hpa_tissue_specificity.pl
#
# Extract selected columns from Human Protein Atlas proteinatlas.tsv
#
# Required columns:
#   Gene
#   Gene description
#   RNA tissue specificity
#   RNA tissue specific nTPM
#
# Usage:
#   perl extract_hpa_tissue_specificity.pl \
#       --input proteinatlas.tsv \
#       --output hpa_tissue_specificity.tsv
# ============================================================

my ($input, $output);

GetOptions(
    "input|i=s"  => \$input,
    "output|o=s" => \$output,
) or die usage();

die usage() unless $input && $output;

open my $IN,  "<", $input  or die "[ERROR] Cannot open input file: $input\n";
open my $OUT, ">", $output or die "[ERROR] Cannot write output file: $output\n";

my $header = <$IN>;
chomp $header;
$header =~ s/\r$//;

my @cols = split /\t/, $header, -1;

# Remove double quotes from column names
for my $c (@cols) {
    $c =~ s/^"//;
    $c =~ s/"$//;
    $c =~ s/^\s+|\s+$//g;
}

my %col_index;
for my $i (0 .. $#cols) {
    $col_index{$cols[$i]} = $i;
}

my @target_cols = (
    "Gene",
    "Gene description",
    "RNA tissue specificity",
    "RNA tissue specific nTPM",
);

for my $col (@target_cols) {
    die "[ERROR] Required column not found: $col\n"
        unless exists $col_index{$col};
}

print $OUT join("\t", @target_cols), "\n";

while (my $line = <$IN>) {
    chomp $line;
    $line =~ s/\r$//;

    next if $line =~ /^\s*$/;

    my @fields = split /\t/, $line, -1;

    my @values;
    for my $col (@target_cols) {
        my $idx = $col_index{$col};
        my $value = defined $fields[$idx] ? $fields[$idx] : "";

        # Remove surrounding quotes if present
        $value =~ s/^"//;
        $value =~ s/"$//;

        push @values, $value;
    }

    print $OUT join("\t", @values), "\n";
}

close $IN;
close $OUT;

print "[INFO] Done.\n";
print "[INFO] Output: $output\n";

sub usage {
    return <<"USAGE";
Usage:
  perl extract_hpa_tissue_specificity.pl --input proteinatlas.tsv --output hpa_tissue_specificity.tsv

Options:
  --input|-i     Input HPA proteinatlas.tsv file
  --output|-o    Output TSV file

Example:
  perl extract_hpa_tissue_specificity.pl \\
      -i proteinatlas.tsv \\
      -o hpa_tissue_specificity.tsv
USAGE
}

