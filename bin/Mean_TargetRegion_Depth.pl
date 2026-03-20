use strict;
use warnings;
use Getopt::Long qw(GetOptions);

# ---------------------------
# Hard-coded gene coordinates
# ---------------------------
my $GENE_CHR   = 'chr18';
my $GENE_START = 33877676;
my $GENE_END   = 34360183;

# ---------------------------
# CLI args
# ---------------------------
my ($depth_file, $del_region);
my $win = 1000;

GetOptions(
  'depth|i=s' => \$depth_file,
  'del=s'     => \$del_region,
  'win=i'     => \$win,
) or die usage();

die usage() unless $depth_file && $del_region;
die "ERROR: --win must be > 0\n" unless $win > 0;

# parse del region like chr18:34232240-34241309
my ($DEL_CHR, $DEL_START, $DEL_END) = parse_region($del_region);
die "ERROR: --del region must be on $GENE_CHR\n" unless $DEL_CHR eq $GENE_CHR;

# clip deletion region to gene bounds (optional safety)
if ($DEL_START < $GENE_START) { $DEL_START = $GENE_START; }
if ($DEL_END   > $GENE_END)   { $DEL_END   = $GENE_END; }
die "ERROR: --del region is outside gene region after clipping\n" if $DEL_START > $DEL_END;

# ---------------------------
# Pre-compute windows
# ---------------------------
my $gene_len = $GENE_END - $GENE_START + 1;

my $nwin = int( ($gene_len + $win - 1) / $win ); # ceil
my @w_sum  = (0) x $nwin;  # sum depth in window
my @w_cnt  = (0) x $nwin;  # number of bases observed in window (from depth file)
# We'll output window genomic coordinates regardless of missing sites.

# ---------------------------
# Accumulators
# ---------------------------
my ($gene_sum, $gene_cnt) = (0, 0);
my ($del_sum,  $del_cnt ) = (0, 0);

open my $fh, '<', $depth_file or die "ERROR: cannot open $depth_file: $!\n";
while (my $line = <$fh>) {
  chomp $line;
  next if $line =~ /^\s*$/;
  my ($chr, $pos, $dep) = split /\t/, $line;
  next unless defined $dep;

  next unless $chr eq $GENE_CHR;
  next if $pos < $GENE_START || $pos > $GENE_END;

  # gene-level
  $gene_sum += $dep;
  $gene_cnt++;

  # del region
  if ($pos >= $DEL_START && $pos <= $DEL_END) {
    $del_sum += $dep;
    $del_cnt++;
  }

  # window index
  my $idx = int( ($pos - $GENE_START) / $win );
  # safety
  next if $idx < 0 || $idx >= $nwin;
  $w_sum[$idx] += $dep;
  $w_cnt[$idx] += 1;
}
close $fh;

# ---------------------------
# Metrics
# ---------------------------
my $gene_mean = $gene_cnt ? ($gene_sum / $gene_cnt) : 0;
my $del_mean  = $del_cnt  ? ($del_sum  / $del_cnt ) : 0;

my $ratio_del_gene = ($gene_mean > 0) ? ($del_mean / $gene_mean) : 0;

# ---------------------------
# Output summary
# ---------------------------
print "# GeneRegion\t$GENE_CHR:$GENE_START-$GENE_END\n";
print "# GeneLength(bp)\t$gene_len\n";
print "# DepthFile\t$depth_file\n";
print "# NOTE: base counts are from depth file lines (missing zero-depth sites unless samtools depth -a is used)\n";
print "\n";

print "=== SUMMARY ===\n";
print join("\t", qw(Item Region BasesObserved MeanDepth RatioToGeneMean)), "\n";

print join("\t",
  "Gene",
  "$GENE_CHR:$GENE_START-$GENE_END",
  $gene_cnt,
  fmt($gene_mean),
  fmt(1.0)
), "\n";

print join("\t",
  "Del",
  "$DEL_CHR:$DEL_START-$DEL_END",
  $del_cnt,
  fmt($del_mean),
  fmt($ratio_del_gene)
), "\n";

print "\n=== WINDOWS ===\n";
print join("\t", qw(WinIndex Region WinSize BasesObserved MeanDepth RatioToGeneMean)), "\n";

for (my $i = 0; $i < $nwin; $i++) {
  my $w_start = $GENE_START + $i * $win;
  my $w_end   = $w_start + $win - 1;
  $w_end = $GENE_END if $w_end > $GENE_END;

  my $w_size = $w_end - $w_start + 1;
  my $mean   = $w_cnt[$i] ? ($w_sum[$i] / $w_cnt[$i]) : 0;
  my $ratio  = ($gene_mean > 0) ? ($mean / $gene_mean) : 0;

  print join("\t",
    $i+1,
    "$GENE_CHR:$w_start-$w_end",
    $w_size,
    $w_cnt[$i],
    fmt($mean),
    fmt($ratio)
  ), "\n";
}

exit 0;

# ---------------------------
# Helpers
# ---------------------------
sub parse_region {
  my ($s) = @_;
  $s =~ s/\s+//g;

  # chr: start-end
  if ($s =~ /^([^:]+):(\d+)-(\d+)$/) {
    my ($c, $a, $b) = ($1, $2, $3);
    die "ERROR: region start > end in $s\n" if $a > $b;
    return ($c, $a, $b);
  }
  die "ERROR: cannot parse region '$s'. Expected format: chr18:33877676-34360183\n";
}

sub fmt {
  my ($x) = @_;
  return sprintf("%.4f", $x);
}

sub usage {
  return <<"USAGE";
Usage:
  perl fhod3_depth_stat.pl --depth FHOD3_25B09089387.depth --del chr18:34232240-34241309 [--win 1000]

Required:
  --depth/-i   samtools depth output (3 columns: chr  pos  depth)
  --del        heterozygous deletion region, e.g. chr18:34232240-34241309

Optional:
  --win        window size in bp (default 1000)

Notes:
  If your depth file is from 'samtools depth' without -a, zero-depth sites are missing and means are over observed sites only.
  For full-length base counting including zeros, generate depth with:
    samtools depth -a -r chr18:33877676-34360183 FHOD3_25B09089387.bam > FHOD3_25B09089387.depth
USAGE
}

