#!/usr/bin/perl -w

use strict;

my $input = shift || undef;

if (!$input) {
  print "Usage: $0 <pseq_lookup_output>\n";
  exit(1);
}

my %outTable;
my %COLUMNS;
my %LOCI;

my $COL_COUNT = 0;
my $LOC_COUNT = 0;

open INPUT, $input or die $!;

my $line;
while ($line = <INPUT>) {
  chomp($line);
  next if ($line =~ m/^\#/);

  my @fields = split(/\s+/, $line);
  die("Invalid line: $line\n") if (scalar(@fields) < 2);

  my $loc = $fields[0];
  my $col = $fields[1];

  my $val = '';
  for (my $i = 2; $i < scalar(@fields); $i++) {
    $val .= " " if ($i > 2);
    $val .= $fields[$i];
  }
  $val = '.' if ($val eq '');

  $COLUMNS{$col} = $COL_COUNT++ if (!defined($COLUMNS{$col}));
  $LOCI{$loc} = $LOC_COUNT++ if (!defined($LOCI{$loc}));

  if (!exists($outTable{$loc}{$col})) {
    $outTable{$loc}{$col} = $val;
  }
  else {
    $outTable{$loc}{$col} .= ",".$val;
  }
}

close(INPUT);

my @COLUMNS = sort byColumn (keys(%COLUMNS));


print "LOCUS";
foreach my $col (@COLUMNS) {
  print "\t$col";
}
print "\n";

my @LOCI = sort byLocus keys(%outTable);
foreach my $loc (@LOCI) {
  print "$loc";

  foreach my $col (@COLUMNS) {
    my $val = '.';
    $val = $outTable{$loc}{$col} if (defined($outTable{$loc}{$col}));

    print "\t$val";
  }

  print "\n";
}


sub byColumn {
  return ($COLUMNS{$a} <=> $COLUMNS{$b});
}

sub byLocus {
  return ($LOCI{$a} <=> $LOCI{$b});
}
