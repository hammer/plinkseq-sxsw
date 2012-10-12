#!/usr/bin/perl -w

use strict;

my $input = shift || undef;

if (!$input) {
  print "Usage: $0 <pseq_lookup_output>\n";
  exit(1);
}

my %outTable;
my %lastRefAlt;

my %COLUMNS;
my %LOCI;
my %TRUE_LOCI;

my $COL_COUNT = 0;
my $LOC_COUNT = 0;

open INPUT, $input or die $!;

my $line;
while ($line = <INPUT>) {
  chomp($line);
  if ($line =~ m/^\#/) {
    if ($line =~ m/^\##([^,]+),.*/) {
      my $col = $1;
      $COLUMNS{$col} = $COL_COUNT++ if (!defined($COLUMNS{$col}));
    }
    next;
  }

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

  # Reset the new alleles:
  if ($col eq "allele_ref") {
	delete $lastRefAlt{$loc};
  }

  my $useLoc = $loc;
  $useLoc = $lastRefAlt{$loc} if (exists($lastRefAlt{$loc}));

  if (!exists($outTable{$useLoc}{$col})) {
        $outTable{$useLoc}{$col} = $val;
  }
  else {
        $outTable{$useLoc}{$col} .= ",".$val;
  }

  # Transfer the data to ref-alt format, if applicable:
  if (exists($outTable{$loc}) && exists($outTable{$loc}{"allele_ref"}) && exists($outTable{$loc}{"allele_alt"}) && ($col eq "allele_ref" || $col eq "allele_alt")) {
	my $locRefAlt = $loc."_".$outTable{$loc}{"allele_ref"}."_".$outTable{$loc}{"allele_alt"};

	$outTable{$locRefAlt} = $outTable{$loc};
	delete $outTable{$loc};

	$LOCI{$locRefAlt} = $LOCI{$loc};
	delete $LOCI{$loc};

	$TRUE_LOCI{$locRefAlt} = $loc;

	$lastRefAlt{$loc} = $locRefAlt;
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
  my $printLoc = $loc;
  $printLoc = $TRUE_LOCI{$loc} if (exists($TRUE_LOCI{$loc}));
  print "$printLoc";

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
