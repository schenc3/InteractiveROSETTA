#!/usr/bin/perl -w

$num_args = $#ARGV + 1;
if ($num_args != 3) {
    print "\nUsage: gen-sequence-profile.pl <input.fasta> <input_checkpoint> <output>\n";
    print "         <input.fasta>: The FASTA file for the protein's whose profile will be generated\n";
    print "         <input_checkpoint>: This is the output of blastpgp, which is a legacy BLAST executable\n";
    print "                             Ex. blastpgp -i input.fasta -j 2 -d nr.15 -C checkpoint.pssm\n";
    print "                             The file checkpoint.pssm corresponds to <input_checkpoint>\n";
    print "         <output>: The name of the output file\n";
    exit;
}
 
my $checkpoint=$ARGV[1];
my $fasta=$ARGV[0];
my $output=$ARGV[2];

my $sequence = read_fasta( $fasta );
print "Sequence: $sequence";
# parse & fortran-ify the checkpoint matrix.
my @checkpoint_matrix;
@checkpoint_matrix = &parse_checkpoint_file("$checkpoint");
@checkpoint_matrix =
  &finish_checkpoint_matrix( $sequence, @checkpoint_matrix );
&write_checkpoint_file( "$output", $sequence,
    @checkpoint_matrix );

sub read_fasta {
    my $fn = shift;

    open( SEQFILE, $fn ) or die "Error opening fasta file $fn!\n";
    my $has_comment = 0;
    my $eof;
    my $seq = '';
    while ( my $line = <SEQFILE> ) {
        $eof = $line;
        $line =~ s/\s//g;
        if ( $line =~ /^>/ ) {
            $has_comment = 1;
        }
        else {
            chomp $line;
            $seq .= $line;
        }
    }

    my $has_eof = ( $eof =~ /\n$/ );
    ( $has_comment && $has_eof )
      or die
      "fasta file (given $fn) must have a comment and end with a new line!\n";
    close SEQFILE or die $!;

    return $seq;
}

## parse_checkpoint_file -- parses a PSI-BLAST binary checkpoint file.
#
# args:  filename of checkpoint file.
# rets:  N x 20 array containing checkpoint weight values, where N
#        is the size of the protein that BLAST thought it saw...

sub parse_checkpoint_file {
    my $filename = shift;
    my $buf;
    my $seqlen;
    my $seqstr;
    my $i;
    my $j;
    my @aa_order = split( //, 'ACDEFGHIKLMNPQRSTVWY' );
    my @altschul_mapping =
      ( 0, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, 18 );
    my @w;
    my @output;

    open( INPUT, $filename ) or die("Couldn't open $filename for reading.\n");

    read( INPUT, $buf, 4 ) or die("Couldn't read $filename!\n");
    $seqlen = unpack( "i", $buf );

    print "Sequence length: $seqlen";

    read( INPUT, $buf, $seqlen ) or die("Premature end: $filename.\n");
    $seqstr = unpack( "a$seqlen", $buf );

    for ( $i = 0 ; $i < $seqlen ; ++$i ) {
        read( INPUT, $buf, 160 ) or die("Premature end: $filename, line: $i\n");
        @w = unpack( "d20", $buf );

        for ( $j = 0 ; $j < 20 ; ++$j ) {
            $output[$i][$j] = $w[ $altschul_mapping[$j] ];
        }
    }

    return @output;
}

## finish_checkpoint_matrix -- "finishes" the parsed PSI-BLAST checkpoint matrix,
##                             by adding pseudo-counts to any empty columns.
#
# args:  1) sequence string  2) array returned by parse_checkpoint_file
# rets:  "finished" array.  suitable for printing, framing, etc.

sub finish_checkpoint_matrix {
    my ( $s, @matrix ) = @_;
    my @sequence = split( //, $s );
    my $i;
    my $j;
    my $sum;
    my @words;
    my @b62;
    my @blos_aa =
      ( 0, 14, 11, 2, 1, 13, 3, 5, 6, 7, 9, 8, 10, 4, 12, 15, 16, 18, 19, 17 );

    my %aaNum = (
        A => 0,
        C => 1,
        D => 2,
        E => 3,
        F => 4,
        G => 5,
        H => 6,
        I => 7,
        K => 8,
        L => 9,
        M => 10,
        N => 11,
        P => 12,
        Q => 13,
        R => 14,
        S => 15,
        T => 16,
        V => 17,
        W => 18,
        Y => 19,
        X => 0
    );

    ( length($s) == scalar(@matrix) )
      or die("Length mismatch between sequence and checkpoint file!\n");

    my $blosum62 = <<BLOSUM;
#  BLOSUM Clustered Target Frequencies=qij
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
   A      R      N      D      C      Q      E      G      H      I      L      K      M      F      P      S      T      W      Y      V
0.0215
0.0023 0.0178
0.0019 0.0020 0.0141
0.0022 0.0016 0.0037 0.0213
0.0016 0.0004 0.0004 0.0004 0.0119
0.0019 0.0025 0.0015 0.0016 0.0003 0.0073
0.0030 0.0027 0.0022 0.0049 0.0004 0.0035 0.0161
0.0058 0.0017 0.0029 0.0025 0.0008 0.0014 0.0019 0.0378
0.0011 0.0012 0.0014 0.0010 0.0002 0.0010 0.0014 0.0010 0.0093
0.0032 0.0012 0.0010 0.0012 0.0011 0.0009 0.0012 0.0014 0.0006 0.0184
0.0044 0.0024 0.0014 0.0015 0.0016 0.0016 0.0020 0.0021 0.0010 0.0114 0.0371
0.0033 0.0062 0.0024 0.0024 0.0005 0.0031 0.0041 0.0025 0.0012 0.0016 0.0025 0.0161
0.0013 0.0008 0.0005 0.0005 0.0004 0.0007 0.0007 0.0007 0.0004 0.0025 0.0049 0.0009 0.0040
0.0016 0.0009 0.0008 0.0008 0.0005 0.0005 0.0009 0.0012 0.0008 0.0030 0.0054 0.0009 0.0012 0.0183
0.0022 0.0010 0.0009 0.0012 0.0004 0.0008 0.0014 0.0014 0.0005 0.0010 0.0014 0.0016 0.0004 0.0005 0.0191
0.0063 0.0023 0.0031 0.0028 0.0010 0.0019 0.0030 0.0038 0.0011 0.0017 0.0024 0.0031 0.0009 0.0012 0.0017 0.0126
0.0037 0.0018 0.0022 0.0019 0.0009 0.0014 0.0020 0.0022 0.0007 0.0027 0.0033 0.0023 0.0010 0.0012 0.0014 0.0047 0.0125
0.0004 0.0003 0.0002 0.0002 0.0001 0.0002 0.0003 0.0004 0.0002 0.0004 0.0007 0.0003 0.0002 0.0008 0.0001 0.0003 0.0003 0.0065
0.0013 0.0009 0.0007 0.0006 0.0003 0.0007 0.0009 0.0008 0.0015 0.0014 0.0022 0.0010 0.0006 0.0042 0.0005 0.0010 0.0009 0.0009 0.0102
0.0051 0.0016 0.0012 0.0013 0.0014 0.0012 0.0017 0.0018 0.0006 0.0120 0.0095 0.0019 0.0023 0.0026 0.0012 0.0024 0.0036 0.0004 0.0015 0.0196
BLOSUM

    $i = 0;
    my @blosum62 = split( /\n/, $blosum62 );

    # read/build the blosum matrix
    foreach my $blosumline (@blosum62) {
        next if ( $blosumline !~ /^\d/ );
        chomp $blosumline;
        @words = split( /\s/, $blosumline );

        for ( $j = 0 ; $j <= $#words ; ++$j ) {
            $b62[ $blos_aa[$i] ][ $blos_aa[$j] ] = $words[$j];
            $b62[ $blos_aa[$j] ][ $blos_aa[$i] ] = $words[$j];
        }

        ++$i;
    }

    # normalize the blosum matrix so that each row sums to 1.0
    for ( $i = 0 ; $i < 20 ; ++$i ) {
        $sum = 0.0;

        for ( $j = 0 ; $j < 20 ; ++$j ) {
            $sum += $b62[$i][$j];
        }

        for ( $j = 0 ; $j < 20 ; ++$j ) {
            $b62[$i][$j] = ( $b62[$i][$j] / $sum );
        }
    }

    # substitute appropriate blosum row for 0 rows
    for ( $i = 0 ; $i <= $#matrix ; ++$i ) {
        $sum = 0;

        for ( $j = 0 ; $j < 20 ; ++$j ) {
            $sum += $matrix[$i][$j];
        }

        if ( $sum == 0 ) {
            for ( $j = 0 ; $j < 20 ; ++$j ) {
                $matrix[$i][$j] = $b62[ $aaNum{ $sequence[$i] } ][$j];
            }
        }
    }

    return @matrix;
}

sub write_checkpoint_file {
    my ( $filename, $sequence, @matrix ) = @_;
    my $row;
    my $col;
    my @seq = split( //, $sequence );

    open( OUTPUT, ">$filename" );

    die("Length mismatch between sequence and checkpoint matrix!\n")
      unless ( length($sequence) == scalar(@matrix) );

    print OUTPUT scalar(@matrix), "\n";

    for ( $row = 0 ; $row <= $#matrix ; ++$row ) {
        print OUTPUT "$seq[$row] ";
        for ( $col = 0 ; $col < 20 ; ++$col ) {
            printf OUTPUT "%6.4f ", $matrix[$row][$col];
        }
        print OUTPUT "\n";
    }

    print OUTPUT "END";

    close(OUTPUT);
}
