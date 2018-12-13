#!/usr/bin/perl -w
        
####################################################################################################
#
#               Sarah B. Kingan
#               Pacific Biosciences
#               3 July 2018
#
#               Title: scorePhase_trioCanuParents.pl
#
#               Project: FALCON-Phase
#       
#               Input:  haplotigPlacementFile, parentAsgn, phaseAsgn, thresholdAccuracy
#
#               Output: STDOUT: 0 or 1 (pass threshold accuracy)
#                               STDERR: primary contig summ stat (length, accuracy, n blocks)
#                               
#                       
#
####################################################################################################

use strict;
use List::Util qw[min max];

my $usage = "scorePhase2.pl haplotig_placement.paf parentAsgn.txt phaseAsgn.txt thresholdAccuracy\n";
        
# placement file
my $placement_file = shift(@ARGV) or die $usage;
#000000F_022    47906   1122    47906   +       000000F 65274471        0       46794   46784   46784   60
#000000F_051    33968   0       33968   +       000000F 65274471        52713   86713   33968   33968   60
#000000F_014    156643  0       156643  +       000000F 65274471        73372   230650  156643  156643  60
#000000F_100    1918733 0       1918733 +       000000F 65274471        257247  2172913 1918733 1918733 60
                        
my $parent_file = shift(@ARGV) or die $usage;
#000032F:2074213-4680387 dam
#000266F:13490-477906    dam
#000038F_040:0-1416539   dam
#000043F:3523180-3623931 dam
#000057F_017:0-360980    sire
        
my $phase_file = shift(@ARGV) or die $usage;
#000000F 000000F_022:0-47906 000000F:0-46794 1.000000 1.2576 0.5385 6 3157
#000000F 000000F_051:0-33968 000000F:52713-86713 0.500202 0.4130 0.4516 19 3158
#000000F 000000F_014:0-156643 000000F:73372-230650 0.500081 0.9529 0.8350 4 3159

my $accuracy = shift(@ARGV) or die $usage;

# make parent hash
my %parent_hash;
my $parent;
open (PAR, $parent_file);
while (my $line = <PAR>) {
        chomp $line;
        my @line_array = split("\t", $line);
        if ($line_array[1] eq 'dam') {
                 $parent = 'mom';
        }
        elsif ($line_array[1] eq 'sire') {
                 $parent = 'dad';
        }
        elsif ($line_array[1] eq 'dad') {
                 $parent = 'dad';
        }
        elsif ($line_array[1] eq 'mom') {
                 $parent = 'mom';
        }
        $parent_hash{$line_array[0]} = $parent;
}

# make phase hash
my %phase_hash;
open (PH, $phase_file);
while (my $line = <PH>) {
        chomp $line;
        my @line_array = split(" ", $line);
        $phase_hash{$line_array[1]} = '0';
        $phase_hash{$line_array[2]} = '1';
}


# header
print STDERR join("\t", qw(primary_contig primary_contig_length unzipped_length n_phase_blocks accuracy)), "\n";
my $A_contigID;
my $B_contigID;
my $pcontig_depth;
my $primary_contig = 'first';
my $primary_contig_length = 0;
my $unzipped_length = 0;
my $n_haplotigs = 0;
my $add_length = 0;
my $total = 0;
my $Aparent;
my $Aphase;
my $Bparent;
my $Bphase;
# dam sire
my @phase0 = (0,0);
my @phase1 = (0,0);
my $phase0_length = 0;
my $phase1_length = 0;
my $dam_length = 0;
my $sire_length = 0;
my $phase_correct = 0;
my $phase_incorrect = 0;
my $tmp;
my @out = ("nan","nan");

my @accuracy_array;
my @length_array;

open (PL, $placement_file);
while (my $line = <PL>) {
        chomp$line;
        my @line_array = split("\t", $line);

# get contig IDs and atrributes from placement file and hashes
        $A_contigID = $line_array[0].":"."0"."-".$line_array[1];
        $B_contigID = $line_array[5].":".$line_array[7]."-".$line_array[8];
# set parents and phase from hashes
        $Aparent = "NA";
        $Aphase = "NA";
        $Bparent = "NA";
        $Bphase = "NA";
        if (exists $parent_hash{$A_contigID}) {
                $Aparent = $parent_hash{$A_contigID};
                if ($Aparent eq 'mom') {
                        $Bparent = 'dad';
                }
                else {
                        $Bparent = 'mom';
                }
        }
        if (exists $parent_hash{$B_contigID}) {
                $Bparent = $parent_hash{$B_contigID};
                if ($Bparent eq 'mom') {
                        $Aparent = 'dad';
                }
                else {
                        $Aparent = 'mom';
                }
        }
        if (exists $phase_hash{$A_contigID}) {
                $Aphase = $phase_hash{$A_contigID};
        }
        if (exists $phase_hash{$B_contigID}) {
                $Bphase = $phase_hash{$B_contigID};
        }
#               print $A_contigID, "\t", $Aparent, "\n";
#               print $B_contigID, "\t", $Bparent, "\n";

        if ($line_array[5] ne $primary_contig) { # new contig
                unless ($primary_contig eq 'first') { # first contig, don't print
                        # process, print, and reset
                        # process...
                        $phase0_length = $phase0[0] + $phase0[1];
                        $phase1_length = $phase1[0] + $phase1[1];
                        $dam_length = $phase0[0] + $phase1[0];
                        $sire_length = $phase0[1] + $phase1[1];                         
                        $total = $phase0_length + $phase1_length;
                        $phase_correct = "nan";
                        $phase_incorrect = "nan";
                        @out = ("nan","nan");
                        if ($total > 0) {
                                $tmp = max(@phase0) + max(@phase1);
                                $phase_correct = $tmp / $total;
                                $tmp = min(@phase0) + min(@phase1);
                                $phase_incorrect = $tmp / $total;
                                $out[0] = sprintf("%.3f", $phase_correct);
                                $out[1] = sprintf("%.3f", $phase_incorrect);
                        }                               
                        # print                 
                        #print join("\t", ($primary_contig, $primary_contig_length, $haplotig_length, $n_haplotigs, $pcontig_depth, $sex_chrom
, $dam_length, $sire_length, $phase0_length, $phase1_length)), "\t";
                        print STDERR join("\t", ($primary_contig, $primary_contig_length, $unzipped_length, $n_haplotigs, $out[0])), "\n";
                        push(@accuracy_array, $out[0]);
                        push(@length_array, $primary_contig_length);

                        # reset
                        $unzipped_length = 0;
                        $n_haplotigs = 0;
                        @phase0 = (0,0);
                        @phase1 = (0,0);
                        @out = ('nan','nan');
                        $total = 0;                     
                }
        }
        if ($Aphase eq '0') {
                if ($Aparent eq 'mom') {
                        $phase0[0] += $line_array[10];
                }
                elsif ($Aparent eq 'dad') {
                        $phase0[1] += $line_array[10];
                }       
        }
        elsif ($Aphase eq '1') {
                if ($Aparent eq 'mom') {
                        $phase1[0] += $line_array[10];
                }
                elsif ($Aparent eq 'dad') {
                        $phase1[1] += $line_array[10];
                }
        }
        $primary_contig = $line_array[5];
        $primary_contig_length = $line_array[6];
        $n_haplotigs++;
        $unzipped_length += $line_array[10];
}

# final print

$phase0_length = $phase0[0] + $phase0[1];
$phase1_length = $phase1[0] + $phase1[1];
$dam_length = $phase0[0] + $phase1[0];
$sire_length = $phase0[1] + $phase1[1];
$total = $phase0_length + $phase1_length;
$phase_correct = "nan";
$phase_incorrect = "nan";
@out = ("nan","nan");
if ($total > 0) {
        $tmp = max(@phase0) + max(@phase1);
        $phase_correct = $tmp / $total;
        $tmp = min(@phase0) + min(@phase1);
        $phase_incorrect = $tmp / $total;
        $out[0] = sprintf("%.3f", $phase_correct);
        $out[1] = sprintf("%.3f", $phase_incorrect);
}
# print                 
print STDERR join("\t", ($primary_contig, $primary_contig_length, $unzipped_length, $n_haplotigs, $out[0])), "\n";
push(@accuracy_array, $out[0]);
push(@length_array, $primary_contig_length);


# STDOUT SCORE
my $weighted_mean = weighted_mean(\@accuracy_array, \@length_array);
my $end = 1;
if ($weighted_mean < $accuracy) {
        $end = 0;
}
print $end, "\n";
print $weighted_mean, "\n";


sub weighted_mean {
        my ($value_ref, $weight_ref) = @_;
        my @value = @{$value_ref};
        my @weight = @{$weight_ref};
        my $prod = 0;
        my $num_sum = 0;
        my $denom_sum = 0;
        for (my $i=0; $i<scalar@value; $i++) {
                if ($value[$i] ne 'nan') {
                        $prod = $value[$i] * $weight[$i];
                        $num_sum += $prod;
                        $denom_sum += $weight[$i];
                }
        }
        my $weighted_mean = $num_sum / $denom_sum,;
        my $out = sprintf("%.4f", $weighted_mean);
        return $out;
}


