#!/usr/bin/perl -w
                
####################################################################################################
#       
#               Sarah B. Kingan
#               Pacific Biosciences
#               7 February 2018
#
#               Title: minimap2parent.pl
#
#               Project: phase unzip
#       
#               Input:  4 paf files from minimap2 (A_dam, A_sire, B_dam, B_sire) AB_pair.txt
#               Output: assignment of haplotig to dam or sire
#                       
####################################################################################################

use strict;

my $usage = "minimap2parents.pl A_dam.paf A_sire.paf B_dam.paf B_sire.paf AB_pair.txt\n";


# file input
my $A_dam_file = shift(@ARGV) or die $usage;
my $A_sire_file = shift(@ARGV) or die $usage;
my $B_dam_file = shift(@ARGV) or die $usage;
my $B_sire_file = shift(@ARGV) or die $usage;
my $pairs = shift(@ARGV) or die $usage;
        
# make hashes   
my $A_dam_hash_ref = paf2hash($A_dam_file);
my %A_dam_hash = %$A_dam_hash_ref;
my $A_sire_hash_ref = paf2hash($A_sire_file);
my %A_sire_hash = %$A_sire_hash_ref;
my $B_dam_hash_ref = paf2hash($B_dam_file);
my %B_dam_hash = %$B_dam_hash_ref;
my $B_sire_hash_ref = paf2hash($B_sire_file);
my %B_sire_hash = %$B_sire_hash_ref;    
                                        
# print header                  
#print "dam sire\n";    
                        
# process pairs                 
my @A;                          
my @B;                          
open (PAIRS, $pairs);   
while (my $line = <PAIRS>) {
        chomp $line;
        my @line_array = split("\t", $line);
        if ($line_array[0] =~ /[0-9]{6}F/) {
                push(@A, $line_array[0]);
        }
        if ($line_array[1] =~ /[0-9]{6}F/) {
                push(@B, $line_array[1]);
        }
}

# check pairs
if (scalar@A != scalar@B) {
        print "ERROR A and B lists are different lengths\n";
        die;
}

# assign parents
my $dam;
my $sire;
my @PID_array = (0,0,0,0); # Adam Asire Bdam Bsire
my %map_asgn_hash;
# key = contigID
# value = dam|sire
for (my $i = 0; $i<scalar(@A); $i++) {
        $dam = 'none';
        $sire = 'none';
# update PIDs with hash
        @PID_array = (0,0,0,0);
        if (exists $A_dam_hash{$A[$i]}{'PID'}) {
                $PID_array[0] = $A_dam_hash{$A[$i]}{'PID'};
        }
        if (exists $A_sire_hash{$A[$i]}{'PID'}) {
                $PID_array[1] = $A_sire_hash{$A[$i]}{'PID'};
        }
        if (exists $B_dam_hash{$B[$i]}{'PID'}) {
                $PID_array[2] = $B_dam_hash{$B[$i]}{'PID'};
        }
        if (exists $B_sire_hash{$B[$i]}{'PID'}) {
                $PID_array[3] = $B_sire_hash{$B[$i]}{'PID'};
        }
# assign parents
        if ($PID_array[0] > $PID_array[1]) { # dam > sire
                $dam = $A[$i]; # initial assignment
        }
        else {
                $sire = $A[$i]; # initial assignment
        }
        if ($PID_array[2] > $PID_array[3]) { # dam > sire
                if ($dam ne 'none') { # don't assign, conflict!
                        $dam = 'none'; # revert back IS THIS RIGHT???
                        next;
                }
                else {
                        $dam = $B[$i]; # assign
                }
        }
        else {
                if ($sire ne 'none') {  # don't assign, conflict!
                        $sire = 'none'; # revert back IS THIS RIGHT???
                        next;
                }
                else {
                        $sire = $B[$i]; # assign
                }
        }
        $map_asgn_hash{$dam} = 'dam';
        $map_asgn_hash{$sire} = 'sire';
}

# print mapping assignments
for my $key (keys %map_asgn_hash) {
        print $key, "\t", $map_asgn_hash{$key}, "\n";
}


sub getPrimary {
        my ($contig) = @_; 
        my $base = 'error';
        if ( $contig =~ /([0-9]{6}F)/ ) {
                $base = $1;
        }
        return $base;
}

sub paf2hash {
# key = contig ID
# value = 'alnL' =>,
#                       'PID' = >

# qryID
# qryL
# qryStart
# qryEnd
# ori
# refID
# refL
# refStart
# refEnd
# NumMatches
# AlnLength
# MapQ
        my ($paf_file) = @_;
        my %hash;
        my $PID;
        open (PAF, $paf_file);
        while (my $line = <PAF>) {
                my @line_array = split("\t", $line);
                # only mapQ 60
                if ($line_array[11] == 60) {
                        if (exists $hash{$line_array[0]}{'alnL'}) {

                                # update with longest alignment
                                if ($hash{$line_array[0]}{'alnL'} < $line_array[10]) {
                                        $hash{$line_array[0]}{'alnL'} = $line_array[10];
                                        $PID = $line_array[9] / $line_array[10];
                                        $hash{$line_array[0]}{'PID'} = $PID;
                                }
                        }
                        else {
                                $hash{$line_array[0]}{'alnL'} = $line_array[10];
                                $PID = $line_array[9] / $line_array[10];
                                $hash{$line_array[0]}{'PID'} = $PID;
                        }
                }
        }
        return \%hash;
}

exit;

