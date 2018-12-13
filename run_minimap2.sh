#!/bin/bash
#$ -S /bin/bash
#$ -N mm2
#$ -cwd
#$ -q default
#$ -pe smp 4

module load minimap2

REF=$1
QRY=$2
OUT=$3

# index ref
#minimap2 -t4 -k19 -w19 -d $REF.mmi $REF 

# contig 2 ref paf output
minimap2 -t 4 -x asm5 $REF $QRY 2> $OUT.stderr > $OUT.paf

# contig 2 ref sam output
#minimap2 -t 4 -aQ -x asm5 $REF $QRY 2> $OUT.stderr > $OUT.sam

# pb read 2 ref sam output
#minimap2 -t 4 -aQ -x map-pb $REF $QRY 2> $OUT.stderr > $OUT.sam

