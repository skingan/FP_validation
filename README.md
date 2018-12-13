# FP_validation
scripts for testing FALCON-Phase results against TrioCanu assemblies

These scripts test the phase assignment of the minced A and B fasta files against two trioCanu assemblies of parents.

1. Run minimap four times for each pairwise combination of A or B minced fasta (QRY) X mom or dad triocanu assemblies (REF).
2. Assign parents to the phase blocks using minimap2parents.pl script.
3. The scoring script computes accuracy of phase assignment weighted by sequence length.
