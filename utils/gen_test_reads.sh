#!/bin/bash
# Thanks Hanquan for kindly make this script
################################################################

# Install wgsim
git clone https://github.com/lh3/wgsim.git
cd wgsim
gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm
# Optional, add excutable to PATH
PATH=$PWD:$PATH

# Retrieve reference sequence
# In this case we use kind of E.coli
wget -O ref.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip ref.fa.gz

# Simulate sequence reads from a reference genome using wgsim
wgsim -N 1000000 -1 75 -2 75 -e 0.01 -r 0.001 ref.fa read1.fq read2.fq
# -N  number of read pairs to generate.
# -1  length of read 1.
# -2  length of read 2.
# -e  base error rate (i.e. sequencing error).
# -r  mutation rate. This instructs how the resulting data will be different from the original reference. 0.001 means there will be 1 mutation in every 1000 base.

################################################################


