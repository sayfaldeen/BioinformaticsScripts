# Bioinformatics scripts

- This is an informal repository of scripts I have created for bioinformatics analysis
- It is definitely a work in progress and mainly an effort at organizing my larger and/or more useful scripts

<br>

## Contents

### CheckSampleComps.sh
- Runs `kraken2` and analyzes which organisms the reads in the sample may be from
- `KrakenPiePlots.py` is used within the RunK2.sh script

<br>

### pull-msa-regions.py
- Exctracts region(s) of interest from an MSA and outputs the regions into a fasta file
- Example usage: `./pull-msa-regions.py -i input-msa.fa -o output.fa -r 100:250`
  - Output the regions between nucleotides 100 and 250
- Example usage: `./pull-msa-regions.py -i input-msa.fa -o output.fa --multi-region "100:250, 400:600, 7000:9500"`
  - Output the alignments within the each of the regions into an fasta file

### TajimasD.py
- Script to calculate Tajima's D from a given MSA file
- Example usage: `./TajimasD.py -f MSA.aln.fa`
