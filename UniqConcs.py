#!/usr/bin/env python3

import gzip
import re
import os
import multiprocessing as mp

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

from tqdm import tqdm


#################### Set up argument parser ####################

import argparse as ap
parser = ap.ArgumentParser(description="""
Script to find the uniquely concordant pairs of reads according to BowTie2
""", formatter_class=ap.RawTextHelpFormatter)

parser.add_argument("-i", "--idir", dest="idir", required=True,
        help="Directory containing the SAM files from the BowTie2 alignments")
parser.add_argument("--pattern", dest="patt", required=True,
        help='The regex pattern matching the files in double-quotes ("")' )
parser.add_argument("-t", "--threads", dest="nthreads",
        default=2, type=int,
        help="Number of threads to use")
parser.add_argument("-o", "--out", dest="out", default=None,
        help="""Specify file to save the same as (image.[png|jpg|pdf|tiff|svg]).
        If no name is given, image will not be created""")

args = parser.parse_args()


#################### Set up the function ####################
def ParseLine(line, outname):
    """Script to parse the reads
    We want to find:
                    1) Total number of reads
                    2) Number of alignments
                    3) Number of concordant pairs
                    4) Number of uniquely aligned concordant pairs with len() > 20
                    5) Store the uniq_concordant alignments in a file

"""
    
    global tot, aligned, concs, uniq_concs
    
    
    if not line.startswith("@"):
        tot += 1

        # Check cigar string; if != '*', then there is a valid alignment
        if line.split("\t")[5] != "*":
            aligned += 1

            # YT:Z:CP flag indicates concordant pair, according to BowTie2
            if "YT:Z:CP" in line:
                concs += 1

                # XS flag indidcates multiple mapping sites
                if "XS" not in line and len(line.split("\t")[9]) > 20: 
                    uniq_concs += 1

                    # Save the line to the file
                    gzip.open(outname, "at").write(line)

def CalculateStats(f):
    
    out = f.split("/")[-1].split(".sam")[0] + "_unique-concordant.sam.gz"
    sample = f.split("/")[-1].split(".paired")[0]
    gzip.open(out, "wt").write("") # Just to overwrite the file
    
    global tot, aligned, concs, uniq_concs
    tot, aligned, concs, uniq_concs = 0,0,0,0
    
    # Iterate through the lines
    [ParseLine(line, out) for line in gzip.open(f, "rt")]
    
    # Print out results
    unalign = tot - aligned
    aligned_perc = round((aligned/tot) * 100,2)
    conc_perc = round(((concs - uniq_concs)/tot) * 100,2)
    uniq_perc = round((uniq_concs/tot) * 100,2)

    print(f"{sample}: TotalAligned: {aligned:,} ({aligned_perc}%) -- Concordant: {concs:,} ({conc_perc}%) -- Uniquely Concordant: {uniq_concs:,} ({uniq_perc}%)")    
    
    # Add the results in
    open("Results.tsv", "a").write("\t".join([str(sample), str(unalign), 
                                              str(aligned - concs - uniq_concs), 
                                              str(concs - uniq_concs), str(uniq_concs)]))

def PlotRes(out=None):
    
    df = pd.read_csv("./Results.tsv", sep="\t",
                    names = ["Unaligned", "Discordant Align", 
                             "Concordant Align", "Unique Concordant Align"])
    
    
    
    df.plot(kind="bar", stacked=True,
        edgecolor="k", linewidth=1,
       figsize=(12,8))


    plt.title("BowTie2 Read Alignment Statistics\n")
    plt.ylabel("Number of reads")
    plt.xlabel("Sample name")

    plt.xticks(rotation=45, ha='right')

    plt.legend(bbox_to_anchor=(1,0,0.25,1),
              facecolor="white", edgecolor="k", frameon=True)
    
    if out:
        plt.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")


#################### Run the functions ####################

# Find the files
pattern = args.patt.replace("*", ".*") # correct pattern to work wtih PyRegex
files = [x for x in os.listdir(args.idir) if re.search(pattern, x)]

# Print out the given arguments
args_dict = args.__dict__
print("\nProvided arguments:")
[print(f"{x}: {args_dict[x]}") for x in args_dict]
print()

print(f"Files ({len(files)}): {files}")
print()


# Run the Calculation function in parallel
pool = mp.Pool(args.nthreads)
open("Results.tsv", "w").write("")
for _ in tqdm(pool.imap_unordered(CalculateStats, files), total=len(files)):
    pass

# Plot the results
if args.out:
    PlotRes(args.out)
