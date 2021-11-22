#!/usr/bin/python3

from tqdm import tqdm
import re

# Set up the argument parser
import argparse as ap
parser = ap.ArgumentParser(description="Script to pull the regions of interest from an MSA.")
parser.add_argument("-i", "--input", type=str, dest="F",
        help="File containing the MSA")
parser.add_argument("-o", "--output", type=str, dest="out",
        help="Output file name.fa")
parser.add_argument("-r", "--region", type=str, dest="region",
        help="The region of interest specified as start:stop (e.g. 100:250)")
parser.add_argument("--multi-regions", type=str, dest="multi_regions",
        help="The regions of interest specified as 'start1:stop1, start2:stop2, etc...'. Must be inside quotation marks")

args = parser.parse_args()

# Check that the required arguments are present
parser_errs = []
if not args.F:
    parser_errs.append("Input file")
if not args.out:
    parser_errs.append("Output file")
if not args.region and not args.multi_regions:
    parser_errs.append("Region of interest")

if len(parser_errs) > 0:
    [print(f"{e} missing!") for e in parser_errs]
    parser.error("See '--help' for options")

# Check if the arguments are formatted the right way
if not re.match("^f*a\Z", args.out.split(".")[-1]):
    print("Output file format is not correct. File extension must be fa/fasta.")
    parser.error("See --help for options")

#################### Create the necessary functions ####################
def FindSeqs(f):
    # Pull the sequences
    seqs = {}
    seq = None
    
    for line in tqdm(open(f, "r"), desc=f"Pulling sequences from {f}"):
        if line.startswith(">"):
            
            if seq:
                seqs[header] = "".join(seq)
            
            header = line.strip()
            seq = []
            
        else:
            seq.append(line.strip())

    seqs[header] = "".join(seq) # The final sequence will not have a header after it
    
    return seqs


#################### Run the function ####################

# Check if it is one region or multi-region
if args.multi_regions:
    regions = [x.strip().split(":") for x in args.multi_regions.split(",")]

else:
    regions = [args.region.split(":")]

# Pull the sequences from the MSA
seqs = FindSeqs(args.F)

# Prepare the output file
o = open(args.out, "w")

# Pull out the sub-regions of interest
for h,s in seqs.items():
    o.write(h + "\n" + "".join([s[int(r[0]):int(r[1])] for r in regions]) + "\n")

print(f"Done! Regions of interest stored in {args.out}")
