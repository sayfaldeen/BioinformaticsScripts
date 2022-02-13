#!/usr/bin/env python3

# Import packages
import numpy as np
from itertools import combinations
import time
from tqdm import tqdm

# Set up argument parser
import argparse as ap
parser = ap.ArgumentParser(description="""
Script to calculate Tajima's D given an MSA file.""")
parser.add_argument("-f", "--file", dest="F",
        help="File containing the MSA")
parser.add_argument("--print-intermediates", dest="print_all", action="store_true",
        help="Option to print out intermediate values (e.g. a1 ... e2)")
args = parser.parse_args()

# Make sure that necessary arguments are present
if not args.F:
    parser.error("File missing! See --help for options")

#################### Set up functions ####################
def PullSeqs(f):
    """
    Pull sequences from a MSA file
    
    Parameters
    ----------
    f: MSA file
    
    Returns
    -------
    Dictionary with {header: sequence}
    """
    
    # Pull the sequences
    seqs = {}
    seq = None
    
    for line in open(f, "r"):
        if line.startswith(">"):
            
            if seq:
                seqs[header] = "".join(seq)
            
            header = line.strip()
            seq = []
            
        else:
            seq.append(line.strip().replace(" ", ""))

    seqs[header] = "".join(seq) # The final sequence will not have a header after it to activate the loop
    
    # print(f"{len(seqs)} sequences found")
    return seqs


def CalcPi(s1, s2):
    return np.sum([1 for n1,n2 in zip(s1, s2) if n1 != n2])

def FindSegSites(s1, s2):
    global S
    return [n for n,nucs in enumerate(zip(s1, s2)) \
        if nucs[0] != nucs[1]]

def TD(f, print_all=False):
    """
    Function to calculate Tajima's D
    
    Parameters
    ----------
    f: File containing MSA
    print_all: Boolean indicating whether to print out intermediate values (a1,a2 .. e1,e2)
    
    
    Returns
    -------
    
    """
    
    t0 = time.time()
    seqs = PullSeqs(f)
    print(f"Time to pull seqs: {round(time.time() - t0, 4)}s")
    
    lseqs = len([k for k in seqs.keys()])
    ns = len(seqs[[x for x in seqs.keys()][0]]) # Number of sites (nucleotides)
    
    t0 = time.time()
    s = [] # List for all segregation sites
    pis = []
    for s1,s2 in tqdm(combinations(seqs.keys(), 2), desc="Calculating Tajima's D"):
        pis.append(np.mean(CalcPi(seqs[s1],seqs[s2])))
        s.extend(FindSegSites(seqs[s1],seqs[s2]))
        
    pi = np.mean(pis)
    del pis
    
    S = len(set(s))
    del s
    print(f"Time to find segregation sites and pi: {round(time.time() - t0, 4)}s")
    
    
#     # Calculate intermediate values for Tajima's D
    t0 = time.time()
    a1 = np.sum([1/i for i in range(1, lseqs)])
    a2 = np.sum([1/i**2 for i in range(1, lseqs)])
    
    b1 = (lseqs + 1) / (3* (lseqs - 1))
    b2 = (2*(lseqs**2 + lseqs + 3)) / (9*lseqs * (lseqs-1))
    
    c1 = b1 - (1/a1)
    c2 = b2 - ((lseqs+2) / (a1 * lseqs)) + (a2 / a1**2)
    
    e1 = c1/a1
    e2 = c2/(a1**2 + a2)
    
    # Calculate Tajima's D
    D = (pi - (S/a1)) / np.sqrt((e1*S) + ((e2 * S) * (S-1)))
    
    print(f"Time to calculate intermediate values: {round(time.time() - t0, 4)}s")
    
    # Print the results
    print(f"""
Number of sequences: {lseqs}
Total number of sites: {ns}

pi: {round(pi,4)}
Segregating sites: {S}
Tajima's D: {round(D,4)}
""")
    
    if args.print_all:
        print(f"""
    a1: {round(a1,4)}
    a2: {round(a2,4)}
    b1: {round(b1,4)}
    b2: {round(b2,4)}
    c1: {round(c1,4)}
    c2: {round(c2,4)}
    e1: {round(e1,4)}
    e2: {round(e2,4)}
    """)

TD(f=args.F)
