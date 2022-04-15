#!/usr/bin/env python3
from Bio import SeqIO
import numpy as np


#################### Create argument parser ####################
import argparse as ap
parser = ap.ArgumentParser(description="""
Simple script to print out MSA alignments.
""", formatter_class=ap.RawTextHelpFormatter)

parser.add_argument("-f", "--file", dest="F", required=True,
        help="MSA file")
parser.add_argument("--color", dest="color", action="store_true",
        help="Option to print-colored outputs. Mainly useful for terminals that supports colors")
args = parser.parse_args()

def FindSeqs(f):
    """
    Function to pull sequences from the multi-FASTA
    """
    
    # Pull the sequences
    full_seqs = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
    
    # Only return what we need; saves us time and memory later on anyways
    seqs = {}
    for k in full_seqs.keys():
        seqs[k] = full_seqs[k].seq
    
    return seqs


# Set up the colors for print-outs
class bcolors:
        MIS = '\033[2;30;41m' #RED
        RESET = '\033[0m' #RESET COLOR

def PrintAligns(r,q, start, end, r_head, q_head, text):
    """
    Function to print out the string alignments wtih colors
    """

    r = np.array(list(r))
    q = np.array(list(q))
    
    # Correct for difference in length of headers
    hdiff = len(r_head) - len(q_head)

    if not text:
        print(str(start) + " " * (end - start + 5) + str(end))
        if hdiff > 0:
        print(f"{r_head}: {''.join(r)}")
        print(f"{q_head}: " + "".join([f"{bcolors.MIS}{q[n]}{bcolors.RESET}" if _ == False \
            else q[n]
            for n,_ in enumerate(q == r)]))
        print()
    else:
        print(str(start) + " " * (end - start + 5) + str(end))
        print(f"{r_head}: " + "".join([f"{r[n]}" if _ == False \
            else "*"
            for n,_ in enumerate(q == r)]))
        print(f"{q_head}: " + "".join([f"{q[n]}" if _ == False \
            else "*"
            for n,_ in enumerate(q == r)]))
    print()


#################### Run the function ####################
seqs = FindSeqs(args.F)

if args.color:
    text = False
else:
    text = True

seq_len = np.max([len(seqs[k]) for k in seqs.keys()])

r_start = 0
for r_end in np.arange(80, seq_len+81, step=80):
    r_head, q_head = [k for k in seqs.keys()] # Pull the headers

    if r_end < seq_len:
        r,q = [seq[r_start:r_end] for seq in seqs.values()]

    else:
        r,q = [seq[r_start:] for seq in seqs.values()]
        r_end = seq_len

    PrintAligns(r,q, r_start, r_end, r_head, q_head, text=text)
    r_start = r_end
