#!/usr/bin/env python3

import numpy as np
import pandas as pd
import random


#################### Set up argument parser ####################
import argparse as ap
parser = ap.ArgumentParser(description="""
Script to generate synthetic metagenomic communities for bacterial species.
""", formatter_class=ap.RawTextHelpFormatter)

parser.add_argument("--n-genomes", dest="nbacs",
        required=False, type=int, default=50,
        help="The number of bacteria to simulate in the community (default: 50)")
parser.add_argument("--n-samples", dest="nsamples",
        required=False, type=int, default=100,
        help="The number of samples to simulate in the community (default: 100)")
parser.add_argument("--read-count", dest="read_count",
        required=False, type=int, default=650_000,
        help="Read counts to simulate in the community for each sample (default: 650,000)")
parser.add_argument("--variable", dest="vary", action="store_true",
        help="Option to have randomly distributed values for the abundances")

args = parser.parse_args()


########################## Take in arguments ###########################
summ_url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"
nbacs = args.nbacs
nsamples = args.nsamples
read_count = args.read_count
vary = args.vary


#################### Download and import the assembly_summary file ####################
acols = ['assembly_accession', 'bioproject', 'biosample', 'wgs_master',
       'refseq_category', 'taxid', 'species_taxid', 'organism_name',
       'infraspecific_name', 'isolate', 'version_status', 'assembly_level',
       'release_type', 'genome_rep', 'seq_rel_date', 'asm_name', 'submitter',
       'gbrs_paired_asm', 'paired_asm_comp', 'ftp_path',
       'excluded_from_refseq', 'relation_to_type_material',
       'asm_not_live_date']

a = pd.read_csv(summ_url, sep="\t",skiprows=2,
        names=acols, low_memory=False)


#################### Pull only the most recent complete genomes ####################
c = a[a.assembly_level == "Complete Genome"]

# Pull out the most recent submission for each species only
rc = c.loc[dict(c.sort_values(by="seq_rel_date", ascending=True).apply(lambda x:(x["species_taxid"], x.name), axis=1).values).values()]
sub = rc.loc[random.sample(list(rc.index), k=nbacs)]

# Add in the strain name
sub["strain_name"] = sub.organism_name.replace(" ", "_") + "_" + sub.infraspecific_name.replace(" ", "_")
sub["path"] = sub.ftp_path.apply(lambda x:"./"+x.split("/")[-1])
sub["abundance"] = 100/len(sub)

######################################## Create the files for nanosim ########################################
out = f"nanosim-{nbacs}_genomes-{nsamples}_samples-{read_count}_reads-"

#################### Write out the GL file ####################
sub[["strain_name", "path"]].to_csv(out +"gl.tsv", sep="\t", index=False)


#################### Make the abundance file ####################
# Make the abundance file
cols = ["Size"]
cols.extend([x for x in range(0,nsamples)])
abund = pd.DataFrame(columns=cols)
abund["Size"] = sub.strain_name.values

if not vary:
    abund[cols[1:]] = 100/len(abund.Size)
else:
    for i in range(0, len(abund.columns)-1):
        abund[i] = random.choices(range(1,1000), k=nbacs)

# Name the columns based off of read_counts
abund.rename(columns=dict([(x, read_count) for x in range(0,nsamples)]), inplace=True)

# Correct the read-counts to make them percentages
abund = abund.apply(lambda x: (x/x.sum()) *100 if x.dtype == "int64" \
           else x)

# Write out the abundance file
abund.to_csv(out + "al.tsv", sep="\t",
             index=None)


#################### Make the dna_type_list file ####################
import wget
import gzip

# Make the proper FTP path
sub["full_ftp_path"] = sub.ftp_path.apply(lambda x:x + "/" + x.split("/")[-1] + "_genomic.fna.gz")

# Download the files
downloads = [wget.download(x) for x in sub.full_ftp_path.values]

dl = sub[["strain_name", "path"]].copy()
dl["header"] = dl.path.apply(lambda x:gzip.open(x, "rt").readline().replace(">", ""))
dl["dna_type"] = "circular"

# Write out the DL file
dl.drop("path", axis=1).to_csv(out + "dl.tsv", sep="\t",
                               index=False, header=False)
