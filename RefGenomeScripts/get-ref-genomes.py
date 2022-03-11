#!/usr/bin/env python3

import pandas as pd
import re
from tqdm import tqdm
import wget
import multiprocessing as mp

#################### Set up argument parser ####################
import argparse as ap

parser = ap.ArgumentParser(
description="""
Script to download reference genomes from RefSeq based on a specified assembly level.
The assembly_summary file is required.
""", formatter_class=ap.RawTextHelpFormatter
)

parser.add_argument("-s", "--summary", dest="summ",
        required=True,
        help="Assembly summary file")

parser.add_argument("-o", "--out", dest="out",
        required=False,
        help="File to write the ftp files out to.")

#parser.add_argument("-c", "--categories", dest="cats",
        #help="categories.dmp file")

parser.add_argument("-l", "--levels", dest="level",
        required=True,
        help="""Minimum genome assembly level. If multiple, separate with a comma and enclose in quotes.
From most complete to least complete:
    1) Complete Genome
    2) Chromosome
    3) Scaffold
    4) Contig
""")

parser.add_argument("-t", "--threads", dest="threads",
        required=False, default=2, type=int,
        help="Number of threads to download files")

parser.add_argument("--name-dict", dest="ndict", action="store_true",
        help="Option to write out a name_dict for 'fixing' reference names downstream")

#parser.add_argument("--domain", dest="domain",
        #choices=["B", "E", "V", "all"],
        #help="""Domain of interest for the genomes
#B: Bacteria,
#E: Eukaryota,
#V: Viruses
#all: B,E,V, unclassified, and other
#""")

args = parser.parse_args()

# Split up the levels
reg = r"complete genome\b|chromosome\b|scaffold\b|contig\b"
levels = re.findall(reg, args.level)

if len(levels) == 0:
    parser.error("No acceptable assembly levels chosen. See --help for options")

#################### Import the assembly_summary file ####################
cols = ['assembly_accession', 'bioproject', 'biosample', 'wgs_master',
       'refseq_category', 'taxid', 'species_taxid', 'organism_name',
       'infraspecific_name', 'isolate', 'version_status', 'assembly_level',
       'release_type', 'genome_rep', 'seq_rel_date', 'asm_name', 'submitter',
       'gbrs_paired_asm', 'paired_asm_comp', 'ftp_path',
       'excluded_from_refseq', 'relation_to_type_material',
       'asm_not_live_date']

a = pd.read_csv(args.summ, sep="\t",skiprows=2,
        names=cols, low_memory=False)

#################### Pull the genomes based on assembly levels ####################
cols = ["assembly_accession", "taxid", "species_taxid", "strain_name",
        "assembly_level", 'ftp_path', "seq_rel_date"]

a["strain_name"] = a.organism_name + " " + a.infraspecific_name.fillna("")
a["strain_name"] = a.strain_name.apply(lambda x: x.strip())

db = a[a.assembly_level == levels[0].title()][cols].sort_values(by="seq_rel_date", ascending=False).copy()

for l in levels[1:]:
    tmp = a[a.assembly_level == l.title()].copy()[cols]
    uniqs = dict([(i,s) for s,i in dict(tmp.sort_values(by="seq_rel_date", ascending=True)["strain_name"]).items()])
    tmp = tmp.loc[uniqs.values()]
    db = pd.concat([db, tmp.loc[[i for i,s in dict(tmp["strain_name"]).items() \
            if s not in db.strain_name.values]]]).reset_index(drop=True)

#################### Print out the results ####################
print(f"There are {'{:,}'.format(len(db.species_taxid.unique()))} unique species in the database")
print(f"There are {'{:,}'.format(len(db.strain_name.unique()))} unique organisms in the database\n")

print(f"There are {'{:,}'.format(len(db[db.assembly_level == 'Complete Genome']))} complete genome assemblies in the database")
print(f"There are {'{:,}'.format(len(db[db.assembly_level == 'Chromosome']))} chromosome assemblies in the database")
print(f"There are {'{:,}'.format(len(db[db.assembly_level == 'Scaffold']))} scaffold assemblies in the database")
print(f"There are {'{:,}'.format(len(db[db.assembly_level == 'Contig']))} contig assemblies in the database\n")

#################### Write out the name dict ####################
if args.ndict:
    import pickle
    name_dict = dict(db.apply(lambda x:(x["ftp_path"].split("/")[-1], "|".join(list(x[["strain_name", "species_taxid", "kingdom", "assembly_level"]])).replace(" ", "_")), 
                                    axis=1).values)
    with open('name_dict.pkl', 'wb') as handle:
        pickle.dump(name_dict, handle)

#################### Download the files ####################
def main(url):
    wget.download(url, bar=None)

#urls = [f"{x}/{x.split('/')[-1]}_genomic.fna.gz\n" for x in db.ftp_path.values]
urls = ["https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/856/265/GCF_000856265.1_ViralProj14966/GCF_000856265.1_ViralProj14966_genomic.fna.gz", 
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/845/205/GCF_000845205.1_ViralProj14549/GCF_000845205.1_ViralProj14549_genomic.fna.gz",
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/866/405/GCF_000866405.1_ViralProj15491/GCF_000866405.1_ViralProj15491_genomic.fna.gz",
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/861/165/GCF_000861165.1_ViralProj15288/GCF_000861165.1_ViralProj15288_genomic.fna.gz",
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/865/725/GCF_000865725.1_ViralMultiSegProj15521/GCF_000865725.1_ViralMultiSegProj15521_genomic.fna.gz"]

pool = mp.Pool(args.threads)
inp = input(f"There are {len(urls)} files to download. Would you like to proceed? [y/n] ")
if inp.lower() == "y":
    for _ in tqdm(pool.imap_unordered(main, urls), total=len(urls)):
        pass
else:
    print("Files were not downloaded!")
