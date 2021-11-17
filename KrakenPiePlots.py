#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Set up parser
import argparse as ap
parser = ap.ArgumentParser()
parser.add_argument("-f", "--file", dest="F",
        help="File containing the kraken report")
parser.add_argument("--labels", dest="labels", default="Name",
        help="The column containing the labels for the slices (default:Name)")
parser.add_argument("--values", dest="values", default="CladeNumFrag",
        help="The column containing the values to turn into a percent for slice sizes")
parser.add_argument("--min-perc", dest="min_perc", default=1, type=int,
        help="Minimum percentage to display in the pie plot")
parser.add_argument("--level", dest="level", default="D",
        help="Level of taxonomy to display (S: species, etc...)")
parser.add_argument("--title", dest="title",
        help="Image title")
parser.add_argument("-o", "--out", dest="out",
        help="File name to save image as (name.png|jpg|tiff)")

args = parser.parse_args()

if not args.F:
    parser.error("Please specify the file containing the kraken2 report. See --help for options.")

# Set up the function
def MakePie(F, labels="Name", values="CladeNumFrag", min_perc=1, level="D", title=None, out=None):
    """
    Function to make a pie plot
    
    Parameters
    ----------
    F: File containing the kraken report
    labels: The column containing the labels for the slices
    values: The column containing the values to turn into a % for slice sizes
    min_perc: Minimum percentage to display in the pie plot
    level: Level of taxonomy to display (S: species, G: genus, F: family, O: order ... C,P,K,D)
    
    """
    
    wdf = pd.read_csv(F, sep="\t",
                header=None,
                names=["CladeFragPerc", "CladeNumFrag", "TaxonFrags",
                       "Rank", "TaxID", "Name"])
    
    # Slice based on the taxonomic level
#     wdf=wdf[(wdf.Rank == level) & (wdf.Rank != "unclassified")]
    wdf = wdf[(wdf.Rank == level) | (wdf.Rank == "U")]
    
    # Fix the spacing
    wdf[labels] = wdf[labels].apply(lambda x:x.strip())
    
    # Group the labels < min_perc into 'Others' category
    wdf[labels] = wdf.apply(lambda x:x[labels] \
                            if x["CladeFragPerc"] > min_perc \
                           else "Other",
                           axis=1)
    
    # Sum up the 'Others' category
    wdf = wdf.groupby(by="Name").sum().reset_index()
    
    # Sort the values for neater pie plots
    wdf.sort_values(by="CladeFragPerc", inplace=True, ascending=False)

    # Make label fontsize normal if not domain-level
    lfs = 10
    
    if level == "D":
        # Set up the colors for domains
        cdict = {"unclassified":"black",
                "Eukaryota":"red",
                "Bacteria":"lightblue",
                "Viruses":"lightgreen",
                "Archaea":"black",
                "Other": "gray"}

        # Make the color list using the color dictionary
        colors = [cdict[x.strip()] for x in wdf[labels]]
        
        # Make label fontsize 0 if domain-level
        lfs=0

    # Create the plot canvas with the right figure size
    fig = plt.figure(1, figsize=(10,10))
    ax = fig.add_axes([0,0,1,1])


    # Only have the unclassified and 'Other' slice explode
    explode = [0.1 if x.strip() == "unclassified" or x.strip() == "Other" \
              else 0 for x in wdf[labels]]
    
    # Create the pie plot
    patches, texts, pct_labels = ax.pie(wdf[values], 
                                labels=wdf[labels], 
                            autopct="%.1f%%", pctdistance=0.8,
#                             colors = colors,
                            labeldistance=1.1, explode=explode,
                           wedgeprops=dict(linewidth=2, edgecolor='k', alpha=0.5),
                           textprops=dict(fontsize=lfs, color="k"))

    # Adjust the labels to make them neater and easier to read
    for l in pct_labels:
        l.set_backgroundcolor("white")
        l.set(bbox=dict(facecolor='w', edgecolor='k', boxstyle="round, pad=0.2"))
        l.set(fontsize=10)
        
    # Adjust the text labels as well
    for t in texts:
        t.set_ha("center")
        t.set_va("center")
        new_name = t._text.replace(" ", "\n")
        t.set_text(new_name)
        t.set_rotation_mode(None)
        t.set_bbox(dict(facecolor='w', edgecolor='k', boxstyle="round, pad=0.2"))
        
        if level == "D":
            t.set_visible(False)
        
        
    # Get the right level name for the legend title
    level_dict = {"S":"Species", "G":"Genera",
                 "F":"Families", "O":"Orders", 
                 "C":"Classes", "P":"Phyla", 
                 "K":"Kingdoms", "D":"Domain"}
    
    if level == "D":
        # Add in the legend
        plt.legend(title=f"{level_dict[level]} > {min_perc}% prevalence\n",
              frameon=True, edgecolor="k", facecolor="white")
        plt.title(title)
    else:
        # Add the title
        plt.title(title + "\n" + f"{level_dict[level]} > {min_perc}% prevalence\n")
    
    # Save the image
    if out:
        plt.savefig(out, dpi=200, bbox_inches="tight")

# Run the function
MakePie(F=args.F, labels=args.labels, values=args.values, 
        min_perc=args.min_perc, level=args.level, 
        title=args.title, out=args.out)
