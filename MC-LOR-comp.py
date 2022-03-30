#!/usr/bin/env python3

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

import random
from collections import Counter
import re

#################### Add in argument parser ####################
import argparse as ap
parser = ap.ArgumentParser(description="""
Script to examine differences in Log-Odds Ratios of specific groups within two different populations.
This script was useful for identifying differences in the proportions of samples 
within a set of groups between two different populations.
""", formatter_class=ap.RawTextHelpFormatter)

# Input file, group, and population argument
parser.add_argument("-f", "--file", dest="F", required=True,
        help="File containing the data. Rows are samples and columns are variables")
parser.add_argument("--pop-col", dest="pop_col", required=True,
        help="Column in given file that should be used to separate the populations")
parser.add_argument("--group-col", dest="group_col", required=True,
        help="Column in given file that should be used to separate the groups within the populations")
parser.add_argument("--pop-order", dest="pop_order",
        help="Order of populations for plotting. First group is the numerator and appears with the first color supplied.")


# Monte Carlo and bootstrapping parameters
parser.add_argument("--niter", dest="niter", type=int, default=500,
        help="Number of iterations to perform the bootstrapping (default: 500)")
parser.add_argument("-k", "--k-samples", dest="k", type=int, default=200,
        help="Number of samples to bootstrap (default: 200)")
parser.add_argument("-a", "--alpha", dest="alpha", type=float, default=0.05,
        help="Significance cut-off for Monte Carlo estimation (default: 0.05)")

# Plot arguments
parser.add_argument("--title", dest="title", default=None,
        help="Title for generated plot")
parser.add_argument("--colors", dest="colors", type=str,
        default="green, gray, steelblue",
        help="""Specify the colors to use for plotting as a string separated by a comma
        NOTE: color 1 is for the elevated groups, color 2 is for groups that are similar, and color 3 is for groups that are elevated in population 2""")
parser.add_argument("-o", "--out", dest="out", default=None,
        help="Specify image name *.png|tiff|jpeg|pdf. If no name specified, then image will not be saved")

args = parser.parse_args()

# Make sure the necessary argument are in the right format
args.colors = [x.strip() for x in args.colors.split(",")]


if args.out and not re.search("png|tiff|pdf|jpe*g", args.out):
    parser.error("Incorrect image format specified. Only accept png|tiff|jpeg|pdf")

if not 0 < args.alpha < 1:
    parser.error("Invalid alpha value. Must be between 0 and 1.")


#################### Set up functions ####################
def MCComp(f, group_col, var_col, niter=500, k=200, sep="\t",
          plot_title="", figsize=(10,10), group_order=None, alpha=0.05, colors=None):
    """
    Script to run Monte Carlo Simulation to compare two groups
    
    Parameters
    ----------
    f: TSV/CSV containing the data
    group_col: Column containing the label to split the samples into two groups
    var_col: Column containing the variable of interest
    niter: Number of iterations to run the simulation
    k: Size of bootstrapped sample population
    sep: Separator for file containing the data
    plot_title: Title for plot
    figsize: Tuple of width,length
    group_order: List of the two groups for comparison. First member will be in the numerator for the comparisons
    xlabel: Label for the x-axis
    colors: List of 3 colors to use for the plots [higher in g1, similar in both, higher in g2]
    """
    
    m = pd.read_csv(f, sep=sep)
    groups = m[group_col].unique()
    if len(groups) == 2:
        # plt.figure(figsize=figsize, facecolor="white", edgecolor="k")
        # plt.title(plot_title)

        main(f, group_col, var_col, niter, k, sep, 
             group_order, alpha, colors)
    else:
        print("More than 2 groups are present!")
        print([x for x in groups])

def main(f, group_col, var_col, niter=100, k=100, sep="\t", 
        group_order=None, alpha=0.05, colors=None, plot_title=None):
    
    # Import the data
    m = pd.read_csv(f, sep=sep, dtype=str)
    
    # Remove the groups with counts < 10 since it messes up the estimation
    low_vars = [k for k,v in dict(m[var_col].value_counts()).items() if v < 10]
    low_index = []
    low_counts = [low_index.extend(list(m[m[var_col] == x].index.values)) \
                  for x in low_vars]
    m.drop(low_index, inplace=True)
    
    print(f"{len(low_vars)} group(s) dropped due to having counts < 10:")
    [print(x) for x in low_vars]
    
    # Split up the groups
    if group_order:
        g1,g2 = group_order
        print(g1,g2)
    else:
        g1,g2 = m[group_col].unique()
        
    n = m[m[group_col] == g1].copy()
    y = m[m[group_col] == g2].copy()
    
    global res
    res = dict([(x,[]) for x in m.dropna()[var_col].unique() if x != "U"])
    nres = dict([(x,[]) for x in m.dropna()[var_col].unique() if x != "U"])

    for i in range(1,niter+1):
        global g1_hat
        g1_hat = dict(Counter(random.choices(n[var_col].dropna().values, k=k)).most_common())
        g2_hat = dict(Counter(random.choices(y[var_col].replace("U", np.nan).dropna().values, k=k)).most_common())

        # Add in the results
        for x in res:
            try:
                g1_hat[x] /= k
            except:
                g1_hat[x] = 0.1
            try:
                g2_hat[x] /= k
            except:
                g2_hat[x] = 0.1

            nres[x].append(np.log2(g1_hat[x]) - np.log2(g2_hat[x]))
            
            if g1_hat[x] > g2_hat[x]:
                res[x].append(g1)
            elif g1_hat[x] < g2_hat[x]:
                res[x].append(g2)
            else:
                res[x].append("Same")

    ret = []
    sigs = []
    ia = 1-alpha # Inverse of the alpha
    for x in res.keys():
        r = dict(Counter(res[x]))
        ret.append([f"{x}: {x1} ({r[x1]}/{niter})" for x1 in r.keys()])
        sigs.extend([x for x1 in r.keys() if r[x1] > ia*niter])
        
    [print(x) for x in ret]

    # Create the values DF
    df = pd.DataFrame(nres)
    vals = []

    for col in df.columns:
        vals.extend([(x,col) for x in df[col].values])

    # global vdf
    vdf = pd.DataFrame(vals, columns=["LOR", "CC"])
    cc_stats = dict([(x, np.mean(vdf[vdf.CC == x]["LOR"])) for x in vdf.CC.unique()])

    vdf["Status"] = vdf.CC.apply(lambda x:f"Higher in {g1}" if cc_stats[x] > 0.2 \
              else (f"Higher in {g2}" if cc_stats[x] < -0.2 \
                   else "Similar"))
    
    # Massage the data
    # global uv,lv,nv
    uv = dict([(cc, [vdf[vdf.CC == cc]["LOR"].median(), \
                vdf[vdf.CC == cc].sort_values(by="LOR", ascending=False).reset_index(drop=True).iloc[int(-0.05 * niter)]["LOR"], \
            vdf[vdf.CC == cc].sort_values(by="LOR", ascending=False).reset_index(drop=True).iloc[int(0.05 * niter)]["LOR"]]) \
            for cc in vdf.CC.unique() if vdf[vdf.CC == cc]["LOR"].median() > 0.2])

    lv = dict([(cc, [vdf[vdf.CC == cc]["LOR"].median(), \
                    vdf[vdf.CC == cc].sort_values(by="LOR", ascending=False).reset_index(drop=True).iloc[int(-0.05 * niter)]["LOR"], \
                vdf[vdf.CC == cc].sort_values(by="LOR", ascending=False).reset_index(drop=True).iloc[int(0.05 * niter)]["LOR"]]) \
                for cc in vdf.CC.unique() if vdf[vdf.CC == cc]["LOR"].median() < -0.2])

    nv = dict([(cc, [vdf[vdf.CC == cc]["LOR"].median(), \
                    vdf[vdf.CC == cc].sort_values(by="LOR", ascending=False).reset_index(drop=True).iloc[int(-0.05 * niter)]["LOR"], \
                vdf[vdf.CC == cc].sort_values(by="LOR", ascending=False).reset_index(drop=True).iloc[int(0.05 * niter)]["LOR"]]) \
                for cc in vdf.CC.unique() if 0.2 > vdf[vdf.CC == cc]["LOR"].median() > -0.2])


    cols = ["MedLOR", "L95", "U95"]
    
    # global udf,ldf,ndf

    udf = pd.DataFrame.from_dict(uv, orient="index", 
                          columns=cols)
    udf["Lerr"] = abs(udf.MedLOR - udf.L95)
    udf["Uerr"] = abs(udf.MedLOR - udf.U95)
    udf["orig"] = "udf"

    ldf = pd.DataFrame.from_dict(lv, orient="index", 
                          columns=cols)
    ldf["Lerr"] = abs(ldf.MedLOR - ldf.L95)
    ldf["Uerr"] = abs(ldf.MedLOR - ldf.U95)
    ldf["orig"] = "ldf"

    ndf = pd.DataFrame.from_dict(nv, orient="index", 
                          columns=cols)
    ndf["Lerr"] = abs(ndf.MedLOR - ndf.L95)
    ndf["Uerr"] = abs(ndf.MedLOR - ndf.U95)
    ndf["orig"] = "ndf"
    
    # Combine the dataframe to sort by median values
    # global tot_df
    tot_df = pd.concat([udf, ndf, ldf]).sort_values(by="MedLOR", ascending=False)
    
    # resplit the sorted data
    udf = tot_df[tot_df.orig == "udf"]
    ndf = tot_df[tot_df.orig == "ndf"]
    ldf = tot_df[tot_df.orig == "ldf"]
    

    # Make the plot
    # plt.style.use("bmh")
    
    if not colors:
        colors = ["green", "gray", "steelblue"]


    # Create the plot
    plt.errorbar(x=udf.index, y=udf.MedLOR, yerr=np.array([udf.Lerr.values, udf.Uerr.values]),
                fmt="_", color=colors[0],
                capsize=10, markersize=30)
    sns.scatterplot(x=udf.index, y=udf.MedLOR, color=colors[0], marker="D", edgecolor="k")

    plt.errorbar(x=ndf.index, y=ndf.MedLOR, yerr=np.array([ndf.Lerr.values, ndf.Uerr.values]), 
            fmt="_", color=colors[1],
                capsize=10, markersize=30)
    sns.scatterplot(x=ndf.index, y=ndf.MedLOR, color=colors[1], zorder=2, marker="D", edgecolor="k")

    plt.errorbar(x=ldf.index, y=ldf.MedLOR, yerr=np.array([ldf.Lerr.values, ldf.Uerr.values]),
                fmt="_", color=colors[2], markersize=30,
                capsize=10)
    sns.scatterplot(x=ldf.index, y=ldf.MedLOR, color=colors[2], zorder=2, marker="D", edgecolor="k")
    
    # Add in the labels and title
    plt.ylabel("Log\N{SUBSCRIPT TWO} Odds Ratio")
    plt.xlabel(f"{var_col}")

    # Add in 'line of no effect'
    plt.axhline(0, color="darkgray", linestyle=":")
    
    # Add in the points for significance
    xpos = 0
    colx = [] # List of xticks to color
    for s in udf.index:
        if s in sigs:
            plt.text(s="*", x=xpos, y=udf.loc[s, "U95"],
                    ha="center", va="bottom",
                    size=14, color=colors[0])
            colx.append(xpos)
        xpos += 1
        
    for s in ndf.index:
        if s in sigs:
            plt.text(s="*", x=xpos, y=ndf.loc[s, "U95"],
                    ha="center", va="bottom",
                    size=14, color=colors[1])
            colx.append(xpos)
        xpos += 1
        
    for s in ldf.index:
        if s in sigs:
            plt.text(s="*", x=xpos, y=ldf.loc[s, "U95"],
                    ha="center", va="bottom",
                    size=14, color=colors[2])
            colx.append(xpos)
        xpos += 1
        
    # Add in a legend for the colors
    custom_lines = [Line2D([0], [0], color=colors[0], lw=2),
                    Line2D([0], [0], color=colors[1], lw=2),
                    Line2D([0], [0], color=colors[2], lw=2)]

    labels = [f"Higher in {g1}",
              "Similar",
              f"Higher in {g2}"]
    
    # Add in legend
    plt.legend(custom_lines, labels,
              facecolor="white", frameon=True, edgecolor="k")
    
    # Color the xticks red if they are significantly different
    for x in colx:
        plt.gca().get_xticklabels()[x].set_color("red")
        
    # Add in the text box for the p-val
    ymin,ymax = plt.ylim()
    plt.text(x=0, y=ymin*0.9, s=f"* indicates an empirical P value < {alpha}",
        va="bottom", ha="left", 
         bbox=dict(facecolor="wheat", edgecolor="k"))
    
    if not plot_title:
        plt.title(f"Comparison of {group_col} in the {g1} and {g2} populations")

#################### Run the functions ####################

if args.F.split(".")[-1].lower() == "csv":
    sep = ","
elif args.F.split(".")[-1].lower() == "tsv":
    sep = "\t"
else:
    parser.error("Only CSV/TSV files are accepted")

plt.figure(figsize=(8,8))
MCComp(f=args.F, sep=sep, group_col="Study", var_col=args.group_col, niter=args.niter, k=args.k,
      plot_title=args.title, alpha=args.alpha, colors=args.colors)

if args.out:
    plt.savefig(args.out, dpi=300, bbox_inches="tight")
else:
    plt.show()
