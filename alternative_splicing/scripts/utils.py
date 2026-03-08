##################################################
#                                                #
#                Import libraries                #
#                                                #
##################################################

import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mpl
import json
import os



##################################################
#                                                #
#                   Parameters                   #
#                                                #
##################################################

#~~~~~~~~~~~~~~~~ rMATS results ~~~~~~~~~~~~~~~~~#

with open("parameters.json") as f:
    config = json.load(f)

only_junctions = config["alternative_splicing"]["only_junctions"]
suffix_files = '.MATS.JC.txt'
if not only_junctions:
    suffix_files = '.MATS.JCEC.txt'

event_types = ["SE", "A3SS", "A5SS", "MXE", "RI"]
columns_to_convert = ['IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_2', 'IncLevel1', 'IncLevel2']
counts_columns = ['IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_2']
inclusion_levels = ['IncLevel1_1', 'IncLevel1_2','IncLevel1_3', 'IncLevel2_1', 'IncLevel2_2', 'IncLevel2_3']


#~~~~~~~~~~~~~~~~~~~ Colors ~~~~~~~~~~~~~~~~~~~~#

color_lavander = '#DBB2D1'
color_magenta = '#ff0090'
color_fuchsia = '#F231FF'
color_violet = '#8D00B8'
color_zaffre = '#0021A3'
color_palblue = '#1741FF'
color_skyblue = '#05C1FF'
color_tiffany = '#9CE8E0'
color_gray95 = '#DEE0E3'
color_event_types = [color_violet, 
                     color_lavander, 
                     color_palblue, 
                     color_skyblue, 
                     color_tiffany]



##################################################
#                                                #
#                   Functions                    #
#                                                #
##################################################

#~~~~~~~~~~~~~~~~ SetPlotParams ~~~~~~~~~~~~~~~~~#

def SetPlotParams(magnification=1.0, 
                  ratio=float(2.2/2.7), 
                  height=None, 
                  width=None, 
                  fontsize=11., 
                  ylabelsize=None, 
                  xlabelsize=None, 
                  lines_w=1.5, 
                  axes_lines_w=0.7, 
                  legendmarker=True, 
                  tex=False, 
                  autolayout=True, 
                  handlelength=1.5):
    
    if (ylabelsize==None):
        ylabelsize = fontsize
    if (xlabelsize==None):
        xlabelsize = fontsize

    ratio = ratio  # usually this is 2.2/2.7
    fig_width = 2.9 * magnification # width in inches
    fig_height = fig_width*ratio  # height in inches
    if height!=None:
        fig_height = height
    if width!=None:
        fig_width = width
    fig_size = [fig_width,fig_height]
    plt.rcParams['figure.figsize'] = fig_size
    plt.rcParams['figure.autolayout'] = autolayout
    plt.rcParams['lines.linewidth'] = lines_w
    plt.rcParams['lines.markeredgewidth'] = 1.
    plt.rcParams['errorbar.capsize'] = 1 #1.5
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['legend.frameon'] = False
    plt.rcParams['legend.numpoints'] = 1
    plt.rcParams['legend.markerscale'] = 1
    plt.rcParams['legend.handlelength'] = handlelength
    plt.rcParams['legend.labelspacing'] = 0.3
    plt.rcParams['legend.columnspacing'] = 0.3
    if legendmarker==False:
        plt.rcParams['legend.numpoints'] = 1
        plt.rcParams['legend.markerscale'] = 0
        plt.rcParams['legend.handlelength'] = 0
        plt.rcParams['legend.labelspacing'] = 0  
    plt.rcParams['legend.fontsize'] = fontsize
    plt.rcParams['axes.facecolor'] = '1'
    plt.rcParams['axes.edgecolor'] = '0.0'
    plt.rcParams['axes.linewidth'] = axes_lines_w
    plt.rcParams['grid.color'] = '0.85'
    plt.rcParams['grid.linestyle'] = '-'
    plt.rcParams['grid.linewidth'] = axes_lines_w
    plt.rcParams['grid.alpha'] = '1.'
    plt.rcParams['axes.labelcolor'] = '0'
    plt.rcParams['axes.labelsize'] = fontsize
    plt.rcParams['axes.titlesize'] = fontsize
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False
    plt.rcParams['xtick.labelsize'] = xlabelsize
    plt.rcParams['ytick.labelsize'] = ylabelsize
    plt.rcParams['xtick.color'] = '0'
    plt.rcParams['ytick.color'] = '0'
    plt.rcParams['xtick.major.size'] = 3.
    plt.rcParams['xtick.major.width'] = axes_lines_w
    plt.rcParams['xtick.minor.size'] = 0
    plt.rcParams['ytick.major.size'] = 3.
    plt.rcParams['ytick.major.width'] = axes_lines_w
    plt.rcParams['ytick.minor.size'] = 0
    plt.rcParams['xtick.major.pad']= 5.
    plt.rcParams['ytick.major.pad']= 5.
    plt.rcParams['text.usetex'] = tex
    mpl.rc('text', usetex = True)
    mpl.rc('text.latex', preamble=r'\usepackage{sfmath}')
    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['axes.spines.top'] = False
    mpl.rcParams['axes.spines.left'] = True
    mpl.rcParams['axes.spines.bottom'] = True


#~~~~~~~~~~~~~~~~ PairedFlagplot ~~~~~~~~~~~~~~~~#

def PairedFlagplot(y1, y2, 
                   name_cond1='cond1', 
                   name_cond2='cond2', 
                   title=None, 
                   xlabel=None, 
                   magnification=1.1, 
                   bar_height=0.8, 
                   ratio=(2.2/4.3), 
                   savefig=None, 
                   pathfig=None):

    SetPlotParams(magnification=magnification, 
                  ratio=ratio, 
                  handlelength=0.8, 
                  fontsize=10)

    fig, ax = plt.subplots()

    conditions = [name_cond1, name_cond2]
    colors = color_event_types
    for i, el_1_2 in enumerate(zip(y1, y2)):
        color = colors[i]
        ax.barh(conditions, width=el_1_2, height=bar_height, color=color, label=event_types[::-1][i])
    if not xlabel:
        ax.set_xlabel('%% events')
    else:
        ax.set_xlabel(xlabel)
    ax.set_xticks([0, 50, 100])
    ax.set_xlim([0, 100])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    if title:
        ax.set_title(title)
    if savefig:
        plt.savefig(pathfig, bbox_inches='tight')# dpi=1000,
    #plt.show()
    plt.close(fig)


#~~~~~~~~~~~~~ SignifEventsBarplot ~~~~~~~~~~~~~~#

def SignifEventsBarplot(l_signif_perc, 
                        labels, 
                        title=None, 
                        magnification=1.2,
                        ratio=(2.2/3.2), 
                        ylabel='\% events', 
                        savefig=None, 
                        pathfig=None):

    SetPlotParams(magnification=magnification, 
                  ratio=ratio, 
                  handlelength=0.8, 
                  fontsize=10)
    
    fig, ax = plt.subplots()

    l1 = [100-x for x in l_signif_perc]
    l2 = l_signif_perc
    
    # Bar positions on the x-axis
    x = np.arange(len(l1))

    # Plot the bars
    ax.bar(x, l2, bottom=l1, label='signif.', color=color_zaffre)
    ax.bar(x, l1, label='not signif.', color=color_gray95)

    # Add labels and title
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45)
    ax.set_ylabel(ylabel)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylim([0, 100])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    
    if title:
        ax.set_title(title)
    if savefig:
        plt.savefig(pathfig, bbox_inches='tight')# dpi=1000,
    #plt.show()
    plt.close(fig)


#~~~~~~~~~~~~~~~~~ VolcanoPlot ~~~~~~~~~~~~~~~~~~#

def VolcanoPlot(x, 
                y_pval, 
                x_th, 
                y_pval_th=0.05, 
                name_cond1='cond1', 
                name_cond2='cond2', 
                title=None, 
                xlabel=None, 
                ylabel=None, 
                magnification=1., 
                ratio=(2.5/2.), 
                s=3, 
                color=color_magenta, 
                color_cond1=color_palblue, 
                color_cond2=color_fuchsia, 
                two_colors=True, 
                y_upper_lim=15, 
                color_nonsignif=color_gray95, 
                savefig=None, 
                pathfig=None):

    SetPlotParams(magnification=magnification, 
                  ratio=ratio, 
                  handlelength=0.8, 
                  fontsize=10)
    
    fig, ax = plt.subplots()

    mask_notnull = pd.notnull(y_pval) & pd.notnull(x)
    y_pval = y_pval[mask_notnull]
    x = x[mask_notnull]
    min_pval = min(y_pval[y_pval>0])
    y_pval_visible_th = 10**(-y_upper_lim)
    y_pval[y_pval<y_pval_visible_th] = 10**(-y_upper_lim+0.1) #min_pval
    mask_nonsignif = abs(x)<x_th
    mask_nonsignif = mask_nonsignif | (y_pval>y_pval_th)
    x_nonsignif = x[mask_nonsignif]
    y_pval_nonsignif = y_pval[mask_nonsignif]
    neg_log_y_pval_nonsignif = -np.log10(y_pval_nonsignif)
    ax.scatter(x_nonsignif, neg_log_y_pval_nonsignif,
               color=color_nonsignif,
               s=s,
               alpha=1)
    x_signif = x[mask_nonsignif==False]
    y_pval_signif = y_pval[mask_nonsignif==False]
    neg_log_y_pval_signif = -np.log10(y_pval_signif)
    if two_colors:
        mask_pos = x_signif>=0
        x_pos = x_signif[mask_pos]
        y_pval_pos = y_pval_signif[mask_pos]
        neg_log_y_pval_pos = -np.log10(y_pval_pos)
        ax.scatter(x_pos, neg_log_y_pval_pos,
                color=color_cond1,
                s=s,
                alpha=1)
        mask_neg = x_signif<0
        x_neg = x_signif[mask_neg]
        y_pval_neg = y_pval_signif[mask_neg]
        neg_log_y_pval_neg = -np.log10(y_pval_neg)
        ax.scatter(x_neg, neg_log_y_pval_neg,
                color=color_cond2,
                s=s,
                alpha=1)
    else:
        ax.scatter(x_signif, neg_log_y_pval_signif,
                color=color,
                s=s,
                alpha=1)

    ax.axhline(-np.log10(y_pval_th), lw=0.7, ls='--', color='black')
    ax.axvline(x_th, lw=0.7, ls='--', color='black')
    ax.axvline(-x_th, lw=0.7, ls='--', color='black')
    ax.set_yticks([0, 5, 10, 15])
    ax.set_ylim([0, y_upper_lim])
    if not xlabel:
        plt.xlabel('$x^{\mathrm{%s}}-x^{\mathrm{%s}}$'%(name_cond1, name_cond2))
    else:
        plt.xlabel(xlabel)
    if not ylabel:
        plt.ylabel('-log(p-adj)')
    else:
        plt.ylabel(ylabel)
    if title:
        plt.title(title)
    #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    if savefig:
        plt.savefig(pathfig, bbox_inches='tight')# dpi=1000,
    #plt.show()
    plt.close(fig)


#~~~~~~~~~~~~~~~~~ export_events_by_gene_type_table ~~~~~~~~~~~~~~~~~~#

def export_events_by_gene_type_table(results_dict_gene_summary, path_gene_summary_output):
    
    # Build long format table
    rows = []
    for event_type, df in results_dict_gene_summary.items():
        for gene_type in df["type"].unique():
            mask = df["type"] == gene_type
            n_sig = df.loc[mask, "num_signif"].sum()
            n_tot = df.loc[mask, "num_events"].sum()
            perc = 100 * n_sig / n_tot if n_tot != 0 else 0

            rows.append({
                "event_type": event_type,
                "subrow": "tot events",
                "gene_type": gene_type,
                "value": "%d"%n_tot
            })
            rows.append({
                "event_type": event_type,
                "subrow": "signif events",
                "gene_type": gene_type,
                "value": "%d"%n_sig
            })
            rows.append({
                "event_type": event_type,
                "subrow": "perc signif",
                "gene_type": gene_type,
                "value": f"{perc:.2f}%"
            })
    df_long = pd.DataFrame(rows)

    # Pivot to MultiIndex rows, single-level columns
    df_pivot = df_long.pivot_table(
        index=["event_type", "subrow"],
        columns="gene_type",
        values="value",
        aggfunc="first"
    ).sort_index()
    df_pivot.index.set_names(["event type", ""], inplace=True)

    # Reorder columns according to desired order
    desired_order = [
        "protein_coding",
        "lncRNA",
        "transcribed_processed_pseudogene",
        "transcribed_unitary_pseudogene",
        "transcribed_unprocessed_pseudogene"
    ]

    # Get the columns that are actually present
    present_cols = df_pivot.columns.tolist()

    # Create final column order: first the desired_order that exist, then the rest
    ordered_cols = [c for c in desired_order if c in present_cols] + \
                [c for c in present_cols if c not in desired_order]

    # Apply MultiIndex columns with top = "gene type", bottom = gene name
    df_pivot = df_pivot[ordered_cols]
    df_pivot.columns = pd.MultiIndex.from_tuples(
        [("gene type", col) for col in df_pivot.columns]
    )

    file_name = path_gene_summary_output / "events_by_gene_type.xlsx"
    df_pivot.to_excel(file_name, merge_cells=True)


#~~~~~~~~~~~~~~~~~ load_or_create_gene_info_df ~~~~~~~~~~~~~~~~~~#

def load_or_create_gene_info_df(gtf_path, path_gene_info):
    """
    Load gene info DataFrame (gene_id, gene_name, gene_type, gene_length, chrom, strand, start, end)
    from CSV if it exists; otherwise create from GTF and save to CSV.

    Parameters:
        gtf_path (str): Path to the GTF file.
        path_gene_info (str): Path to the CSV file to load/save gene info.

    Returns:
        pd.DataFrame: Gene info with columns:
            ['symbol', 'id', 'type', 'chrom', 'strand', 'start', 'end', 'length']
    """
    if os.path.exists(path_gene_info):
        df = pd.read_csv(path_gene_info)
        print(f"Loaded gene info from {path_gene_info}")
    else:
        df = get_gene_info_df(gtf_path)
        df.to_csv(path_gene_info, index=False)
        print(f"Created gene info file and saved to {path_gene_info}")

    return df


#~~~~~~~~~~~~~~~~~ get_gene_info_df ~~~~~~~~~~~~~~~~~~#

def get_gene_info_df(gtf_path):
    """
    Parse a GTF file and return a DataFrame with gene_id, gene_name, gene_type,
    gene_length, and genomic position (chromosome, strand, start, end).

    Parameters:
        gtf_path (str): Path to the GTF file.

    Returns:
        pd.DataFrame: Columns = [
            'gene_id', 'gene_name', 'gene_type', 'gene_length',
            'chrom', 'strand', 'start', 'end'
        ]
    """
    # Load GTF (skip comments)
    gtf = pd.read_csv(
        gtf_path,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "chrom", "source", "feature", "start", "end",
            "score", "strand", "frame", "attribute"
        ]
    )

    # Keep only 'gene' entries
    genes = gtf[gtf["feature"] == "gene"].copy()

    # Extract attributes using regex
    genes["symbol"] = genes["attribute"].str.extract('gene_name "([^"]+)"')
    genes["id"]   = genes["attribute"].str.extract('gene_id "([^"]+)"')
    genes["type"] = genes["attribute"].str.extract('gene_type "([^"]+)"')

    # Fallback to gene_biotype if gene_type is missing
    missing_mask = genes["type"].isna()
    if missing_mask.any():
        genes.loc[missing_mask, "type"] = genes.loc[missing_mask, "attribute"].str.extract('gene_biotype "([^"]+)"')

    # Compute gene length
    genes["length"] = genes["end"] - genes["start"] + 1

    # Clean gene_id (remove version suffix like ".1", ".2", etc.)
    genes["id"] = genes["id"].str.replace(r"\.\d+$", "", regex=True)

    # Drop duplicates — keep the longest version if multiple
    genes = genes.sort_values("length", ascending=False).drop_duplicates("id")

    # Keep relevant columns
    genes = genes[["symbol", "id", "type", "chrom", "strand", "start", "end", "length"]]

    # Reset index
    genes.reset_index(drop=True, inplace=True)

    return genes


#~~~~~~~~~~~~~~~~~ summarize_splicing_per_gene ~~~~~~~~~~~~~~~~~~#

def summarize_splicing_per_gene(rmats_result, avg_reads_columns, event_type, avg_reads_threshold=10, 
                                  FDR_threshold=0.05, abs_delta_psi_threshold=0.1, gene_info_df=None):
    """
    Summarize splicing events per gene, including significance filtering.

    Parameters:
        rmats_result (pd.DataFrame): rMATS results table.
        avg_reads_columns (list of str): columns with read counts to compute average coverage.
        event_type (str): one of 'SE', 'RI', 'A5SS', 'A3SS', 'MXE'
        avg_reads_threshold (float): minimum average reads to consider significant.
        FDR_threshold (float): maximum FDR to consider significant.
        abs_delta_psi_threshold (float): minimum |ΔPSI| to consider significant.
        gene_info_df (pd.DataFrame, optional): mapping gene_name -> gene length and position.

    Returns:
        pd.DataFrame: Summary per gene with metrics and lists of ΔPSI for significant and non-significant events.
    """
    df = rmats_result.copy()

    # Compute average reads
    avg_reads = df[avg_reads_columns].mean(axis=1)

    # Replace zero FDRs to avoid -inf in log
    df["FDR"] = df["FDR"].replace(0, 1e-300)

    # Absolute delta PSI
    df["abs_psi"] = df["IncLevelDifference"].abs()

    # Weighted score
    df["weighted_score"] = df["abs_psi"] * -np.log10(df["FDR"])

    # Determine significance
    mask = pd.Series(True, index=df.index)
    mask &= avg_reads > avg_reads_threshold
    mask &= df["FDR"] < FDR_threshold
    mask &= df["abs_psi"] > abs_delta_psi_threshold
    df["significant"] = mask

    # Group by geneSymbol
    #summary = df.groupby("geneSymbol").apply(summarize_gene_group).reset_index()
    cols_to_use = [c for c in df.columns if c != "geneSymbol"]
    summary = df.groupby("geneSymbol")[cols_to_use].apply(lambda g: summarize_gene_group(g, event_type)).reset_index()

    # Add gene length if provided
    if gene_info_df is not None:
        # Clean up potential version suffixes in gene IDs for matching
        gene_info_df["symbol"] = gene_info_df["symbol"].astype(str)
        summary["geneSymbol"] = summary["geneSymbol"].astype(str)
        summary = gene_info_df.merge(
            summary,
            how="inner",
            left_on="symbol",
            right_on="geneSymbol"
        )
        # Drop redundant column after merge
        summary.drop(columns=["geneSymbol"], inplace=True)

    return summary


#~~~~~~~~~~~~~~~~~ summarize_gene_group ~~~~~~~~~~~~~~~~~~#

def summarize_gene_group(g, event_type):
    """
    Summarize a single gene group (DataFrame).
    Returns a pd.Series with counts, ΔPSI lists, significance, direction, and genome positions.
    """
    signif = g[g["significant"]]
    nonsig = g[~g["significant"]]

    delta_psi_signif = list(signif["IncLevelDifference"])
    num_signif = len(signif)

    # Gene-level significance
    is_significant = num_signif >= 1

    # Sums and averages
    if num_signif == 0:
        sum_abs_psi_signif = None
        avg_abs_psi_signif = None
        sum_weighted_score_signif = None
        avg_weighted_score_signif = None
        direction_delta_psi_signif = None
    else:
        sum_abs_psi_signif = signif["abs_psi"].sum()
        avg_abs_psi_signif = signif["abs_psi"].mean()
        sum_weighted_score_signif = signif["weighted_score"].sum()
        avg_weighted_score_signif = signif["weighted_score"].mean()
        if all(v > 0 for v in delta_psi_signif):
            direction_delta_psi_signif = "positive"
        elif all(v < 0 for v in delta_psi_signif):
            direction_delta_psi_signif = "negative"
        else:
            direction_delta_psi_signif = "mixed"

    # Genome positions per event type (use original column names)
    if event_type == "SE":
        position_cols = ["exonStart_0base", 
                         "exonEnd", 
                         "upstreamES", 
                         "upstreamEE", 
                         "downstreamES", 
                         "downstreamEE"]
    elif event_type == "RI":
        position_cols = ["riExonStart_0base", 
                         "riExonEnd", 
                         "upstreamES", 
                         "upstreamEE", 
                         "downstreamES", 
                         "downstreamEE"]
    elif event_type in ["A5SS", "A3SS"]:
        position_cols = ["longExonStart_0base", 
                         "longExonEnd", 
                         "shortES", 
                         "shortEE", 
                         "flankingES", 
                         "flankingEE"]
    elif event_type == "MXE":
        position_cols = ["1stExonStart_0base", 
                         "1stExonEnd", 
                         "2ndExonStart_0base", 
                         "2ndExonEnd", 
                         "upstreamES", 
                         "upstreamEE", 
                         "downstreamES", 
                         "downstreamEE"]
    else:
        position_cols = []

    genome_positions = []
    if len(position_cols)>0:
        for col_name in position_cols:
            genome_positions.append(list(signif[col_name]))

    # Build result series
    res = {
        "num_events": len(g),
        "num_signif": num_signif,
        "num_non_signif": len(nonsig),
        "sum_abs_psi_signif": sum_abs_psi_signif,
        "avg_abs_psi_signif": avg_abs_psi_signif,
        "sum_weighted_score_signif": sum_weighted_score_signif,
        "avg_weighted_score_signif": avg_weighted_score_signif,
        "delta_psi_signif": delta_psi_signif,
        "significant": is_significant,
        "direction_delta_psi_signif": direction_delta_psi_signif
    }

    # Add genome positions to the result
    for i, col_name in enumerate(position_cols):
        res[col_name] = genome_positions[i]

    return pd.Series(res)