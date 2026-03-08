##################################################
#                                                #
#                Import libraries                #
#                                                #
##################################################

from pathlib import Path
import numpy as np
import pandas as pd
import json
import argparse
import os
import sys

from utils import (
    # Functions
    PairedFlagplot,
    SignifEventsBarplot,
    VolcanoPlot,
    export_events_by_gene_type_table,
    load_or_create_gene_info_df,
    summarize_splicing_per_gene,
    # Variables
    event_types,
    inclusion_levels,
    suffix_files,
    columns_to_convert
)


##################################################
#                                                #
#                   Parameters                   #
#                                                #
##################################################

with open("parameters.json") as f:
    config = json.load(f)

FDR_threshold = config["alternative_splicing"]["FDR_threshold"]
abs_delta_psi_threshold = config["alternative_splicing"]["abs_delta_psi_threshold"]
avg_reads_threshold = config["alternative_splicing"]["avg_reads_threshold"]

path_figures = Path(config["paths"]["figures_dir"])
path_gene_summary = Path(config["paths"]["gene_summary_dir"])
path_rmats_output_root = Path(config["paths"]["rmats_output_dir"])
path_gtf = Path(config["paths"]["gtf_file"])
path_gene_info = Path(config["paths"]["gene_info_file"])


##################################################
#                                                #
#                  Import file                   #
#                                                #
##################################################

parser = argparse.ArgumentParser(description="Assign command-line arguments to script variables.")
parser.add_argument("--line_name", type=str, required=True, help="Name of the line (e.g., OVCAR8)")
parser.add_argument("--cond1", type=str, required=True, help="Condition 1 (e.g., PTresSR)")
parser.add_argument("--cond2", type=str, required=True, help="Condition 2 (e.g., PTresU)")
args = parser.parse_args()
line_name = args.line_name
cond1 = args.cond1
cond2 = args.cond2
comparison = cond1 + 'vs' + cond2

path_rmats_output = path_rmats_output_root / line_name / comparison
if os.path.isdir(path_rmats_output):
    print('results OK')
    path_figures_subfolder = path_figures / line_name / comparison
    os.makedirs(path_figures_subfolder, exist_ok=True)

    results_dict = {}
    all_counts_columns = {}
    for event_type in event_types:
        all_counts_columns[event_type] = []
        file_path = path_rmats_output / Path(event_type + suffix_files)
        df = pd.read_csv(file_path, delimiter='\t', quotechar='"')
        for column in columns_to_convert:
            df[column] = df[column].apply(lambda x: [float(i) if i != 'NA' else None for i in x.split(',')])        
            max_length = max(df[column].apply(len))        
            for i in range(max_length):
                df[f'{column}_{i+1}'] = df[column].apply(lambda x: x[i] if i < len(x) else None)
                if 'SAMPLE' in column:
                    all_counts_columns[event_type].append(column+'_%d'%(i+1))
            df.drop(column, axis=1, inplace=True)
        results_dict[event_type] = df
else:
    print('results NOT found for %s %s\n'%(line_name, comparison))
    sys.exit(1)

gene_info = load_or_create_gene_info_df(path_gtf, path_gene_info)



##################################################
#                                                #
#               Significant events               #
#                                                #
##################################################

results_dict_significant = {}
for event_type in event_types:
    df = results_dict[event_type]
    avg_reads = np.mean(df[all_counts_columns[event_type]].values, axis=1)
    mask = avg_reads>avg_reads_threshold
    FDR = df['FDR']
    mask = mask & (FDR<FDR_threshold)
    abs_psi = abs(df['IncLevelDifference'])
    mask = mask & (abs_psi>abs_delta_psi_threshold)
    results_dict_significant[event_type] = df.loc[mask, :]

perc_signif_events = []
for event_type in results_dict_significant.keys():
    perc_signif_events.append(100*results_dict_significant[event_type].shape[0]/results_dict[event_type].shape[0])

figure_name = 'barplot_significant_events_%s_%s.pdf'%(comparison, line_name)
fig_export_path = path_figures_subfolder / figure_name
SignifEventsBarplot(perc_signif_events, 
                    results_dict_significant.keys(),
                    magnification=1.2,
                    ratio=(2.2/3.2),
                    title='%s vs %s'%(cond1, cond2),
                    ylabel='\% events - '+line_name, 
                    savefig=True, 
                    pathfig=fig_export_path)



##################################################
#                                                #
#               Psi - volcano plot               #
#                                                #
##################################################

inclusion_level_cond1 = [col for col in inclusion_levels if 'Level1' in col]
inclusion_level_cond2 = [col for col in inclusion_levels if 'Level2' in col]

for event_type in event_types:

    title = '%s - %s' % (line_name, event_type)

    df = results_dict[event_type]
    x = df['IncLevelDifference'].values
    y = df['FDR'].values
    xlabel = '$\Psi^{\mathrm{%s}}-\Psi^{\mathrm{%s}}$' % (cond1, cond2)
    figure_name = 'volcano_deltapsi_%s_%s_%s.pdf'%(comparison, event_type, line_name)
    fig_export_path = path_figures_subfolder / figure_name
    VolcanoPlot(x=x, 
                y_pval=y, 
                x_th=abs_delta_psi_threshold, 
                y_pval_th=FDR_threshold, 
                name_cond1=cond1, 
                name_cond2=cond2, 
                title=title, 
                xlabel=xlabel, 
                savefig=True, 
                pathfig=fig_export_path)



##################################################
#                                                #
#               delta_psi flagplot               #
#                                                #
##################################################

'''Flagplot: n. events delta_psi>threshold'''

number_diffspliced_events_cond1 = []
number_diffspliced_events_cond2 = []

for event_type in event_types:
    df = results_dict_significant[event_type]
    mask1 = df['IncLevelDifference'].values>abs_delta_psi_threshold
    n_events1 = np.sum(mask1)
    average_psi1 = np.mean(df.loc[mask1,inclusion_level_cond1].values)
    mask2 = df['IncLevelDifference'].values<-abs_delta_psi_threshold
    n_events2 = np.sum(mask2)
    average_psi2 = np.mean(df.loc[mask2,inclusion_level_cond2].values)
    #print('%s avg. psi: %s=%.2f, %s=%.2f'%(event_type, cond1, average_psi1, cond2, average_psi2))
    number_diffspliced_events_cond1.append(n_events1)
    number_diffspliced_events_cond2.append(n_events2)
number_diffspliced_events_cond1 = number_diffspliced_events_cond1/sum(number_diffspliced_events_cond1)
number_diffspliced_events_cond2 = number_diffspliced_events_cond2/sum(number_diffspliced_events_cond2)
proportions_events_cond1 = np.cumsum(number_diffspliced_events_cond1)
proportions_events_cond2 = np.cumsum(number_diffspliced_events_cond2)

title = line_name
xlabel = '\% events'+' ($\Delta\Psi>%.1f$)' % abs_delta_psi_threshold
y1 = proportions_events_cond1[::-1]*100
y2 = proportions_events_cond2[::-1]*100

figure_name = 'events_proportions_%s_%s.pdf'%(comparison, title)
fig_export_path = path_figures_subfolder / figure_name
PairedFlagplot(y1, 
               y2, 
               name_cond1=cond1, 
               name_cond2=cond2, 
               title=title, 
               xlabel=xlabel, 
               savefig=True, 
               pathfig=fig_export_path)

print('   -> All figures generated\n')



##################################################
#                                                #
#              gene splicing summary             #
#                                                #
##################################################

path_gene_summary_output = path_gene_summary / line_name / comparison

os.makedirs(path_gene_summary_output, exist_ok=True)

results_dict_gene_summary = {}
for event_type in event_types:
    df = results_dict[event_type]
    results_dict_gene_summary[event_type] = summarize_splicing_per_gene(
                                                rmats_result=df,
                                                avg_reads_columns=all_counts_columns[event_type],
                                                event_type=event_type,
                                                avg_reads_threshold=avg_reads_threshold,
                                                FDR_threshold=FDR_threshold,
                                                abs_delta_psi_threshold=abs_delta_psi_threshold,
                                                gene_info_df=gene_info
                                            )
    
    file_name = path_gene_summary_output / Path(event_type + '.csv')
    df_to_export = results_dict_gene_summary[event_type].sort_values(by='num_signif', ascending=False)
    df_to_export.to_csv(file_name, sep=";", float_format="%.4f", index=False)

export_events_by_gene_type_table(results_dict_gene_summary, path_gene_summary_output)

print('   -> Gene events summary completed\n')
print('Finished!\n')