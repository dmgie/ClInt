from cProfile import label
import json
import matplotlib.pyplot as plt
import numpy as np
import re
from natsort import natsorted
from statsmodels.stats.multitest import multipletests

variant_data = {}

with open("./output_all_variants.json", "r") as json_file:
    variant_data = json.load(json_file)


## Extract keys and p_values
keys = list(variant_data.keys())
keys_sorted = natsorted(keys, key=lambda x: int(re.search(r'chr(\d+)', x).group(1)))

## Extract p-values to list
p_values = [value['p_value'] for key, value in variant_data.items()]

## FDR correction by Benjamini Hochberg method
reject, corrected_p_values, _, _ = multipletests(p_values, method='fdr_bh')

for i, (key, value) in enumerate(variant_data.items()):
    value['fdr_corrected'] = corrected_p_values[i]
    value['reject'] = reject[i]


unique_chromosomes = np.unique([chrom.split(":")[0] for chrom in keys])

## Create a color mapping dictionary
color_mapping = {chrom: plt.cm.tab10(i) for i, chrom in enumerate(unique_chromosomes)}

## Create figure
fig, axs = plt.subplots(2, 1, figsize=(8, 6))

for chrom_label in keys_sorted:
    p_value = variant_data[chrom_label]["p_value"]
    fdr_value = variant_data[chrom_label]["fdr_corrected"]
    chrom, value = chrom_label.split(":")        
 
    axs[0].scatter(chrom_label, -1*np.log10(p_value), color=color_mapping[chrom], label=chrom, s=25)
    axs[1].scatter(chrom_label, -1*np.log10(fdr_value), color=color_mapping[chrom], label=chrom, s=25)

## Bonferroni
axs[0].axhline(y=-1*np.log10(0.05/len(keys)), color='red', linestyle='--', label='bonferroni correction')
axs[0].set_xlabel('Variant positions per chromosome')
axs[0].set_ylabel('-log10 p-values')
axs[0].set_xticks([])
axs[0].set_title("Bonferroni corrected p-values")

## Benjamini Hochberg
axs[1].axhline(y=-1*np.log10(0.05), color='red', linestyle='--', label='bonferroni correction')
axs[1].set_xlabel('Variant positions per chromosome')
axs[1].set_ylabel('-log10 p-values')
axs[1].set_xticks([])
axs[1].set_title("Benjamini-Hochberg corrected p-values")


legend_handles = [plt.Line2D([0], [0], marker='o', color='w', label=chrom, markerfacecolor=color_mapping[chrom], markersize=10) for chrom in natsorted(unique_chromosomes)]
plt.legend(handles=legend_handles, title='Chromosomes', loc='center left', bbox_to_anchor=(1.05, 1.2), borderaxespad=0.0)
plt.subplots_adjust(hspace=0.3)
plt.suptitle('Variant association with altered gene expression', fontsize=16)
plt.savefig('manhattan_plot.png',bbox_inches='tight')