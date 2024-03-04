import json
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import  fisher_exact
from statsmodels.stats.multitest import multipletests

def get_condition_counts(strings, condition):
    ## Calculate number of patients in certain condition
    total_num = len(strings)
    num = strings.count(condition)

    cond_counts = {
        condition:num,
        "Rest":total_num-num
    }

    return cond_counts

def fisher(vec1, vec2):
    odds_ratio, p_value = fisher_exact([vec1, vec2])
    return p_value

with open('./output_all_variants.json', 'r') as file:
    data = json.load(file)

alpha = 0.05
l = len(data)

## get list of p-values
## FDR correct by Benjamini Hochberg

p_values = [value['p_value'] for key, value in data.items()]
reject, corrected_p_values, _, _ = multipletests(p_values, method='fdr_bh')

for i, (key, value) in enumerate(data.items()):
    value['fdr_corrected'] = corrected_p_values[i]

x = []
y = []
texts = []

significance_dict = {}

## Determine variant association with certain severity condition
#### -> Calculate Fisher test for significant variants (in terms of altered expression)
for key, entry in data.items():
    # if entry["p_value"] <= alpha/l:
    if entry["fdr_corrected"] <= alpha:

        significance_dict[key] = {}

        for condition in ["COVIDsevere", "COVIDmoderate", "COVIDmild", "postCOVID", "COVIDasymptomatic"]:

            conditions_dict_has = get_condition_counts(entry["has_var_cond"], condition)
            conditions_dict_hasnt = get_condition_counts(entry["hasnt_var_cond"], condition)

            counts_has = [conditions_dict_has[condition], conditions_dict_has["Rest"]]
            counts_hasnt = [conditions_dict_hasnt[condition], conditions_dict_hasnt["Rest"]]

            p = fisher(counts_has, counts_hasnt)

            significance_dict[key][condition] = p

        x.append(-1*entry["log_fc"])
        y.append(-1*math.log(entry["p_value"]))
        texts.append(key)


data_matrix = np.array([[significance_dict[row][col] for col in (significance_dict[row])] for row in sorted(significance_dict)])


## Heatmap plot

fig, axs = plt.subplots(1,2, figsize=(10, 6), gridspec_kw={'width_ratios': [0.8, 1]})
axs[0].scatter(x, y)
axs[0].grid(True)
axs[0].set_xlabel('log2 fold-change has variant / has no variant')
axs[0].set_ylabel('-log10 p-value')

scatter = axs[0].scatter(x, y, c='blue', alpha=0.75, edgecolors='w', s=100)

for i, txt in enumerate(texts):
    axs[0].text(x[i] + 0.1, y[i], txt, fontsize=6, rotation=10)

## Red for significant blue hue for rest
colors_below_p = 'red'
colors_above_p = 'Blues'

from matplotlib.colors import ListedColormap
cmap_above_p = plt.cm.get_cmap(colors_above_p)
colors = [*cmap_above_p(np.linspace(1, 0, num=20))]
colors[0] = colors_below_p
cmap_new = ListedColormap(colors)

heatmap = axs[1].imshow(data_matrix, cmap=cmap_new, interpolation='nearest')
axs[1].set_title('Heatmap')
axs[1].set_xticks(np.arange(data_matrix.shape[1]))
axs[1].set_xticklabels((txts.replace("COVID", "") for txts in ["COVIDsevere", "COVIDmoderate", "COVIDmild", "postCOVID", "COVIDasymptomatic"]), rotation=45, ha='right')
axs[1].set_yticks(np.arange(data_matrix.shape[0]))
axs[1].set_yticklabels(sorted(significance_dict.keys()))
plt.colorbar(heatmap, ax=axs[1], label='p-values')
   	
axs[0].set_title('Significance of expression\nfold change differences')
axs[1].set_title("Significance for association of\nvariant presence with disease condition")

plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('volcano_plot.png', dpi=300, bbox_inches='tight')