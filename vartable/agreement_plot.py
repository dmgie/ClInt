import os
# from turtle import color
import pandas as pd
import matplotlib.pyplot as plt
import pickle
from adjustText import adjust_text
import numpy as np


data_path = './vartable_output_grch37'
data = []
all_data = {}

def flatten(xss):
    return [x for xs in xss for x in xs]

def extract_matching_count(file):
    with open(file, 'r') as file:
        lines = file.readlines()

    desired_line = None
    for line in lines:
        if 'MATCHING' in line:
            desired_line = line
            break

    matching_count = None
    if desired_line:
        matching_count = int(desired_line.split('MATCHING: ')[1].split()[0])

    return matching_count

def extract_agreement(file, patient_number):
    df = pd.read_csv(file, delimiter='\t', skiprows=1)
    agreement_dict_list = []

    for index, row in df.iterrows():
            agreement_dict = {}
            agreement_rate_dna = row['agreement_rate_dna']
            agreement_rate_rna = row['agreement_rate_rna']
            agreement_dict[patient_number+str(index)] = {'DNA-RNA': float(agreement_rate_dna), 
                                                        'RNA-DNA': float(agreement_rate_rna), 
                                                        "COUNT_MATCH": extract_matching_count(file)}
            agreement_dict_list.append(agreement_dict)
    
    return agreement_dict_list


for patient_folder in os.listdir(data_path):
    
    patient_files = os.path.join(data_path, patient_folder)
    if os.path.isdir(patient_files):
        patient_number = patient_folder.split('_')[-1]

        for file in os.listdir(patient_files):
            if 'agreement' in file:
                condition = file.split('_')[0]
                file_path = os.path.join(patient_files, file)
                agreement_data = extract_agreement(file_path, patient_number+"_"+condition.replace("COVID", ""))
                data.append(agreement_data)


data = flatten(data)

for d in data:
    all_data.update(d)

df = pd.DataFrame(all_data).transpose()


runtimes = {}

with open('./vartable_output_grch37/runtimes.pkl', 'rb') as handle:
    runtimes = pickle.load(handle)

rt_counts_dict = {}

for pat in runtimes:
    pat_trimmed = pat.replace("COVID", "").replace(" ", "")+"0"
    if pat_trimmed in df.index:
        rt_counts_dict[pat_trimmed] = {"rt": runtimes[pat],
                                       "counts": df.loc[pat_trimmed, "COUNT_MATCH"]}
    else:
        print(pat_trimmed, "NOT FOUND")
    


rt_counts_df = pd.DataFrame(rt_counts_dict).T


fig, axs = plt.subplots(1,2, figsize=(10, 6))

axs[1].set_title('DNA-RNA vs. RNA-DNA agreement rates')
axs[1].scatter(df['DNA-RNA'], df['RNA-DNA'])
axs[1].grid(True)
axs[1].set_xlabel('DNA-RNA agreement rate')
axs[1].set_ylabel('RNA-DNA agreement rate')


# for i, patient in enumerate(data):
#     for txt in patient:
#         axs[1].annotate(txt, (df['DNA-RNA'].loc[txt], df['RNA-DNA'].loc[txt]), xytext=(5, -5), textcoords='offset points', fontsize=3)


# texts = [axs[1].text(df['DNA-RNA'].iloc(i), df['RNA-DNA'].iloc(i), 'Text%s' %i) for i in range(len(df['DNA-RNA']))]
# adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'))




axs[0].set_title('Matching counts vs. runtime')
axs[0].scatter(rt_counts_df['counts'], rt_counts_df['rt'])
axs[0].grid(True)
axs[0].set_xlabel('Number of DNA/RNA matching variants')
axs[0].set_ylabel('VarTable runtime in s')

y_rts = rt_counts_df['rt'].values
x_cts = rt_counts_df['counts'].values

coefficients = np.polyfit(x_cts, y_rts, 1)
poly = np.poly1d(coefficients)

print(coefficients)
axs[0].plot(x_cts, poly(x_cts), color="grey")

plt.savefig('agreement_rates_without_labels.png', dpi=300, bbox_inches='tight')
