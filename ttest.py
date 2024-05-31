import pandas as pd
from scipy import stats
import os
import numpy as np

control_dir = '/home/jasper/Thesis/CONTROL/output'
exposed_dir = '/home/jasper/Thesis/PET-EXPOSED/output'

control_data = {}
exposed_data = {}
significant_plastics = []

def read_data(directory, data):
    for file in os.listdir(directory):
        if file.endswith('.tsv'):
            df = pd.read_csv(os.path.join(directory, file), sep='\t')
            for i, row in df.iterrows():
                plastic = row['plastic name']
                rpkm = row['rpkm']
                if plastic not in data:
                    data[plastic] = []
                data[plastic].append(rpkm)


read_data(control_dir, control_data)
read_data(exposed_dir, exposed_data)


for plastic in control_data.keys():
    if plastic in exposed_data:
        control_rpkm = ' '.join([f'{val:.2f}' for val in control_data[plastic]])
        exposed_rpkm = ' '.join([f'{val:.2f}' for val in exposed_data[plastic]])
        print(f"{plastic} \n {control_rpkm} \n {exposed_rpkm}")


# Calculate the Bonferroni corrected significance level
num_tests = len(control_data.keys())
bonferroni_corrected_alpha = 0.05 / num_tests

print(f"Bonferroni corrected significance level: {bonferroni_corrected_alpha}")

# Perform a t-test on each plastic's rpkm values
for plastic in control_data.keys():
    if plastic in exposed_data:
        control_rpkm = control_data[plastic]
        exposed_rpkm = exposed_data[plastic]
        if not any(pd.isnull(control_rpkm)) and not any(pd.isnull(exposed_rpkm)):
            if len(set(control_rpkm)) > 1 and len(set(exposed_rpkm)) > 1:
                t_statistic, p_value = stats.ttest_ind(control_rpkm, exposed_rpkm)
                control_avg = np.mean(control_rpkm)
                exposed_avg = np.mean(exposed_rpkm)
                print(f"{plastic}: t-statistic = {t_statistic}, p-value = {p_value}, control average = {control_avg}, exposed average = {exposed_avg}")
                if p_value < bonferroni_corrected_alpha:
                    significant_plastics.append((plastic, p_value, control_avg, exposed_avg))
            else:
                print(f"{plastic}: All values are identical, t-test not performed.")
        else:
            print(f"{plastic}: Missing values detected, t-test not performed.")

for plastic, p_value, control_avg, exposed_avg in significant_plastics:
    print(f"{plastic}: p-value = {p_value}, control average = {control_avg}, exposed average = {exposed_avg}")

