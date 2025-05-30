# -*- coding: utf-8 -*-
"""CCR10 Expression in Monocytes - Scanpy analysis Zenodo.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1E4vOgshTxu7Fg4UdsP7mdbICqTaTkPLb
"""

from google.colab import drive
drive.mount('/content/drive')

!pip install scanpy
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.io

file_path = "/content/drive/My Drive/scRNAseq data/panImmune.h5ad"

adata = sc.read_h5ad(file_path)

print(adata)

adata.var_names[:5]  # Show first 5 gene names

adata.obs_names[:5]  # Show first 5 cell IDs

print(adata.shape)  # (cells, genes)

# Check unique cell types in the 'lineage' column
print(adata.obs['lineage'].unique())

unique_lineages = adata.obs['lineage'].unique().tolist()
print(unique_lineages)

for lineage in adata.obs['lineage'].unique():
    print(lineage)

# Check unique cell types in the 'sample_id' column
print(adata.obs['sample_id'].unique())

# Check unique cell types in the 'tissue' column
print(adata.obs['tissue'].unique())

import re

# Define the patterns for each category
mm_patterns = r'Foster_2024.MM4_neg|Foster_2024.MM4_pos|Foster_2024.MM5_pos|Foster_2024.MM6_pos|Foster_2024.MM6_neg|' \
              r'Foster_2024.MM5_neg|Foster_2024.MM3_pos|Foster_2024.MM3_neg|Foster_2024.MM7_neg|Foster_2024.MM7_pos|' \
              r'Bailur_2019.mm2|Bailur_2019.mm5|Bailur_2019.mm4|Bailur_2019.mm3|Bailur_2019.mm11|Bailur_2019.mm10|' \
              r'Bailur_2019.mm8|Bailur_2019.mm6|Bailur_2019.mm1|Bailur_2019.mm7|Bailur_2019.mm9|' \
              r'Zavidij_2020.MM-5|Zavidij_2020.MM-1|Zavidij_2020.MM-8|Zavidij_2020.MM-6|Zavidij_2020.MM-3|' \
              r'Zavidij_2020.MM-7|Zavidij_2020.MM-4|Zheng_2021.P20181219-T|Zheng_2021.P20190122-P|' \
              r'Zheng_2021.P20190122-T|Zheng_2021.P20190322-P|Zheng_2021.P20190322-T|' \
              r'Liu_2021.77570_Primary|Liu_2021.59114_Primary|Liu_2021.83942_Primary|Liu_2021.58408_Primary|' \
              r'Liu_2021.47491_Primary|Liu_2021.57075_Primary|Liu_2021.37692_Primary|Liu_2021.60359_Primary|' \
              r'Liu_2021.56203_Primary|Maura_2023.PT11|Maura_2023.PT12|Maura_2023.PT18|Maura_2023.PT13|' \
              r'Maura_2023.PT15|Maura_2023.PT55|Maura_2023.PT78|Maura_2023.PT85|Maura_2023.PT03|' \
              r'Maura_2023.PT30|Maura_2023.PT32|Maura_2023.PT33|Maura_2023.PT39|Maura_2023.PT49|' \
              r'Maura_2023.PT58|Maura_2023.PT59|Maura_2023.PT63|Foster_2024.MM2|Foster_2024.MM1'
smm_patterns = r'Foster_2024.SMM9_neg|Foster_2024.SMM7_pos|Foster_2024.SMM3_neg|Foster_2024.SMM7_neg|' \
              r'Foster_2024.SMM9_pos|Foster_2024.SMM12_pos|Foster_2024.SMM6_neg|Foster_2024.SMM8_pos|' \
              r'Foster_2024.SMM8_neg|Foster_2024.SMM6_pos|Foster_2024.SMM12_neg|Foster_2024.SMM4_pos|' \
              r'Foster_2024.SMM5_pos|Foster_2024.SMM5_neg|Foster_2024.SMM4_neg|Foster_2024.SMM11_neg|' \
              r'Foster_2024.SMM10_pos|Foster_2024.SMM2_neg|Foster_2024.SMM1_pos|Foster_2024.SMM1_neg|' \
              r'Foster_2024.SMM2_pos|Foster_2024.SMM10_neg|Foster_2024.SMM11_pos|Foster_2024.SMM14_neg|' \
              r'Foster_2024.SMM14_pos|Foster_2024.SMM13_neg|Foster_2024.SMM13_pos|' \
              r'Zavidij_2020.SMMh-4|Zavidij_2020.SMMl-1|Zavidij_2020.SMMh-10|Zavidij_2020.SMMh-3|' \
              r'Zavidij_2020.SMMh-7|Zavidij_2020.SMMh-2|Zavidij_2020.SMMh-9|Zavidij_2020.SMMh-6|' \
              r'Zavidij_2020.SMMl-3|Zavidij_2020.SMMh-8|Zavidij_2020.SMMl-2|' \
              r'Liu_2021.47491_SMM|Liu_2021.58408_SMM|Liu_2021.37692_SMM'
mgus_patterns = r'Foster_2024.MGUS1_pos|Foster_2024.MGUS1_neg|' \
                r'Bailur_2019.mgus14|Bailur_2019.mgus13|Bailur_2019.mgus4|Bailur_2019.mgus3|' \
                r'Bailur_2019.mgus2|Bailur_2019.mgus5|Bailur_2019.mgus12|Bailur_2019.mgus10|' \
                r'Bailur_2019.mgus9|Bailur_2019.mgus7|Bailur_2019.mgus6|Bailur_2019.mgus1|' \
                r'Bailur_2019.mgus8|Bailur_2019.mgus11|' \
                r'Zavidij_2020.MGUS-1|Zavidij_2020.MGUS-4|Zavidij_2020.MGUS-3|' \
                r'Zavidij_2020.MGUS-2|Zavidij_2020.MGUS-6'
hd_patterns = r'Foster_2024.HD3_neg|Foster_2024.HD2_pos|Foster_2024.HD2_neg|' \
              r'Foster_2024.HD3_pos|Foster_2024.HD1_pos|Foster_2024.HD1_neg|' \
              r'Oetjen_2018.G|Oetjen_2018.S_T1|Oetjen_2018.U|Oetjen_2018.T|' \
              r'Oetjen_2018.W|Oetjen_2018.L|Oetjen_2018.B|Oetjen_2018.C_T1_S1|' \
              r'Oetjen_2018.O|Oetjen_2018.P|Oetjen_2018.M|Oetjen_2018.N|' \
              r'Oetjen_2018.Q|Oetjen_2018.J|Oetjen_2018.C_T1_S2|Oetjen_2018.A|' \
              r'Oetjen_2018.F|Oetjen_2018.C_T2|Oetjen_2018.S_T2|Oetjen_2018.K|' \
              r'Oetjen_2018.E|Oetjen_2018.H|Oetjen_2018.R|' \
              r'Bailur_2019.hd1|Bailur_2019.hd6|Bailur_2019.hd8|Bailur_2019.hd7|' \
              r'Bailur_2019.hd5|Bailur_2019.hd2|Bailur_2019.hd3|Bailur_2019.hd4|' \
              r'Zavidij_2020.NBM-10|Zavidij_2020.NBM-5|Zavidij_2020.NBM-9|' \
              r'Zavidij_2020.NBM-11|Zavidij_2020.NBM-2|Zavidij_2020.NBM-1|' \
              r'Zavidij_2020.NBM-8|Kfoury_2021.BMM4|Kfoury_2021.BMM8|Kfoury_2021.BMM3|' \
              r'Kfoury_2021.BMM9|Kfoury_2021.BMM5|Kfoury_2021.BMM6|Kfoury_2021.BMM2|' \
              r'Granja_2019.BMMC_D1T1|Granja_2019.BMMC_D1T2|Conde_2022.A29_BMA|' \
              r'Conde_2022.A31_BMA|Conde_2022.A36_BMA|Conde_2022.A36_BLD|' \
              r'Conde_2022.A35_BMA|Conde_2022.A35_BLD|Conde_2022.A37_BMA|' \
              r'Conde_2022.621B_BMA|Conde_2022.637C_BMA|Conde_2022.637C_BLD|' \
              r'Conde_2022.640C_BMA|Conde_2022.D503_BMA|Conde_2022.D503_BLD|' \
              r'Conde_2022.D496_BMA|Conde_2022.D496_BLD|Stephenson_2021.newcastle65|' \
              r'Stephenson_2021.MH8919226|Stephenson_2021.MH8919333|Stephenson_2021.MH8919332|' \
              r'Stephenson_2021.MH8919227|Stephenson_2021.MH8919283|Stephenson_2021.MH8919178|' \
              r'Stephenson_2021.MH8919177|Stephenson_2021.MH8919176|Stephenson_2021.MH8919179|' \
              r'Stephenson_2021.newcastle74|Stephenson_2021.MH8919282|Stephenson_2021.CV0904|' \
              r'Stephenson_2021.CV0902|Stephenson_2021.CV0911|Stephenson_2021.CV0929|' \
              r'Stephenson_2021.CV0915|Stephenson_2021.CV0917|Stephenson_2021.CV0939|' \
              r'Stephenson_2021.CV0926|Stephenson_2021.CV0934|Stephenson_2021.CV0940|' \
              r'Stephenson_2021.CV0944'

# Function to categorize sample_id
def categorize_sample(sample):
    if re.search(mm_patterns, sample, re.IGNORECASE):
        return "MM"
    elif re.search(smm_patterns, sample, re.IGNORECASE):
        return "SMM"
    elif re.search(mgus_patterns, sample, re.IGNORECASE):
        return "MGUS"
    elif re.search(hd_patterns, sample, re.IGNORECASE):
        return "HD"
    else:
        return "Unknown"  # Default to Unknown if no match is found

# Extract sample_ids from adata.obs
sample_ids = adata.obs['sample_id'].unique()

# Create the DataFrame
df = pd.DataFrame(sample_ids, columns=["sample_id"])

# Apply the function to categorize the sample_id
df["sample_category"] = df["sample_id"].apply(categorize_sample)

# Display the DataFrame
print(df[['sample_id', 'sample_category']])

# Apply the function to categorize the sample_id and update adata.obs
adata.obs['sample_category'] = adata.obs['sample_id'].apply(categorize_sample)

# Check the result
print(adata.obs[['sample_id', 'sample_category']])

import pandas as pd

# Example: Create a DataFrame for sample_id and sample_category comparison
sample_ids = adata.obs['sample_id'].unique()
df = pd.DataFrame(sample_ids, columns=["sample_id"])

# Apply the categorization function to the sample_id column
df["sample_category"] = df["sample_id"].apply(categorize_sample)

# Save the DataFrame to a CSV file
df.to_csv('/content/sample_comparison.csv', index=False)

# Provide a download link
from google.colab import files
files.download('/content/sample_comparison.csv')

# Replace sample_id with sample_category
adata.obs['sample_id'] = adata.obs['sample_category']

# Drop the sample_category column if no longer needed
adata.obs.drop('sample_category', axis=1, inplace=True)

# Check the result
print(adata.obs[['sample_id']])

# Filter adata for Myeloid cells and BM tissue
My_bm_adata = adata[(adata.obs['lineage'] == 'Myeloid') & (adata.obs['tissue'] == 'BM')].copy()

# Check filtered data
print(My_bm_adata)
print(My_bm_adata.obs['sample_id'].value_counts())  # Verify category counts

print(My_bm_adata.obs['scr_doublet'].dtype)
print(My_bm_adata.obs['scr_doublet'].unique())

print("Total cells before filtering:", My_bm_adata.shape[0])

print("Passing n_genes > 500:", My_bm_adata[My_bm_adata.obs['n_genes'] > 500].shape[0])

print("Passing n_genes <= 50000:", My_bm_adata[My_bm_adata.obs['n_genes'] <= 5000].shape[0])

print("Passing n_counts >= 1000:", My_bm_adata[My_bm_adata.obs['n_counts'] >= 1000].shape[0])

print("Passing n_counts <= 10000:", My_bm_adata[My_bm_adata.obs['n_counts'] <= 10000].shape[0])

print("Passing pct_counts_mt < 10:", My_bm_adata[My_bm_adata.obs['pct_counts_mt'] < 10].shape[0])

print("Passing scr_doublet == True:", My_bm_adata[My_bm_adata.obs['scr_doublet'] == True].shape[0])

# Apply QC filters
filtered_adata = My_bm_adata[
    (My_bm_adata.obs['n_genes'] > 500) &
    (My_bm_adata.obs['n_genes'] <= 5000) &
    (My_bm_adata.obs['n_counts'] >= 1000) &
    (My_bm_adata.obs['n_counts'] <= 10000) &
    (My_bm_adata.obs['pct_counts_mt'] < 10)  # Mitochondrial percentage filter
].copy()

# Check filtered data
print(filtered_adata)
print(filtered_adata.obs['sample_id'].value_counts())

### NEW ANALYSIS IN Myeloid CELLS

hd_mm_adata = filtered_adata[filtered_adata.obs['sample_id'].isin(['HD', 'MGUS', 'SMM', 'MM'])].copy()

import scanpy as sc
import matplotlib.pyplot as plt

# Make sure CCR10 is in the data
# Add the CCR10 expression to the plotting
sc.pl.umap(hd_mm_adata, color='CCR10', title='UMAP of CCR10 Expression', show=True)

import scanpy as sc
import matplotlib.pyplot as plt


# Define conditions to subset by
conditions = ['MM', 'MGUS', 'SMM', 'HD']

# Loop through each condition to create a UMAP for that condition
for condition in conditions:
    # Subset the AnnData object to the current condition
    condition_adata = hd_mm_adata[hd_mm_adata.obs['sample_id'] == condition]

    # Plot UMAP for the current condition, specifying color for CCR10 from the gene expression matrix
    # Use condition_adata.raw[:, 'CCR10'].X to access gene expression if you had previously normalized the data
    sc.pl.umap(condition_adata, color=['CCR10'], title=f'UMAP of CCR10 Expression for {condition}', show=True, use_raw=False, vmin=0, # Set minimum expression value
    vmax=5)   # Set maximum expression value # Edited line. Set 'use_raw' to False and provide 'color' as a list

import seaborn as sns
import matplotlib.pyplot as plt

# Fetch the CCR10 expression values from the .X matrix using the gene name
CCR10_expression = hd_mm_adata[:, 'CCR10'].X.toarray().flatten()

# Add CCR10 expression to the hd_mm_adata.obs DataFrame
hd_mm_adata.obs['CCR10'] = CCR10_expression

# Create a boxplot showing CCR10 expression across the groups
plt.figure(figsize=(8, 6))
sns.boxplot(x='sample_id', y='CCR10', data=hd_mm_adata.obs)  # Now 'CCR10' is in hd_mm_adata.obs
plt.title('CCR10 Expression Across Groups (HD, MGUS, SMM, MM)')
plt.ylabel('CCR10 Expression')
plt.xlabel('Group')
plt.show()

