import pandas as pd

dir = 'public-data/preterm/metadata/'

# Read the two CSV files into DataFrames
df_samples = pd.read_csv(f'{dir}infant-samples.csv')
print(df_samples)
df_timepoints = pd.read_csv(f'{dir}time-points.csv')

df_sra = pd.read_csv(f'{dir}SraRunInfo.csv')
df_run_samples = pd.DataFrame()
df_run_samples['Run'] = df_sra['Run']
df_run_samples['Sample'] = df_sra['LibraryName']
df_run_samples['Size_Mb'] = df_sra['size_MB']
print(df_run_samples)

# Perform a merge (join) based on the 'Sample' column
df_merged = pd.merge(df_samples, df_timepoints, on='Sample', how='inner')

df_final = pd.merge(df_merged, df_run_samples, on='Sample', how='inner')

# Save the merged DataFrame to a new CSV file if needed
df_final.to_csv(f'{dir}sample-metadata.csv', index=False)

df_infant = df_final[~df_final['Infant'].duplicated(keep='first')]
df_infant.to_csv(f'{dir}unique-infants.csv')


# Infants who had antibiotics before birth
df_infant = df_infant[df_infant['ABCourses'] != 'OnlyEarlyAB']
df_infant = df_infant[df_infant['T'] == 'T1']
df_infant.to_csv(f'{dir}useable-patients.csv')
print(df_infant['Infant'])

print(df_final)

infants = [
    'P320', 'P284', 'P358', 'P346',
    'P072', 'P128', 'P184', 'P533',
    'P375', 'P420', 'P241', 'P209',
    'P076', 'P096', 'P094', 'P664',
    'P721', 'P691', 'P414', 'P541',
    'P277', 'P379', 'P399'
]



# Get all the SRR run accessions:
df_final = df_final[df_final['Infant'].isin(infants)]
df_runs = pd.DataFrame()
df_runs['Infant'] = df_final['Infant']
df_runs['Run'] = df_final['Run']
df_runs.to_csv(f'{dir}sra-accs.csv')