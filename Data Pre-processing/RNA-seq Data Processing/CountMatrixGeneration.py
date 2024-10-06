import pandas as pd
import os

# Counts files from HT-seq
archivos = ['C1.txt', 'C2.txt', 'C3.txt', 'C4.txt', 'C5.txt', 'C6.txt', 'C7.txt', 'C8.txt', 'T1.txt', 'T2.txt', 'T3.txt', 'T4.txt', 'T5.txt', 'T6.txt', 'T7.txt', 'T8.txt']

# Column delimiter
delimitador = '\t'  # Adjust to own data

# Final Count Matriz dataframe
df_final = pd.DataFrame()

for i, archivo in enumerate(archivos):
    # Read as dataframe
    df = pd.read_csv(archivo, sep=delimitador, header=None, names=['col1', os.path.splitext(archivo)[0]])
    
    # Copy first file both columns
    if i == 0:
        df_final = df.copy()
    else:
        # Only second column from other files
        # Merge with final dataframe to align first column terms
        df_final = pd.merge(df_final, df, on='col1', how='outer')
        
        # Replace NaN (Not apearing terms) with 0
        df_final.fillna(0, inplace=True)

# Save
df_final.to_csv('archivo_final.tsv', sep='\t', index=False)
