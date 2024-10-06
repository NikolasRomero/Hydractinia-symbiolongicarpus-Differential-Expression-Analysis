import pandas as pd
import os

## Counts files from HT-seq output reeplace TXT_FILES with your own data
filer = [TXT_FILES]

## Column delimiter
delimiter = '\t'  # Adjust to your own data

## Final Count matrix dataframe
df_final = pd.DataFrame()

for i, archivo in enumerate(files):
    ## Read as a dataframe
    df = pd.read_csv(archivo, sep=delimiter, header=None, names=['col1', os.path.splitext(file)[0]])
    
    ## Copy first file's both columns
    if i == 0:
        df_final = df.copy()
    else:
        ## Only second column from the other files
        ## Merge in final dataframe to align with the first column terms
        df_final = pd.merge(df_final, df, on='col1', how='outer')
        
        ## Replace NaN with 0
        df_final.fillna(0, inplace=True)

## Save
df_final.to_csv('final_file.tsv', sep='\t', index=False)
