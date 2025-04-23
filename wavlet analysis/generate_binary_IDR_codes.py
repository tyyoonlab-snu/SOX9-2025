"""
Author: Tae Gyun Kim
This script reads an Excel file containing intrinsically disordered region (IDR) sequences and associated metadata, 
and generates binary-encoded representations of these sequences based on aromatic residues (Y, F, W).

The results are saved in an Excel file.
"""

#%% Imports
import pandas as pd

#%% Read input Excel file
input_path = "example_input_path.xlsx"  # Example path to input Excel file

data = pd.read_excel(input_path)

# Extract relevant columns
uniprot_ids = data['UniprotID'].tolist()
gene_names = data['Name'].tolist()
idrs = data['IDR'].tolist()
differences = data['Difference'].tolist()
pvalues = data['-log10 pvalue'].tolist()
fold_changes = data['log2 fold change'].tolist()

#%% Define aromatic amino acids and non-aromatics
aromatics = ['Y', 'F', 'W']
aa_all = list('ACDEFGHIKLMNPQRSTVWY')
non_aromatics = [aa for aa in aa_all if aa not in aromatics]

#%% Generate binary-coded data for aromatics
ones_map = str.maketrans(''.join(aromatics), '1' * len(aromatics))
zeros_map = str.maketrans(''.join(non_aromatics), '0' * len(non_aromatics))

binary_df = pd.DataFrame(columns=[
    'UniprotID', 'Name', 'IDR', 'binary code', 'Difference',
    '-log10 pvalue', 'log2 fold change']
)

for i in range(len(uniprot_ids)):
    idr_seq = idrs[i]
    temp_code = idr_seq.translate(zeros_map).translate(ones_map)

    entry = {
        'UniprotID': uniprot_ids[i],
        'Name': gene_names[i],
        'IDR': idrs[i],
        'binary code': temp_code,
        'Difference': differences[i],
        '-log10 pvalue': pvalues[i],
        'log2 fold change': fold_changes[i]
    }
    binary_df = pd.concat([binary_df, pd.DataFrame([entry])], ignore_index=True)

print(f"{i+1} proteins processed for aromatic residues")

#%% Save output to Excel
output_path = input_path.replace('.xlsx', '_aromatic_binary_code.xlsx')
binary_df.to_excel(output_path, sheet_name='Aromatic', index=False)
