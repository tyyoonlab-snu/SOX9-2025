"""
Author: TaeGyunKim

This script extracts significant proteins from mass-spectrometry data based on user-defined p-value and fold-change thresholds.
It retrieves matching protein sequences from UniProt, predicts intrinsic disorder regions (IDRs) using MetaPredict,
and compiles the results into a structured Excel file with separate sheets per comparison.
"""

#%% Load reviewed human FASTA sequences from UniProt
import re
import requests
import metapredict as meta
import pandas as pd

# Download UniProt reviewed human proteome in FASTA format
uniprot_url = 'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=((reviewed:true)%20AND%20(organism_id:9606))'
all_fastas = requests.get(uniprot_url).text
fasta_list = re.split(r'\n(?=>)', all_fastas)

#%% Load mass-spec result Excel file (edit the path as necessary)
input_excel = 'C:/temp/file.xlsx' #Path to input Excel file 
sheet_names = ['C9_D_LFQ', 'C9_B8_LFQ', 'B8_D_LFQ']
sheet_data = {name: pd.read_excel(input_excel, sheet_name=name) for name in sheet_names}

#%% Threshold settings
# All p-values and fold-changes are retained; modify below thresholds if filtering is required
pval_threshold = 0   # Corresponds to 1e-0 = 1 (example only — modify as needed)
fc_threshold = 0     # Log2 fold-change cutoff (example only — modify as needed)

# Dictionary to hold output DataFrames per comparison
output_dfs = {}

#%% Main processing loop per sheet
for sheet_name, df in sheet_data.items():
    contrast = sheet_name.replace('_LFQ', '')
    pval_col = f"N: -Log Student's T-test p-value {contrast}"
    fc_col = f"N: Student's T-test Difference {contrast}"

    if pval_col not in df.columns or fc_col not in df.columns:
        print(f"[Warning] Missing columns in sheet {sheet_name}. Skipping...")
        continue

    pval_raw = df[pval_col].tolist()
    fc_raw = df[fc_col].tolist()
    uniprot_ids = df['T: Protein IDs'].tolist()
    gene_names = df['T: Gene names'].tolist()

    # Apply threshold filters
    filt_ids, filt_genes, filt_pvals, filt_fcs = [], [], [], []
    for i in range(len(uniprot_ids)):
        pval = 10 ** (-pval_raw[i])
        fc = fc_raw[i]
        if pval <= 10 ** (-pval_threshold) and abs(fc) >= fc_threshold:
            filt_ids.append(uniprot_ids[i])
            filt_genes.append(gene_names[i])
            filt_pvals.append(pval_raw[i])
            filt_fcs.append(fc)

    # Normalize multiple Uniprot ID entries per row
    ids_clean, genes_clean, pvals_clean, fcs_clean = rearrange_db(filt_pvals, filt_fcs, filt_ids, filt_genes)

    # Retrieve corresponding sequences from UniProt
    final_ids, final_pvals, final_fcs, final_seqs = get_seq_from_uniprot(pvals_clean, fcs_clean, ids_clean)

    # Predict IDRs and record results
    output_rows = []
    for i in range(len(final_ids)):
        if i >= len(genes_clean):
            continue
        idrs = meta.predict_disorder_domains(final_seqs[i]).disordered_domains
        disorder_pct = meta.percent_disorder(final_seqs[i])

        for region in idrs:
            row = {
                'UniprotID': final_ids[i],
                'Name': genes_clean[i],
                'IDR': region,
                '-log10 pvalue': final_pvals[i],
                'log2 fold change': final_fcs[i],
                'Difference': contrast,
                'disorder percent': disorder_pct
            }
            output_rows.append(row)
        print(f"{contrast} - {i+1}/{len(final_ids)} protein processed")

    output_dfs[contrast] = pd.DataFrame(output_rows)

#%% Save final results to Excel with each comparison as a separate sheet
from datetime import datetime
date = datetime.today().strftime('%y%m%d')

outname_full = f"IDR_analysis_MassSpec_MetaPredict_{date}.xlsx"
with pd.ExcelWriter(outname_full, engine='openpyxl') as writer:
    for sheet, df in output_dfs.items():
        df.to_excel(writer, sheet_name=sheet, index=False)

#%% Helper functions
def rearrange_db(pvals, fold_changes, ids, names):
    """Split entries containing multiple IDs into individual rows."""
    split_ids, split_names, split_pvals, split_fcs = [], [], [], []
    for i in range(len(ids)):
        if ';' in ids[i]:
            id_list = ids[i].split(';')
            try:
                name_list = names[i].split(';')
            except AttributeError:
                name_list = [names[i]] * len(id_list)
            if len(name_list) != len(id_list):
                name_list = [names[i]] * len(id_list)

            split_ids.extend(id_list)
            split_names.extend(name_list)
            split_pvals.extend([pvals[i]] * len(id_list))
            split_fcs.extend([fold_changes[i]] * len(id_list))
        else:
            split_ids.append(ids[i])
            split_names.append(names[i])
            split_pvals.append(pvals[i])
            split_fcs.append(fold_changes[i])
    return split_ids, split_names, split_pvals, split_fcs

def get_seq_from_uniprot(pvals, fold_changes, ids):
    """Retrieve protein sequences from UniProt based on Uniprot IDs."""
    sequences, matched_ids, matched_pvals, matched_fcs = [], [], [], []
    for i in range(len(ids)):
        for fasta in fasta_list:
            if ids[i] in fasta:
                seq_start = fasta.find('SV=')
                seq = fasta[seq_start+4:].replace('\n','')
                sequences.append(seq)
                matched_ids.append(ids[i])
                matched_pvals.append(pvals[i])
                matched_fcs.append(fold_changes[i])
                break
    return matched_ids, matched_pvals, matched_fcs, sequences
