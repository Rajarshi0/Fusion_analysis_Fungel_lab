import pandas as pd
import numpy as np
import os
from pathlib import Path

# ================================================================= #
# === CONFIGURATION: Just change the paths in this section === #
# ================================================================= #

# 1. Set the directory where your sample folders are.
SEARCH_DIRECTORY = "arriba" 

# 2. Set the directory where you want to save the converted files.
OUTPUT_DIRECTORY = "cicero_converted_outputs"

# ================================================================= #
# === You do not need to change anything below this line === #
# ================================================================= #

def convert_single_file(arriba_file_path, output_file_path):
    """
    Converts a single Arriba fusions file to the specified Cicero format.
    """
    try:
        # The exact Cicero header.
        final_cicero_header = [
            'sample', 'geneA', 'chrA', 'posA', 'ortA', 'featureA', 'geneB', 'chrB',
            'posB', 'ortB', 'featureB', 'sv_ort', 'readsA', 'readsB', 'matchA', 'matchB',
            'repeatA', 'repeatB', 'coverageA', 'coverageB', 'ratioA', 'ratioB',
            'qposA', 'qposB', 'total_readsA', 'total_readsB', 'contig', 'type', 'score',
            'rating', 'medal', 'functional effect', 'frame', 'sv_refseqA',
            'sv_refseqA_codon', 'sv_refseqA_exon', 'sv_refseqA_anchor_type',
            'sv_refseqA_coding_base_number', 'sv_refseqA_last_coding_base_number',
            'sv_refseqA_AA_index', 'sv_refseqA_contig_index', 'sv_refseqB',
            'sv_refseqB_codon', 'sv_refseqB_exon', 'sv_refseqB_anchor_type',
            'sv_refseqB_coding_base_number', 'sv_refseqB_last_coding_base_number',
            'sv_refseqB_AA_index', 'sv_refseqB_contig_index', 'sv_AA', 'sv_desc',
            'sv_processing_exception', 'sv_general_info', 'sv_interstitial_AA',
            'sv_frame_index'
        ]

        # Read the Arriba file
        arriba_df = pd.read_csv(arriba_file_path, sep='\t')

        # --- Aggregation to group fusions by gene pair ---
        aggregation_rules = {
            'breakpoint1': 'first', 'breakpoint2': 'first',
            'strand1(gene/fusion)': 'first', 'strand2(gene/fusion)': 'first',
            'site1': 'first', 'site2': 'first', 'type': 'first',
            'split_reads1': 'sum', 'split_reads2': 'sum',
            'coverage1': 'first', 'coverage2': 'first',
            'peptide_sequence': lambda x: ','.join(x.dropna().unique()) if not x.dropna().empty else '.',
            'fusion_transcript': lambda x: ','.join(x.dropna().unique()),
            'reading_frame': 'first',
            'confidence': lambda x: x.mode()[0] if not x.empty else 'low'
        }
        grouped_df = arriba_df.groupby(['#gene1', 'gene2'], as_index=False).agg(aggregation_rules)

        # Create the final DataFrame, starting empty
        final_df = pd.DataFrame()

        # --- Perform Only Direct Mapping and Logical Inference from Arriba to Cicero format ---
        final_df['sample'] = Path(arriba_file_path).parent.name
        final_df['geneA'] = grouped_df['#gene1']
        final_df['chrA'] = grouped_df['breakpoint1'].str.split(':').str[0]
        final_df['posA'] = pd.to_numeric(grouped_df['breakpoint1'].str.split(':').str[1])
        final_df['ortA'] = grouped_df['strand1(gene/fusion)'].str.split('/').str[0]
        final_df['featureA'] = grouped_df['site1']
        final_df['geneB'] = grouped_df['gene2']
        final_df['chrB'] = grouped_df['breakpoint2'].str.split(':').str[0]
        final_df['posB'] = pd.to_numeric(grouped_df['breakpoint2'].str.split(':').str[1])
        final_df['ortB'] = grouped_df['strand2(gene/fusion)'].str.split('/').str[0]
        final_df['featureB'] = grouped_df['site2']

        final_df['readsA'] = grouped_df['split_reads1']
        final_df['readsB'] = grouped_df['split_reads2']
        
        final_df['total_readsA'] = grouped_df['coverage1']
        final_df['total_readsB'] = grouped_df['coverage2']
        final_df['coverageA'] = grouped_df['coverage1']
        final_df['coverageB'] = grouped_df['coverage2']

        total_reads = final_df['readsA'] + final_df['readsB']
        final_df['ratioA'] = np.divide(final_df['readsA'], total_reads, where=total_reads!=0, out=np.zeros_like(final_df['readsA'], dtype=float)).round(4)
        final_df['ratioB'] = np.divide(final_df['readsB'], total_reads, where=total_reads!=0, out=np.zeros_like(final_df['readsB'], dtype=float)).round(4)

        final_df['type'] = grouped_df['type']
        final_df['score'] = grouped_df['confidence']
        confidence_mapping = {'high': 'HQ', 'medium': 'MQ', 'low': 'LQ'}
        final_df['rating'] = grouped_df['confidence'].map(confidence_mapping)
        final_df['frame'] = grouped_df['reading_frame']
        final_df['functional effect'] = grouped_df['reading_frame']

        final_df['contig'] = grouped_df['peptide_sequence']
        final_df['sv_AA'] = grouped_df['peptide_sequence']
        final_df['sv_desc'] = grouped_df['fusion_transcript']

        # Ensure final structure has all columns from the Cicero header, filling missing ones
        output_df = pd.DataFrame(final_df, columns=final_cicero_header)

        # Write the final file
        output_df.to_csv(output_file_path, sep='\t', index=False, na_rep='.')
        
        return True

    except Exception as e:
        print(f"    ❌ An error occurred processing {Path(arriba_file_path).name}: {e}")
        return False

# --- Main Batch Processing Logic ---
search_path = Path(SEARCH_DIRECTORY)
output_path = Path(OUTPUT_DIRECTORY)

# Create the output directory if it doesn't exist
output_path.mkdir(parents=True, exist_ok=True)

print(f"🚀 Starting batch conversion...")
print(f"   Searching for 'arriba_fusions.tsv' in: '{search_path.resolve()}'")
print(f"   Saving converted files to: '{output_path.resolve()}'")

# Use glob to find all matching Arriba files recursively
arriba_files = list(search_path.rglob('fusions.tsv'))

if not arriba_files:
    print(f"⚠️ No files named 'arriba_fusions.tsv' were found. Please check the SEARCH_DIRECTORY path.")
else:
    success_count = 0
    for arriba_file in arriba_files:
        sample_name = arriba_file.parent.name
        output_file = output_path / f"{sample_name}_cicero_converted.tsv"
        
        print("--------------------------------")
        print(f"➡️  Processing sample: {sample_name}")
        
        if convert_single_file(arriba_file, output_file):
            print(f"   ✅ Success! Saved to: {output_file.name}")
            success_count += 1
        else:
            print(f"   ❌ Failed to convert.")

    print("--------------------------------")
    print(f"🎉 All done. {success_count} of {len(arriba_files)} files converted successfully.")
