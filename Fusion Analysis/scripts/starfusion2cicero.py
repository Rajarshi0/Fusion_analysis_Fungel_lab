import pandas as pd
from pathlib import Path
import sys

# ================================================================= #
# === CONFIGURATION: Paths corrected based on your file structure === #
# ================================================================= #

# 1. Directory containing all your sample folders (P173PB, etc.).
SEARCH_DIRECTORY = "/Fungel/rajarshi/02_processed_data_08-1-2025/fusion_analysis_output/star_fusion"

# 2. Directory where your final, fully formatted Cicero files will be saved.
OUTPUT_DIRECTORY = "SF2cicero_converted_batch1"

# ================================================================= #
# === You do not need to change anything below this line        === #
# ================================================================= #

def convert_starfusion_to_cicero(file_path):
    """
    Reads a STAR-Fusion file, extracts key info, and formats it directly
    into the full Cicero specification.
    """
    try:
        # === Step 1: Robustly parse the STAR-Fusion file ===
        with open(file_path, 'r') as f:
            header = f.readline().strip().lstrip('#').split('\t')
            lines = f.readlines()

        if not lines:
            print(f"    - INFO: No fusions found in {file_path.name}. Skipping.")
            return None

        data = [line.strip().split('\t') for line in lines]
        df = pd.DataFrame(data, columns=header)
        
        required_cols = ['FusionName', 'JunctionReadCount', 'SpanningFragCount']
        if not all(col in df.columns for col in required_cols):
            print(f"    - ERROR: File {file_path.name} is missing required columns. Skipping.")
            return None

        # === Step 2: Create the full Cicero DataFrame ===
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
        
        cicero_df = pd.DataFrame(columns=final_cicero_header)

        # === Step 3: Map the extracted data into the Cicero format ===
        gene_pairs = df['FusionName'].str.split('--', n=1, expand=True)
        
        cicero_df['sample'] = file_path.parent.name
        cicero_df['geneA'] = gene_pairs[0]
        cicero_df['geneB'] = gene_pairs[1]
        cicero_df['readsA'] = df['JunctionReadCount']
        cicero_df['readsB'] = df['SpanningFragCount']
        
        # Fill all other columns with the placeholder '.'
        cicero_df.fillna('.', inplace=True)
        
        return cicero_df

    except Exception as e:
        print(f"    - ERROR: Could not process file {file_path.name}. Reason: {e}")
        return None

def main():
    """
    Main function to run the batch processing.
    """
    search_path = Path(SEARCH_DIRECTORY)
    output_path = Path(OUTPUT_DIRECTORY)

    output_path.mkdir(parents=True, exist_ok=True)

    print(f"🚀 Starting STAR-Fusion to Cicero Conversion...")
    print(f"   Searching for input files in: '{search_path.resolve()}'")
    print(f"   Saving final Cicero files to: '{output_path.resolve()}'")
    print("-" * 40)

    starfusion_files = list(search_path.rglob('star-fusion.fusion_predictions.tsv'))

    if not starfusion_files:
        print(f"⚠️ No 'star-fusion.fusion_predictions.tsv' files found. Please check your SEARCH_DIRECTORY path.")
        sys.exit(1)

    success_count = 0
    for file in starfusion_files:
        sample_name = file.parent.name
        print(f"Processing sample: {sample_name}...")
        
        final_cicero_df = convert_starfusion_to_cicero(file)
        
        if final_cicero_df is not None:
            output_file = output_path / f"{sample_name}_SF2cicero_converted.tsv"
            final_cicero_df.to_csv(output_file, sep='\t', index=False)
            print(f"   ✅ Success! Saved to {output_file.name}")
            success_count += 1
            
    print("-" * 40)
    print(f"🎉 All done. {success_count} of {len(starfusion_files)} files with fusions were processed and converted successfully.")

if __name__ == "__main__":
    main()
