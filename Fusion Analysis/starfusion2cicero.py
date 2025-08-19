import pandas as pd
from pathlib import Path
import sys

# ================================================================= #
# === CONFIGURATION: Just change the paths in this section      === #
# ================================================================= #

# 1. Directory to search for your STAR-Fusion output folders.
SEARCH_DIRECTORY = "/Fungel/rajarshi/02_processed_data_08-1-2025/fusion_analysis_output/star_fusion/"

# 2. Directory where you want to save the final count files.
OUTPUT_DIRECTORY = "SF2cicero_output"

# ================================================================= #
# === You do not need to change anything below this line        === #
# ================================================================= #

def process_starfusion_file(file_path):
    """
    Reads a star-fusion.fusion_predictions.tsv file, extracts key columns,
    and returns a clean pandas DataFrame.
    """
    try:
        # === ROBUST MANUAL PARSING LOGIC ===
        with open(file_path, 'r') as f:
            # Read the header line and split it by tabs to get column names.
            header = f.readline().strip().lstrip('#').split('\t')
            
            # Read all subsequent lines.
            lines = f.readlines()

        # If there are no data lines, the file has no fusions.
        if not lines:
            print(f"    - INFO: No fusions found in {file_path.name}. Skipping.")
            return None

        # Create the DataFrame manually from the remaining lines.
        # Use a list of lists where each inner list is a row.
        data = [line.strip().split('\t') for line in lines]
        df = pd.DataFrame(data, columns=header)

        # Define the columns we absolutely need.
        required_cols = ['FusionName', 'JunctionReadCount', 'SpanningFragCount']
        if not all(col in df.columns for col in required_cols):
            print(f"    - ERROR: File {file_path.name} is missing required columns after manual parsing. Skipping.")
            return None

        # Select only the columns we need.
        subset_df = df[required_cols].copy()

        # Parse GeneA and GeneB from the #FusionName column.
        subset_df[['geneA', 'geneB']] = subset_df['FusionName'].str.split('--', n=1, expand=True)

        return subset_df

    except Exception as e:
        print(f"    - ERROR: Could not process file {file_path.name}. Reason: {e}")
        return None


def main():
    """
    Main function to run the batch processing.
    """
    search_path = Path(SEARCH_DIRECTORY)
    output_path = Path(OUTPUT_DIRECTORY)

    # Create the output directory if it doesn't exist.
    output_path.mkdir(parents=True, exist_ok=True)

    print(f"🚀 Starting Fusion Count Extraction...")
    print(f"   Searching in: '{search_path.resolve()}'")
    print(f"   Saving results to: '{output_path.resolve()}'")
    print("-" * 40)

    # Find all the STAR-Fusion result files.
    starfusion_files = list(search_path.rglob('star-fusion.fusion_predictions.tsv'))

    if not starfusion_files:
        print("⚠️ No 'star-fusion.fusion_predictions.tsv' files found. Please check your SEARCH_DIRECTORY.")
        sys.exit(1)

    success_count = 0
    for file in starfusion_files:
        sample_name = file.parent.name
        print(f"Processing sample: {sample_name}...")

        # Process the file to get the data.
        processed_df = process_starfusion_file(file)

        # If processing was successful and returned data...
        if processed_df is not None:
            # Add the sample name as the first column.
            processed_df.insert(0, 'sample', sample_name)

            # Define the final, exact header and column order.
            final_header = ['sample', 'geneA', 'geneB', 'JunctionReadCount', 'SpanningFragCount']
            final_df = processed_df[final_header]

            # Define the output file path.
            output_file = output_path / f"{sample_name}_starfusion2cicero_converted.tsv"

            # Save the result.
            final_df.to_csv(output_file, sep='\t', index=False)
            print(f"   ✅ Success! Saved to {output_file.name}")
            success_count += 1

    print("-" * 40)
    print(f"🎉 All done. {success_count} of {len(starfusion_files)} files with fusions were processed successfully.")


if __name__ == "__main__":
    main()
