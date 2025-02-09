#!/usr/bin/env python3
import os
import re
import argparse
import pandas as pd
import logging

def parse_fa(lipid):
    """
    Parse a string like 'FA(5:1)' to extract two integers (e.g., (5, 1)).
    If the string doesn't match the pattern, return (999999, 999999).
    """
    pattern = r"FA\((\d+):(\d+)\)"
    match = re.match(pattern, lipid)
    if match:
        major = int(match.group(1))
        minor = int(match.group(2))
        return (major, minor)
    else:
        return (999999, 999999)

def main():
    parser = argparse.ArgumentParser(
        description="Process a parquet file by filtering OzESI_Intensity, grouping by Lipid and Sample, "
                    "selecting the top rows per group, and sorting by the numeric values in the Lipid column. "
                    "Also saves rows with a minor value of 0 in Lipid to a separate file if requested, both as parquet and CSV."
    )
    parser.add_argument("input_file", help="Path to the input parquet file.")
    parser.add_argument("output_file", help="Path to the output parquet file for non-0 values.")
    parser.add_argument("--how_many", type=int, default=2, help="Number of top rows to keep per group (default: 2).")
    parser.add_argument("--threshold", type=float, default=1000, help="Minimum OzESI_Intensity to include (default: 1000).")
    parser.add_argument("--zero_output_dir", default="", help="Directory to save rows with Lipid minor value 0.")
    
    args = parser.parse_args()

    # Set up logging to show debug messages on the console.
    logging.basicConfig(level=logging.DEBUG, format="%(levelname)s: %(message)s")

    # Read the input parquet file
    df = pd.read_parquet(args.input_file)
    logging.debug("Input data contains %d rows", len(df))
    
    # Log the unique Lipid values before filtering (if desired)
    unique_lipids = df["Lipid"].unique()
    logging.debug("Unique Lipid values before filtering: %s", unique_lipids)
    
    # Filter rows based on the threshold for OzESI_Intensity
    df_filtered = df[df["OzESI_Intensity"] >= args.threshold]
    logging.debug("After filtering, %d rows remain (threshold: %s)", len(df_filtered), args.threshold)
    
    # *** STEP: Separate rows with Lipid minor value == 0 ***
    # Create a temporary parsed column
    df_filtered["fa_parsed"] = df_filtered["Lipid"].apply(parse_fa)
    # Separate the rows based on the minor value.
    df_zero = df_filtered[df_filtered["fa_parsed"].apply(lambda t: t[1] == 0)]
    df_nonzero = df_filtered[df_filtered["fa_parsed"].apply(lambda t: t[1] != 0)]
    logging.debug("Number of rows with Lipid minor value 0: %d", len(df_zero))
    logging.debug("Number of rows with non-zero minor value: %d", len(df_nonzero))
    
    # Remove the temporary 'fa_parsed' column from non-zero rows
    df_nonzero = df_nonzero.drop(columns=["fa_parsed"])
    
    # Process non-zero rows: Group by Lipid and Sample and select top rows.
    groups_original = df.groupby(["Lipid", "Sample"]).ngroups
    groups_filtered = df_nonzero.groupby(["Lipid", "Sample"]).ngroups
    logging.debug("Number of groups before filtering: %d", groups_original)
    logging.debug("Number of groups after filtering (non-zero): %d", groups_filtered)
    
    grouped = df_nonzero.groupby(["Lipid", "Sample"])
    top_rows_list = []
    for group_key, group_df in grouped:
        lipid, sample = group_key
        logging.debug("Processing group Lipid=%s, Sample=%s with %d rows", lipid, sample, len(group_df))
        if len(group_df) < args.how_many:
            logging.debug("Group Lipid=%s, Sample=%s has only %d rows, expected %d.",
                          lipid, sample, len(group_df), args.how_many)
        top_grp = group_df.nlargest(args.how_many, "OzESI_Intensity")
        logging.debug("For group Lipid=%s, Sample=%s, selected row indices: %s",
                      lipid, sample, top_grp.index.tolist())
        top_rows_list.append(top_grp)
    
    if top_rows_list:
        df_top = pd.concat(top_rows_list)
    else:
        df_top = pd.DataFrame()
        logging.debug("No groups remaining after filtering.")
    
    # Sort the non-zero DataFrame by the numeric values in the Lipid column.
    df_top["fa_sort_tuple"] = df_top["Lipid"].apply(parse_fa)
    df_top = df_top.sort_values(by=["fa_sort_tuple"]).drop(columns=["fa_sort_tuple"])
    
    # Log the unique Lipid values in the output to see which lipids remain.
    unique_output_lipids = df_top["Lipid"].unique()
    logging.debug("Unique Lipid values in non-zero output: %s", unique_output_lipids)
    
    # Ensure the output directory for non-zero rows exists.
    output_dir = os.path.dirname(args.output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    # Save the resulting non-zero DataFrame to a new parquet file.
    df_top.to_parquet(args.output_file, index=False)
    logging.debug("Non-zero output parquet file written to %s", args.output_file)
    
    # *** NEW STEP: Also save non-zero rows as a CSV file in a subdirectory "csv_dir" ***
    nonzero_csv_dir = os.path.join(output_dir, "csv_dir")
    os.makedirs(nonzero_csv_dir, exist_ok=True)
    nonzero_csv_file = os.path.join(nonzero_csv_dir, os.path.splitext(os.path.basename(args.output_file))[0] + ".csv")
    df_top.to_csv(nonzero_csv_file, index=False)
    logging.debug("Non-zero output CSV file written to %s", nonzero_csv_file)
    
    # *** NEW STEP: Process and save the zero minor value rows if a zero_output_dir is provided ***
    if args.zero_output_dir:
        os.makedirs(args.zero_output_dir, exist_ok=True)
        # For the zero values, group by Lipid and select the top 1 row by OzESI_Intensity for each Lipid.
        grouped_zero = df_zero.groupby("Lipid", group_keys=False).apply(lambda grp: grp.nlargest(1, "OzESI_Intensity"))
        # Optionally, sort the zero values as well based on numeric values in the Lipid column.
        grouped_zero["fa_sort_tuple"] = grouped_zero["Lipid"].apply(parse_fa)
        grouped_zero = grouped_zero.sort_values(by=["fa_sort_tuple"]).drop(columns=["fa_sort_tuple"])
        # Remove the temporary 'fa_parsed' column if it exists.
        if "fa_parsed" in grouped_zero.columns:
            grouped_zero = grouped_zero.drop(columns=["fa_parsed"])
        # Construct the output file path in the zero directory using the same basename.
        zero_output_file = os.path.join(args.zero_output_dir, os.path.basename(args.output_file))
        grouped_zero.to_parquet(zero_output_file, index=False)
        logging.debug("Zero minor value output parquet file written to %s", zero_output_file)
        
        # Also save the zero values as a CSV file in a subdirectory "csv_dir" within the zero_output_dir.
        zero_csv_dir = os.path.join(args.zero_output_dir, "csv_dir")
        os.makedirs(zero_csv_dir, exist_ok=True)
        zero_csv_file = os.path.join(zero_csv_dir, os.path.splitext(os.path.basename(args.output_file))[0] + ".csv")
        grouped_zero.to_csv(zero_csv_file, index=False)
        logging.debug("Zero minor value output CSV file written to %s", zero_csv_file)

if __name__ == "__main__":
    main()
