#!/usr/bin/env python3
import os
import re
import argparse
import pandas as pd

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
        description="Process a parquet file by filtering OzESI_Intensity, grouping by Lipid, Isomer, Sample, "
                    "selecting the top rows per group, and sorting by the numeric values in the Lipid column."
    )
    parser.add_argument("input_file", help="Path to the input parquet file.")
    parser.add_argument("output_file", help="Path to the output parquet file.")
    parser.add_argument("--how_many", type=int, default=2, help="Number of top rows to keep per group (default: 2).")
    parser.add_argument("--threshold", type=float, default=1000, help="Minimum OzESI_Intensity to include (default: 1000).")
    
    args = parser.parse_args()

    # Read the input parquet file
    df = pd.read_parquet(args.input_file)

    # Filter rows based on the threshold for OzESI_Intensity
    df = df[df["OzESI_Intensity"] >= args.threshold]

    # Group by Lipid, Isomer, and Sample, then select the top 'how_many' rows by OzESI_Intensity in each group
    df_top = (
        df.groupby(["Lipid", "Isomer", "Sample"], group_keys=False)
          .apply(lambda grp: grp.nlargest(args.how_many, "OzESI_Intensity"))
    )

    # Create a temporary column to sort by the numeric values in the Lipid column (e.g., FA(5:0) before FA(10:0))
    df_top["fa_sort_tuple"] = df_top["Lipid"].apply(parse_fa)
    df_top = df_top.sort_values(by=["fa_sort_tuple"]).drop(columns=["fa_sort_tuple"])

    # Ensure the output directory exists
    output_dir = os.path.dirname(args.output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Save the resulting DataFrame to a new parquet file
    df_top.to_parquet(args.output_file, index=False)

if __name__ == "__main__":
    main()
