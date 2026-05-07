import sys
import argparse
import pandas as pd
import numpy as np
import os

# This script filters tables using pandas query

def main():
    parser = argparse.ArgumentParser(description='Filter table according to an expression (pandas query)')
    parser.add_argument('--input', type=str, help='Input file')
    parser.add_argument('--sep', type=str, default=',', help='Separator used in the input files')
    parser.add_argument('--expression', type=str, help='Expression according to which the table should be filtered (compatible with pandas query)')
    args = parser.parse_args()

    print(f'File used: {args.input}', file = sys.stderr)

    # Load table
    df = pd.read_csv(args.input, sep = args.sep, index_col=False)
    print(f'Original file: {df.shape[0]} rows, {df.shape[1]} columns.', file = sys.stderr)

    # Filter
    df_filtered = df.query(args.expression)
    print(f'Filtered file: {df_filtered.shape[0]} rows, {df_filtered.shape[1]} columns.', file = sys.stderr)
    # Print filtered output
    df_filtered.to_csv(sys.stdout, index = False, sep = args.sep) 

if __name__ == "__main__":
    main()