import sys
import argparse
import pandas as pd
import os

# This script combines multiple tables and combines them using the pands.concat() function
# The first column is assumed to be the index and the first row is assumed to be the header

def concatenate_files(files, sep, axis):
  # List to store DataFrames
  df_list = []
  nfiles = 0
  
  print(f'Reading {len(files)} arguments...', file = sys.stderr)
  # Loop through each file path provided as arguments
  for path in files:
      df = pd.read_csv(path, sep = sep, index_col=0)
      nfiles += 1
      print(f'Reading file #{nfiles}, {df.shape[0]} rows, {df.shape[1]} columns.', file = sys.stderr)
      df_list.append(df)
      # If concatenating by column, get the first files index name to use for output (otherwise csv first column doesnt have a header)
      if nfiles==1:
          ind = df.index.name
          print(ind, file = sys.stderr)
  
  # Concatenate all dataframes
  result = pd.concat([df_list[i] for i in range(len(df_list))], axis=axis)
  # get first column name back
  result = result.reset_index().rename(columns={'index': ind})
  print(f'Resulting table contains {result.shape[0]} rows, {result.shape[1]} columns', file = sys.stderr)
  return result

def main():
    parser = argparse.ArgumentParser(description='Concatenate multiple tables, either by column or by row.')
    parser.add_argument('--sep', type=str, default=',', help='Separator used in the input files')
    parser.add_argument('--axis', type=str, default="columns", help='Axis to concatenate along ("rows" or "columns")')
    parser.add_argument('files', metavar='F', type=str, nargs='+', help='Files to concatenate')
    
    args = parser.parse_args()
    
    if args.axis == "rows":
      axis = 0
    elif args.axis == "columns":
      axis = 1
    else:
      raise Exception('--axis can be either "rows" or "columns"')
    
    print(f'Combining {args.axis} from {len(args.files)} files. The separator is: "{args.sep}".', file = sys.stderr)
    
    result = concatenate_files(args.files, args.sep, args.axis)
    result.to_csv(sys.stdout, index = False, sep = args.sep)

if __name__ == "__main__":
    main()
