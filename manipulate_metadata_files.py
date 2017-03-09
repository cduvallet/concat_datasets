"""
Read in all metadata files in given directory.
Add 'dataset_id' column to all individual files.
Check that there are no duplicate samples across all studies.
    Print error if there is (but don't fix it... TBD how to do that)
Concatenate all metadata files into one large metadata file.
Write that file.
"""

import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('in_dir', help='directory with all metadata files to be read')
parser.add_argument('out_dir', help='directory to write big metadata file to')
parser.add_argument('out_file', help='output file name')
args = parser.parse_args()

files = os.listdir(args.in_dir)

all_metas = []
all_samples = {}
for f in files:
    df = pd.read_csv(os.path.join(args.in_dir, f), sep='\t', index_col=0)
    # metadata files are named dataset_id.metadata.txt
    dataset = f.split('.')[0]
    df['dataset_id'] = dataset
    df.index = [dataset + '--' + str(i) for i in df.index]
    all_samples[f.split('.')[0]] = list(df.index)
    all_metas.append(df)

# Check if any two datasets overlap samples
datasets = all_samples.keys()
for d1 in datasets:
    for d2 in datasets[datasets.index(d1)+1:]:
        common_samples = [i for i in all_samples[d1] if i in all_samples[d2]]
        if len(common_samples) > 0:
            print('WARNING: {} and {} have {} overlapping sample IDs'.format(d1, d2, len(common_samples)))

# Concatenate all metadata and write to file
all_metas_df = pd.concat(all_metas)
all_metas_df.to_csv(os.path.join(args.out_dir, args.out_file), sep='\t')
