"""
Relabel sequences in *.raw_trimmed.fasta to have
datasetID--sampleID_N
"""

import os
import argparse
import util

parser = argparse.ArgumentParser()
parser.add_argument('trimmed_dir', help='directory with *.raw_trimmed.fasta files')
args = parser.parse_args()

files = [os.path.join(args.trimmed_dir, i) for i in os.listdir(args.trimmed_dir) if i.endswith('.raw_trimmed.fasta')]

for fname in files:
    dataset = fname.split('/')[-1].split('.')[0]
    with open(fname + '.relabeled', 'w') as fnew:
        for sid, seq in util.iter_fst(fname):
            sid = '>' + dataset + '--' + sid[1:]
            fnew.write('\n'.join([sid, seq]) + '\n')
