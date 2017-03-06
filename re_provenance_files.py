"""
This file does a lot of provenancing business.

"""

import argparse
import util

parser = argparse.ArgumentParser()
parser.add_argument('map_file', help='dereplication map (seqID    sample:count sample2:count)')
parser.add_argument('fasta_in', help='raw_dereplicated fasta corresponding to map')
parser.add_argument('fasta_out', help='relabled fasta file')
args = parser.parse_args()

## Parse the dereplication map and get total size per sequence
print('Parsing dereplication map: {}'.format(args.map_file))
seq_sizes = {}

with open(args.map_file, 'r') as f:
    lines = f.readlines()

for line in lines:
    line = line.strip().split('\t')
    seqID = line[0]
    total_size = sum([int(i.split('size=')[1].split(':1')[0]) for i in line[1].split(' ')])
    seq_sizes[seqID] = total_size

## Read in the entire fasta file into a dict {seqID: sequence}
print('Reading fasta file: {}'.format(args.fasta_in))
fasta = {}
for sid, seq in util.iter_fst(args.fasta_in):
    # sid in the fasta is something like >444;size=8
    # sid.split(';')[0][1:] returns 444, which is a key in seq_sizes
    sid = sid.split(';')[0][1:]
    newsid = '>' + sid + ';size=' + str(seq_sizes[sid])
    fasta[sid] = {}
    fasta[sid]['new_sid'] = newsid
    fasta[sid]['seq'] = seq

## Get list of sequence IDs in descending size (i.e. largest first)
ordered_seqs = sorted(seq_sizes, key=lambda k: seq_sizes[k], reverse=True)

## Write new fasta file in descending size
print('Writing sorted and relabled fasta: {}'.format(args.fasta_out))
with open(args.fasta_out, 'w') as f:
    f.write('\n'.join([fasta[s]['new_sid'] + '\n' + fasta[s]['seq'] for s in ordered_seqs]))
