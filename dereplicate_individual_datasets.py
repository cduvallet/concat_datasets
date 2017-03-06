"""
This script
1. dereplicates each dataset's raw_trimmed.fasta 
1. relabels the dereplicated seqIDs to track the original dataset/seqID combination
...

"""

import os
import subprocess
import argparse

### Arguments
parser = argparse.ArgumentParser()
parser.add_argument('in_dir', help='directory with raw_trimmed fasta files')
parser.add_argument('out_dir', help='out directory (where to write resulting files)')
parser.add_argument('-d', help='dereplicate reads', action='store_true', default=False)
parser.add_argument('-l', help='relabel map files and raw_dereplicated fastas', action='store_true', default=False)

args = parser.parse_args()

### Get all the raw_trimmed fasta files
files = os.listdir(args.in_dir)
trimmed_fastas = [i for i in files if i.endswith('.raw_trimmed.fasta')]

### Dereplicate them
if args.d:
    print('Dereplicating reads')
    # Prepare commands
    fnamefun = lambda f: (os.path.join(args.in_dir, f),
                          os.path.join(args.out_dir, f.replace('.raw_trimmed.fasta', '') + '.map'),
                          os.path.join(args.out_dir, f.replace('.raw_trimmed.fasta', '') + '.raw_dereplicated.fasta'),
                          os.path.join(args.out_dir, f.replace('.raw_trimmed.fasta', '') + '.proc_summary.txt'))
    cmdfun = lambda f: "python ~/scripts/3.dereplicate.py -f {} -s _ -o {} -d {} -P {}".format(*fnamefun(f))
    commands = [cmdfun(f) for f in trimmed_fastas]

    # Run in parallel
    processes = [subprocess.Popen(cmd, shell=True) for cmd in commands]
    # Make sure they've all finished
    for p in processes:
        p.wait()

if args.l:
    print('Relabeling output files from dereplication')
    ### Relabel sequences with dataset ID
    # Note: dataset ID shouldn't have underscores in them!
    ## This needs to happen in the first column of the provenance map (*.map)
    # Current sequences are labeled seqID
    # Turn them into seqID-dataset;size=whatever
    files = os.listdir(args.out_dir)
    mapfiles = [os.path.join(args.out_dir, i) for i in files if i.endswith('.map')]
    derepfiles = [os.path.join(args.out_dir, i) for i in files if i.endswith('.raw_dereplicated.fasta')]
    
    for fname in mapfiles:
        print(fname)
        dataset = '-'.join(fname.split('/')[-1].split('.')[0].split('_'))
        with open(fname, 'r') as f:
            lines = f.readlines()
            lines = [l.strip().split('\t') for l in lines]
            lines = [l[0] + '--' + dataset + '\t' +  l[1] for l in lines]
            with open(fname + '.relabled', 'w') as fnew:
                fnew.write('\n'.join(lines))
                fnew.write('\n')

    ## And in the sequence headers of the raw_dereplicated fasta
    # Current sequences are labeled >seqID;size
    # Turn them into >seqID-dataset_1
    # You need to have (1) no underscores in the dataset ID and
    # (2) an underscore after seqID-dataset, so that downstream
    # dereplication sees the sequence handle as the "sample ID"
    for fname in derepfiles:
        print(fname)
        dataset = '-'.join(fname.split('/')[-1].split('.')[0].split('_'))
        lines = []
        with open(fname, 'r') as f:
            for line in f.readlines():
                if line.startswith('>'):
                    # Change >seqID;size=132 to >seqID--dataset-id;size=132_1
                    line = line.split(';')[0] + '--' + dataset + ';' + line.split(';')[1].strip() + '_1'
                lines.append(line.strip())
        with open(fname + '.relabeled', 'w') as f:
            f.write('\n'.join(lines))
            f.write('\n')
