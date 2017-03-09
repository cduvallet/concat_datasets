"""
This file does a lot of re-provenancing.

First, it reads in the usearch clustering results which
clustered the doubly dereplicated reads (i.e. each individual
dataset's dereplicated reads, concatenated together, and then
re-dereplicated). This results in a {seqID: OTU_ID} map.

It also reads the corresponding dereplication map to identify
which of the original dataset's dereplicated sequences
correspond to these doubly dereplicated seqIDs. In other words,
it returns {seqID: {dataset: origID}}.

Then it combines this into {otuID: {dataset: [list of origIDs]}}

Finally, it parses each individual dataset's original 
dereplication map and uses the above mappings to make a
{OTU_ID: {dataset--sample1: counts, dataset--sample2: counts}}
dictionary for each master OTU.
"""
import os
import argparse
import pandas as pd

def parse_clustering_results(cluster_file):
    """
    Parses a clustering_results output from usearch into a
    dictionary mapping original sequences to OTUs.

    Parameters
    ----------
    cluster_file : str
        file path for clustering results. File should have
        five columns:
            seqID    otu    *    *    *
            seqID    match  99.3 *    match_seqID

    Returns
    -------
    otudict : dict
        {seqID: OTU_ID}
    """
    with open(cluster_file, 'r') as f:
        lines = f.readlines()

    lines = [l.strip().split('\t') for l in lines]
    otudict = {}

    for line in lines:
        seqid = line[0].split(';')[0]
        if line[1] == 'otu':
            otuid = seqid
            otudict[otuid] = [otuid]
        elif line[1] == 'match':
            otuid = line[4]
            otudict[otuid].append(seqid)

    return otudict

def parse_master_derep_map(derep_map):
    """
    Parse the dereplication map which came from dereplicating
    the concatenated dataset-wise dereplicated reads.

    Parameters
    ----------
    derep_map : str
        file with format:
        seqID    origID--datasetID;size=123:1 origID--datasetID;size=142:1

    Returns
    -------
    derep_dict : dict
        dictionary with {seqID: {dataset: origID, dataset2: origID}}
    """

    with open(derep_map, 'r') as f:
        lines = f.read().splitlines()

    derep_dict = {}
    for line in lines:
        line = line.split('\t')
        seqid = line[0]
        datasets = {_split_dataset_derep(s): s.split('--')[0] for s in line[1].split(' ')}
        derep_dict[seqid] = datasets

    return derep_dict

def _split_dataset_derep(s):
    """
    Splits a string in format 'origID--datasetID;size=123:1'
    and returns datasetID, with '-' replaced by underscores
    """
    return '_'.join(s.split('--')[1].split(';')[0].split('-'))

def update_derep_dict(otudict, derepdict):
    """
    Combines otudict and derepdict to return a dictionary
    with {OTU_ID: {dataset: [list_of_orig_IDs]}

    Parameters
    ----------
    otudict : dict
        dictionary with {seqID: OTU_ID}

    derepdict : dict
        dictionary with {seqID: {dataset: orig_ID, dataset: orig_ID}}
    
    Returns
    -------
    newdict : dict
        {OTU_ID: {dataset: [list_of_orig_IDs]}
    """
    newdict = {}

    for otu in otudict:
        if otu not in newdict:
            newdict[otu] = {}
        
        for seq in otudict[otu]:
            # derepdict[seq][dataset] is a string, orig_ID
            for dataset in derepdict[seq]:
                if dataset in newdict[otu]:
                    # derepdict[seq] is a dict where dataset is a key
                    newdict[otu][dataset].append(derepdict[seq][dataset])
                else:
                    newdict[otu][dataset] = [derepdict[seq][dataset]]

    return newdict

def parse_one_map_file(map_file):
    """
    Parse one normal dereplication map

    Parameters
    ----------
    map_file : str
        file with format:
        seqID    s1:counts s2:counts s4:counts

    Returns
    -------
    derep_dict : dict
        {seqID: {s1: counts, s2:counts, s4: counts}}
    """

    with open(map_file, 'r') as f:
        derepmap = f.readlines()
    derepmap = [l.strip().split('\t') for l in derepmap]

    return {l[0]:  {i.split(':')[0]: i.split(':')[1] for i in l[1].split(' ')} for l in derepmap}


def read_dataset_derep_maps(derep_dir):
    """
    Reads all datasetID.map files in derep_dir and parses into
    a master dereplication dictionary.

    Parameters
    ----------
    derep_dir : str
        path to directory containing all dereplication maps,
        which should be labeled dataset.map

    Returns
    -------
    dataset_maps : dict
        {dataset:  {origID: {s1: count, s2: count},
                    origID2: {s1: count, s2: count}
                   },
         dataset2: {origID: {s1: count, s2: count},                                                                                                                                                                                             origID2: {s1: count, s2: count}
                   }
        }
    """
    map_files = [i for i in os.listdir(derep_dir) if i.endswith('.map')]
    # make sure the dereped_datasets_concated.map file doesn't get read in
    map_files = [i for i in map_files if not i.startswith('dereped')]

    dataset_maps = {}
    for map_file in map_files:
        dataset = map_file.split('.')[0]
        dataset_maps[dataset] = parse_one_map_file(os.path.join(derep_dir, map_file))

    return dataset_maps

def collapse_derep_map(dataset, list_of_orig_IDs, dataset_derep_map):
    """
    Given a dataset, a list of original sequence IDs in that dataset,
    and a dereplication dict for that dataset, return the total counts
    in each sample (i.e. sum of counts in all seqs in list_of_orig_IDs).
    Also rename the samples to dataset--sample_name

    Parameters
    ----------
    dataset : str
        name of dataset, will be prepended to all sample IDs
    list_of_orig_IDs : list
        list of sequence IDs (should be keys in dataset_derep_map)
        to sum counts for across all samples
    dataset_derep_map : dict
        {seq1: {s1: counts, s2: counts, ...},
         seq2: {s2: counts, s10: counts, ...}, ...}

    Returns
    -------
    collapsed : dict
        Single-level dictionary with relabeled sample IDs (dataset + '--' + original
        sample ID) and total counts across all seqs in list_of_orig_IDs. e.g.:
            {dataset--s1: counts, dataset--s2: counts, dataset--s3: counts, ...}
    """
    
    collapsed = {}
    for seq in list_of_orig_IDs:
        for sample in dataset_derep_map[seq]:
            try:
                collapsed[sample] += float(dataset_derep_map[seq][sample])
            except:
                collapsed[sample] = float(dataset_derep_map[seq][sample])
    
    # Relabel samples (i.e. keys) in collapsed
    collapsed = {dataset + '--' + k: collapsed[k] for k in collapsed}
    return collapsed


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('cluster_file', help='path to cluster_results.tab file')
    parser.add_argument('derep_map', help='dereplication map, indicating which datasets each of the sequences in cluster_file were found in')
    parser.add_argument('derep_dir', help='directory with dataset-wise dereplication maps, labeled datasetID.map')
    parser.add_argument('table_out', help='file name for output OTU table')
    args = parser.parse_args()
        
    ## Parse clustering results to {seqID: OTU_ID}
    print("Parsing clustering results...")
    otudict = parse_clustering_results(args.cluster_file)

    ## Get the OTU_ID --> dataset --> original_seq_IDs map
    # Parse dereplication map to {seqID: {dataset: orig_ID, dataset: orig_ID}}
    print("Parsing master dereplication map..")
    derepdict = parse_master_derep_map(args.derep_map)

    # Combine these two dicts to get {OTU_ID: {dataset: [list_of_orig_IDS], dataset: [list_orig_IDs]}}
    derepOTUdict = update_derep_dict(otudict, derepdict)

    ## Convert that dict into one with {OTU_ID: {dataset: {s1: total_otu_counts, s2: total_otu_counts...}}}
    
    ## Read in all of the dataset dereplication maps into a dictionary
    # dataset_maps is {dataset: {seq: {s1: counts, s2: counts}}}
    print("Reading all dataset dereplication maps...")
    dataset_maps = read_dataset_derep_maps(args.derep_dir)

    # Use that dict to map each OTU_ID to sample abundances in each dataset
    # i.e. return {OTU_ID: {dataset--s1: total_counts, dataset--s2: total_counts}}
    print("Collapsing..."),
    finaldict = {}
    for otu in derepOTUdict:
        for dataset in derepOTUdict[otu]:
            try:
                finaldict[otu].update(collapse_derep_map(dataset, 
                                                         derepOTUdict[otu][dataset], # [list_of_orig_IDS]
                                                         dataset_maps[dataset])) # {seq: {s1: counts, s2: counts}}
            except:
                finaldict[otu] = (collapse_derep_map(dataset, 
                                                     derepOTUdict[otu][dataset], # [list_of_orig_IDS]
                                                     dataset_maps[dataset])) # {seq: {s1: counts, s2: counts}}

    # This takes a long time. Is there a better way to save this final dict??
    df = pd.DataFrame.from_dict(finaldict)
    df.to_csv(args.table_out, sep='\t')

