"""
This script reads in all the OTU tables and metadata files from a folder
with the clean tables, converts to relative abundance, collapses to genus
level, relabels sample IDs, and concatenates into bigmeta and bigdf files.
"""
import feather
import pandas as pd
import os
import argparse

import copy

def collapse_taxonomic_contents_df(OTU_table, taxonomic_level):
    """
    Collapses OTU table to given taxonomic level by string-matching.

    Adapated from Thomas Gurry's code at
    https://github.com/thomasgurry/amplicon_sequencing_pipeline

    Parameters
    ----------
    OTU_table : pandas dataframe
        OTUs in columns, samples in rows.
        Taxonomic levels in OTU strings should be semicolon-delimited,
        starting with kingdom level.
        Unannotated taxonomic levels should end with '__' (e.g. ''...;g__Roseburia;s__')
    taxonomic_level : str
        kingdom, phylum, class, order, family, genus, or species

    Returns
    -------
    newdf : pandas dataframe
        OTUs in columns, samples in rows.
        OTUs are collapsed to the given taxonomic level.
        Matching values (for annotated taxa) are summed for each sample.
        Values corresponding to unannotated taxa are discarded.
    """

    OTU_IDs = list(OTU_table.columns)

    # Collapse to the right level
    if(taxonomic_level == "kingdom"):
        OTU_taxa = [OTU_ID.split(';')[0] for OTU_ID in OTU_IDs]
    if(taxonomic_level == "phylum"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:2]) for OTU_ID in OTU_IDs]
    if(taxonomic_level == "class"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:3]) for OTU_ID in OTU_IDs]
    if(taxonomic_level == "order"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:4]) for OTU_ID in OTU_IDs]
    if(taxonomic_level == "family"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:5]) for OTU_ID in OTU_IDs]
    if(taxonomic_level == "genus"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:6]) for OTU_ID in OTU_IDs]
    if(taxonomic_level == "species"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:7]) for OTU_ID in OTU_IDs]

    # Get indices of each unique taxon
    taxa_indices = {}
    for i in range(len(OTU_taxa)):
        if (OTU_taxa[i] not in taxa_indices) and (OTU_taxa[i][(len(OTU_taxa[i])-2):] != "__"):
            taxa_indices[OTU_taxa[i]] = []
        if (OTU_taxa[i][(len(OTU_taxa[i])-2):] != "__"):
            taxa_indices[OTU_taxa[i]].append(i)
    # Make new empty df with the same samples as original df and taxa in taxa_indices.keys()
    newdf = pd.DataFrame(index=OTU_table.index, columns=taxa_indices.keys(), data=0)

    # Get sample contents for each taxa of the chosen level and put into newdf
    for key in taxa_indices:
        indices = taxa_indices[key]
        newcol = OTU_table.iloc[:, indices].sum(axis=1)
        newdf[key] = copy.copy(newcol)

    return newdf

p = argparse.ArgumentParser()
p.add_argument('clean_dir', help='path to directory with clean '
    + 'OTU and metadata tables in feather format. Files should be labeled '
    + ' dataset.otu_table.clean.feather or dataset.metadata.clean.feather.')
p.add_argument('otu_out', help='file name of output OTU table')
p.add_argument('meta_out', help='file name of output metadata file')
args = p.parse_args()

files = os.listdir(args.clean_dir)
# Files must have feather suffix
files = [i for i in files if i.endswith('.feather')]
datasets = list(set([i.split('.')[0] for i in files]))

alldfs = []
allmetas = []

for d in datasets:
    print(d)
    fnotu = d + '.otu_table.clean.feather'
    fnmeta = d + '.metadata.clean.feather'

    df = feather.read_dataframe(os.path.join(args.clean_dir, fnotu))
    # Feather format does not support index names, first column has index
    df.index = df.iloc[:,0]
    df = df.iloc[:, 1:]

    meta = feather.read_dataframe(os.path.join(args.clean_dir, fnmeta))
    meta.index = meta.iloc[:, 0]
    meta = meta.iloc[:, 1:]

    # Convert to relative abundance
    df = df.divide(df.sum(axis=1), axis=0)
    df = collapse_taxonomic_contents_df(df, 'genus')

    # Relabel sample IDs
    df.index = [str(i) + '/' + d for i in df.index]
    meta.index = [str(i) + '/' + d for i in meta.index]

    alldfs.append(df)
    allmetas.append(meta)

bigdf = pd.concat(alldfs, axis=0)
bigmeta = pd.concat(allmetas, axis=0)

bigdf = bigdf.fillna(0.0)

bigdf.to_csv(args.otu_out, sep='\t')
bigmeta.to_csv(args.meta_out, sep='\t')
