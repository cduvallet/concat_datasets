#!/usr/bin/env bash

## Usage stuff
display_help() {
    echo
    echo "Usage: $0 S3_file.txt"
    echo
    echo "   S3 file has a list of directories to download from mbit.storage.bucket1"
    echo "   each of these should have a corresponding summary file in almlab.bucket/duvallet/summary_files/ "
    echo "   S3_file.txt should end in '.txt'"
    echo
    exit 1
}

if [[ "$1" == "-h" || "$1" == "--help" ]] ; then 
    display_help
fi 


### Options
#download='False'    # Whether to re-download and re-process datasets from S3
#catall='False'      # concatenate all raw_trimmed.fastas and metadata files
#processall='False'  # Process the big ol' concatenated raw_trimmed
#
#while getopts 'f:cp' flag; do
#  case "${flag}" in
#    f) download="${OPTARG}" ;;
#    c) catall='True' ;;
#    p) processall='True' ;;
#    *) error "Unexpected option ${flag}" ;;
#  esac
#done


## For each given S3 bucket, download and process the dataset to raw_trimmed fastas
while read bucket;
do
    # Download data from S3
    aws s3 cp s3://mbit.storage.bucket1/${bucket} tmp_data/ --recursive

    # Download summary file from almlab S3
    aws s3 cp s3://almlab.bucket/duvallet/summary_files/${bucket}.summary_file.txt tmp_data/summary_file.txt

    # Call python script to edit summary file
    python update_summary_file.py tmp_data/summary_file.txt
    dataset=$(cut -f 2 tmp_data/summary_file.txt | head -n 1)
    
    # Run Master.py
    python ~/scripts/Master.py -i ~/users/duvallet/pop_study_data/tmp_data/

    # Move raw_trimmed.fasta and summary file somewhere safe
    mv ~/proc/${dataset}_proc_16S/${dataset}.raw_trimmed.fasta data/raw_trimmed/

    # Also grab the metadata file, put it here too
    aws s3 cp s3://almlab.bucket/duvallet/metadata_files/${dataset}.metadata.txt data/metadata/

    # Delete the data, proc, and processing_results folders
    rm -r ~/proc/${dataset}_proc_16S
    rm -r ~/processing_results/${dataset}_results
    rm tmp_data/*
done < $1

## Concatenate all of the raw_trimmed.fasta files
cat data/raw_trimmed/* > data/for_pipeline/${1%.txt}.all_raw_trimmed.fasta

## Concatenate all the metadata and add dataset_id column
# Also check that there are no duplicate sample IDs.
# If there are duplicate IDs, this code just prints that to screen and keeps going...
python manipulate_metadata_files.py data/metadata/ data/for_pipeline/ ${1%.txt}.all_metadata.txt

## Dereplicate that huge raw_trimmed.fasta file
## Or maybe just run it through the pipeline??
python ~/scripts/Master.py -i ~/users/duvallet/pop_study_data/data/for_pipeline/


