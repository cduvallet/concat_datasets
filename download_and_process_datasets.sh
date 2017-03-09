#!/usr/bin/env bash

## Usage stuff
display_help() {
    echo
    echo "Usage: -f S3_file.txt -d"
    echo
    echo "   S3 file has a list of directories to download from mbit.storage.bucket1"
    echo "   each of these should have a corresponding summary file in almlab.bucket/duvallet/summary_files/ "
    echo "   S3_file.txt should end in '.txt'"
    echo
    echo " -d flag indicates whether to dereplicate each individual dataset. Default is False"
    echo
    exit 1
}

if [[ "$1" == "-h" || "$1" == "--help" ]] ; then 
    display_help
fi 


## Options
download='False'    # Whether to re-download and re-process datasets from S3
derepall='False'    # Whether to de-replicate each individual dataset
#catall='False'      # concatenate all raw_trimmed.fastas and metadata files
#processall='False'  # Process the big ol' concatenated raw_trimmed
#
while getopts 'f:d' flag; do
  case "${flag}" in
    f) download="${OPTARG}" ;;
    d) derepall='True' ;;
#    c) catall='True' ;;
#    p) processall='True' ;;
    *) error "Unexpected option ${flag}" ;;
  esac
done


## For each given S3 bucket, download and process the dataset to raw_trimmed fastas
if [ "$download" != 'False' ]; then
    echo -e "Downloading all datasets in ${download} from S3"
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
    done < $download
fi

### Concatenate all of the raw_trimmed.fasta files
#cat data/raw_trimmed/* > data/for_pipeline/${1%.txt}.all_raw_trimmed.fasta

## Concatenate all the metadata and add dataset_id column
# Also check that there are no duplicate sample IDs.
# If there are duplicate IDs, this code just prints that to screen and keeps going...
# TODO: rename samples in metadata to match what they are in big OTU table (dataset--sampleID)
echo -e "Concatenating metadata files..."
python manipulate_metadata_files.py data/metadata/ data/for_pipeline/ ${1%.txt}.all_metadata.txt

### Dereplicate that huge raw_trimmed.fasta file
### Or maybe just run it through the pipeline??
#python ~/scripts/Master.py -i ~/users/duvallet/pop_study_data/data/for_pipeline/

## Dereplicate each dataset individually
# TODO: add flag in dereplicate...py for min_count (i.e. min reads in one dataset)
trimmed_dir=data/raw_trimmed
derep_dir=data/derep_datasets
if [ "$derepall" != 'False' ]; then
    echo -e "Dereplicating individual datasets..."
    python dereplicate_individual_datasets.py -d -l $trimmed_dir $derep_dir
fi

## Concatenate all the relabeled raw_dereplicated files
concat_dir=data/derep_concat
concat_trimmed=${concat_dir}/dereped_datasets_concated.raw_trimmed.fasta
cat $(find ${derep_dir}/*.raw_dereplicated.fasta) > $concat_trimmed

## Dereplicate the concatenated dereplicated files
# TODO: add flag or variable to change min_count here (i.e. min # of datasets)
echo -e "Dereplicating the concatenated dereplicated files..."
raw_trimmed=${concat_dir}/dereped_datasets_concated.raw_trimmed.fasta
derep_map=${concat_dir}/dereped_datasets_concated.map
derep_fasta=${concat_dir}/dereped_datasets_concated.raw_dereplicated.fasta
proc_summary=${concat_dir}/dereped_datasets_concated.proc_summary.txt
python ~/scripts/3.dereplicate.py -f $raw_trimmed -s _ -o $derep_map -d $derep_fasta -P $proc_summary -M 1

# Update the re-dereplicated concatenated file with total sequence size
# And write in descending size order (bc that's what usearch wants)
map_in=$derep_map
fasta_in=$derep_fasta
fasta_out=${concat_dir}/dereped_datasets_concated.raw_dereplicated.fasta.relabled_and_sorted
python update_concated_derep_fasta.py $map_in $fasta_in $fasta_out

## Cluster the doubly-dereplicated fasta with usearch
echo -e "Clustering the overall dereplicated reads with usearch..."
fasta_in=$fasta_out
otu_seqs_fasta=${concat_dir}/dereped_datasets_concated.otu_seqs.fasta
clustering_results=${concat_dir}/dereped_datasets_concated.clustering_results.tab
usearch8 -cluster_otus $fasta_in -otus $otu_seqs_fasta -otu_radius_pct 0.97 -sizein -uparseout $clustering_results

## Parse the usearch results and re-provenance everything to make massive OTU table
echo -e "Reprovenance everything to make the overall OTU table..."
otu_table=v4_datasets.otu_table.txt
python reprovenance_all_files.py $clustering_results $derep_map $derep_dir $otu_table
