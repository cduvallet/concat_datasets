#!/usr/bin/env bash

## Usage stuff
display_help() {
    echo
    echo "Usage: $0 -p -c -f S3_file.txt -d"
    echo
    echo " -p          instructs script to concatenate each dataset's *.raw_trimmed.fasta"
    echo "             and process this file as a RAW_FASTA_FILE input to the pipeline."
    echo "             If this option is chosen, make sure there's a correct summary_file.txt"
    echo "             in ./data/for_pipeline/. Right now, this script calls the RAW_FASTA_FILE"
    echo "             'all_raw_trimmed.fasta'"
    echo
    echo " -c          instructs script to dereplicate each dataset individually, concatenate"
    echo "             these and re-dereplicate, cluster the concatenated re-dereplicated sequences"
    echo "             and reprovenance original dataset reads to these OTUs."      
    echo
    echo " -f          Indicates file to use to download from S3 and reprocess each individual dataset."
    echo "             S3 file has a list of directories to download from mbit.storage.bucket1."
    echo "             Each of these should have a corresponding summary file in almlab.bucket/duvallet/summary_files/. "
    echo "             If -f is not specified, script assumes that data is already downloaded and processed,"
    echo "             and that corresponding raw_trimmed.fasta files are in ./data/raw_trimmed"
    echo
    echo " -d          indicates whether to dereplicate each individual dataset. Default is False"
    echo "             If this flag is not used, individual dataset reads should be dereplicated in"
    echo "             ./data/derep_datasets directory"
    echo
    exit 1
}

if [[ "$1" == "-h" || "$1" == "--help" ]] ; then 
    display_help
fi 


## Options
download='False'         # Whether to re-download and re-process datasets from S3
derepall='False'         # Whether to de-replicate each individual dataset
pipeline_proc='False'    # Whether to do processing via cat *.raw_trimmed > RAW_FASTA_FILE
derep_proc='False'       # Whether to do processing via derep -> cat derep -> re-derep -> cluster

while getopts 'pcf:d' flag; do
  case "${flag}" in
    f) download="${OPTARG}" ;;
    d) derepall='True' ;;
    p) pipeline_proc='True' ;;
    c) derep_proc='True' ;;
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
    	rm -r tmp_data/*
    done < $download
fi


## Concatenate all the metadata and add dataset_id column
# Also check that there are no duplicate sample IDs.
# If there are duplicate IDs, this code just prints that to screen and keeps going...
# TODO: rename samples in metadata to match what they are in big OTU table (dataset--sampleID)
echo -e "Concatenating metadata files..."
python manipulate_metadata_files.py data/metadata/ data/for_pipeline/ all_metadata.txt

## Dereplicate each dataset individually
# TODO: add flag in dereplicate...py for min_count (i.e. min reads in one dataset)
trimmed_dir=data/raw_trimmed
derep_dir=data/derep_datasets
if [ "$derepall" != 'False' ]; then
    echo -e "Dereplicating individual datasets..."
    python dereplicate_individual_datasets.py -d -l $trimmed_dir $derep_dir
fi

if [ "$derep_proc" == 'True']; then
    ## Concatenate all the relabeled raw_dereplicated files
    echo -e "Concatenating individual raw_dereplicate fasta files..."
    concat_dir=data/derep_concat
    concat_trimmed=${concat_dir}/dereped_datasets_concated.raw_trimmed.fasta
    cat $(find ${derep_dir}/*.raw_dereplicated.fasta.relabeled) > $concat_trimmed

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
fi

if [ "$pipeline_proc" == 'True']; then
    ## Concatenate all of the raw_trimmed.fasta files
    echo -e "Concatenating all of the *raw_trimmed.fasta* files"
    cat data/raw_trimmed/* > data/for_pipeline/all_raw_trimmed.fasta

    ## Dereplicate that huge raw_trimmed.fasta file
    echo -e "Running the concatenated raw_trimmed file through pipeline"
    python ~/scripts/Master.py -i ~/users/duvallet/pop_study_data/data/for_pipeline/
fi
