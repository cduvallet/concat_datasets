# pop_study_data

This directory contains work related to downloading and processing
all V4 datasets together.

Target applications are Fangqiong's population estimation study and
Eric's batch effects business.

# What's the plan?

1. Write down the S3 datasets that correspond to all the V4 regions (manually, using
`datasets_info.txt` from my meta-analysis work).
2. Download the dataset from S3 mbit bucket
3. Download the summary file from my almlab bucket
   - or modify it so they're all the same?
   - I think modify it, yeah!
4. Process the dataset
5. Grab relevant files...
   - raw_dereplicated.fasta
   - dereplication map
6. Concatenate the raw_dereplicated fastas
7. See how big that file is. Hopefully less that 4 G. If not... dereplicate that again!?!
8. ???
9. Fix batch effects. Profit.

# Processing parameters

## Goal
We'll be making 97% OTUs and assigning taxonomy with RDP.

I think the goal for processing is to find the lowest common denominator
that can work for all datasets. This means that our parameters will likely
be super conservative.

## Datasets
Okay. There are 11 V4 datasets.

2 of these started out with fasta's.

1 used max expecte errors filtering (maxee = 2) (t1d_alkanani)

The minimum trim length was 150.

1 was split by barcodes and had primers removed (t1d_mejialeon)

## Plan
So, I think we should:
use max expected errors = 2
trim length = 150

And while we're at it let's make min count = 10 again

# Processing workflow

Just run `./download_and_process_datasets.sh S3_V4_datasets.txt` and make sure you've copied
`master_summary_file.txt` to `data/for_pipeline`.

The `download_and_process_datasets.sh` script:
1. downloads each dataset from S3
1. updates its summary file with the right parameters
1. processes it through the pipeline
1. copies the resulting `*.raw_trimmed.fasta` to `data/raw_trimmed/`
1. also downloads the corresponding metadata from the `almlab` bucket, puts it in `data/metadata/`
1. concatenates all of the metadata files, into `data/for_pipeline/`
1. dereplicates each raw_trimmed file into a raw_dereplicated file, in `data/derep_data/`
   - this uses the default `min_count = 10`, i.e. throws out sequences which had fewer than 10 reads in each dataset
1. concatenates all of the `raw_trimmed` fastas into `data/derep_data/dereped_datasets_concated.raw_trimmed.fasta`
1. re-dereplicates this file into `data/derep_data/dereped_datasets_concated.raw_dereplicated.fasta`
   - this dereplication step uses `min_count = 1`, i.e. it was in at least one dataset
1. reorders and relabels the sequence IDs according to total size across all studies in `data/derep_data/dereped_datasets_concated.raw_dereplicated.fasta.relabeled_and_sorted`
1. clusters this super de-replicated fasta with usearch

The `download_and_process_datasets.sh` script calls the other two Python scripts:
* `update_summary_file.py` updates a dataset's summary file with the processing parameters (as above)
* `manipulate_metadata_files.py` reads through all the metadata files, checks if there are duplicate 
sample IDs, concatenates all the metadata files, and writes that metadata to `data/for_pipeline`.

 