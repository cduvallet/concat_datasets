import argparse
from SummaryParser import *

parser = argparse.ArgumentParser()
parser.add_argument('summary_file', help='summary_file.txt', default='summary_file.txt')
args = parser.parse_args()

## Read summmary file
sumobj = SummaryParser(args.summary_file)
sumobj.ReadSummaryFile()

## Set 16S attributes
# atts is a dict
atts = sumobj.attribute_value_16S
atts['OTU_SIMILARITY'] = "97"
atts['TRIM_LENGTH'] = "150"
atts['MIN_COUNT'] = "10"

# remove quality trim bc otherwise it's the default
if 'QUALITY_TRIM' in atts:
    del atts['QUALITY_TRIM']
atts['MAX_ERRORS'] = "2"

atts['PROCESSED'] = "False"

## Update summary file and write
sumobj.attribute_value_16S = atts
sumobj.WriteSummaryFile()
