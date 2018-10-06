#!/usr/bin/env python

import argparse
## Parse command line arguments
parser = argparse.ArgumentParser(description='Count allele frequencies for single cells.')
parser.add_argument("-v", "--verbose", help="set verbose mode", action="store_true")
parser.add_argument("--snps", help="SNPs file, tab separated file with no header and 4 columns: chr, pos (0-indexed), reference base, alternate base",required=True);
parser.add_argument("--barcodes", help="barcodes file containing only barcodes of interest from the CB tag in the bam file; one barcode per line",required=True);
parser.add_argument("--output-format", help="output format for resulting tables; mm for Matrix Market, hdf for HDF", default="mm", choices=['mm'], dest="output_format"); # TODO: add hdf
parser.add_argument("bamfile", help="sorted and indexed bamfile from single cell experiment")
parser.add_argument("--max-depth", help="max_depth argument for fileup", type=int, default=99999999, dest="maxdepth")
parser.add_argument("--output-prefix", help="prefix of output files", default="", dest="prefix")
args = parser.parse_args()

## This takes a while, delays help if at the top
import pandas as pd
import pysam
import sys
import scipy
import scipy.io
import scipy.sparse
import numpy
import os.path

## Read positions
if args.verbose:
    sys.stdout.write("reading snps file... ")
    sys.stdout.flush();
if not os.path.exists(args.snps):
    print("Error: snps file does not exist");
    sys.exit(1);
with open(args.snps, 'r') as f:
    postable = pd.read_table(f, header=None, sep='\t', names=['chr','pos','ref','alt'])
    postable['pos'] = postable['pos'].astype(int)
if args.verbose:
    sys.stdout.write("done\n")
    sys.stdout.flush();
postable["posindex"] = postable["chr"].map(str) + ':' +  postable["pos"].map(str)
    
if args.verbose:
    sys.stdout.write("building positions index... ")
    sys.stdout.flush();
pos_index = {}
for index, row in postable.iterrows():
    pos_index[row['posindex']] = index;
if args.verbose:
    sys.stdout.write("done\n")
    sys.stdout.flush();

## Read barcodes
if args.verbose:
    sys.stdout.write("reading barcodes file... ");
    sys.stdout.flush();
if not os.path.exists(args.barcodes):
    print("Error: barcodes file does not exist");
    sys.exit(1);
with open(args.barcodes,'r') as f:
    barcodetable = pd.read_table(f, header=None, sep=' ', names=['barcode'])
if args.verbose:
    sys.stdout.write("reading barcodes file... ");
    sys.stdout.write("done\n")

if args.verbose:
    sys.stdout.write("building barcodes index... ");
    sys.stdout.flush()
barcode_index = {}
for index, row in barcodetable.iterrows():
    barcode_index[row['barcode']] = index;
if args.verbose:
    sys.stdout.write("done\n");
    
pd.DataFrame({'index': list(barcode_index.values()),
              'barcode': list(barcode_index.keys())}).to_csv(args.prefix + '-barcode_index.csv')

pd.DataFrame({'index': list(pos_index.values()),
              'region': list(pos_index.keys())}).to_csv(args.prefix + '-region_index.csv')



## Generate dok arrays to store the data
# cells x positions
ncells = barcodetable.shape[0];
npositions = postable.shape[0];
# dok arrays allow fast access
refarray = scipy.sparse.dok_matrix((ncells,npositions),'uint32')
altarray = scipy.sparse.dok_matrix((ncells,npositions),'uint32')
covarray = scipy.sparse.dok_matrix((ncells,npositions),'uint32')

#####################3
    # cells x positions
    # ncells = barcodetable.shape[0];
    # npositions = postable.shape[0];
    #postable["posindex"]
    #barcodetable["barcode"]
# for debug


# import pdb; pdb.set_trace()


# import csv

# # write it
# with open('rownames.csv', 'w') as csvfile:
#     writer = csv.writer(csvfile)
#     [writer.writerow(r) for r in postable]

# import pdb; pdb.set_trace()

####################


# start pileup
if args.verbose:
    sys.stdout.write("perfoming pileup...");
    sys.stdout.flush();

# open the bam file (must be sorted and indexed)
print("Creating alignment file")
samfile = pysam.AlignmentFile(args.bamfile,"rb")
print("Done")

print(pos_index.items()[0:10])


# do pileups
i=0;
for pileupcolumn in samfile.pileup(max_depth=args.maxdepth):
    i += 1;
    
    # Everything is 0 based
    posindex = str(pileupcolumn.reference_name) +':' + str(pileupcolumn.reference_pos); 
    
    if ( args.verbose and i % 10000000 == 0):
        sys.stdout.write('At region ' + posindex + "\n");
        sys.stdout.flush();
    
    if posindex in pos_index:

        posi = pos_index[posindex];
        refallele = postable.iloc[posi,]['ref']
        altallele = postable.iloc[posi,]['alt']
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                try:
                    cbtagvalue  = pileupread.alignment.get_tag('CB')
                    if cbtagvalue in barcode_index:
                        barcodei = barcode_index[cbtagvalue];
                        readbase = pileupread.alignment.query_sequence[pileupread.query_position];
                        covarray[barcodei,posi] += 1;
                        if (refallele == readbase):
                            refarray[barcodei,posi] += 1;
                        elif (altallele == readbase):
                            altarray[barcodei,posi] += 1;
                except:
                    pass

if args.verbose:
    sys.stdout.write("done\n");

## close the sam file
samfile.close()

# Save coverage as triplets
if args.verbose:
    sys.stdout.write('saving matrices...');
    sys.stdout.flush()

if args.output_format == 'mm':
    # cells x positions
    # ncells = barcodetable.shape[0];
    # npositions = postable.shape[0];
    #postable["posindex"]
    #barcodetable["barcode"]
    
    # todo save the names
    
    scipy.io.mmwrite(target=args.prefix + "covmat.mtx", a=covarray, field='integer');
    scipy.io.mmwrite(target=args.prefix + "refmat.mtx", a=refarray, field='integer');
    scipy.io.mmwrite(target=args.prefix + "altmat.mtx", a=altarray, field='integer');
else:
    sys.stdout.write('unsupported output format');
    sys.stdout.flush()
    
if args.verbose:
    sys.stdout.write('done\n');
    sys.stdout.flush()

# terminate
sys.exit(0)
