#!/usr/bin/env python

'''
Quick script to pull out the GC content and overall CDS length for
each feature in a VcfAnnotator object.
'''

import sys
import os
from vcf_walker.annotator import VcfAnnotator
from vcf_walker.utils import flexi_open
from Bio.SeqUtils import GC

def extract_metadata(vcfobj, outfh):

  for (chrom, chrdat) in vcfobj.chromosomes.iteritems():
    for feature in chrdat.features:

      seq     = feature.extract(chrdat.seq.tostring())
      seqlen  = len(seq)
      gcfract = GC(seq)/100

      if 'transcript_id' in feature.qualifiers:
        trnsid = feature.qualifiers['transcript_id'][0]
      elif hasattr(feature, 'id'):
        trnsid = getattr(feature, 'id')
      else:
        trnsid = 'NA'

      outfh.write(",".join([str(x) for x
                            in (trnsid, seqlen, gcfract)]) + "\n")

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(\
    description="Retrieve CDS metadata from a pickled VcfAnnotator object.")
  
  PARSER.add_argument('-o', '--output', dest='outfile', type=str, required=True,
                      help='The name of the output file in which statistics'
                      + ' are to be dumped (CSV format).')

  PARSER.add_argument('-p', '--pickle', dest='pickle', type=str, required=True,
                      help='The name of a pickle file containing the saved'
                      + ' VcfAnnotator object to reload.')

  ARGS = PARSER.parse_args()

  if os.path.exists(ARGS.outfile):
    raise StandardError("Will not overwrite pre-existing output file %s"
                        % ARGS.outfile)

  VCFANN = VcfAnnotator(picklefile = ARGS.pickle)

  with open(ARGS.outfile, 'w') as outfh:
    extract_metadata(VCFANN, outfh)
