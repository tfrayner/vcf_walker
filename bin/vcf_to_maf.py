#!/usr/bin/env python
#
# $Id: vcf_to_maf.py 4006 2017-03-31 09:47:13Z tfrayner $

'''
Script converts a VCF to MAF, using a pickled VcfGeneModel object
to infer variant effects and affected gene IDs.
'''

import sys
import logging

from vcf_walker.vcf2maf import Vcf2Maf, LOGGER

################################################################################

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(\
    description="Convert a VCF to MAF as best we can.")

  PARSER.add_argument('-i', '--input-file', dest='infile', type=str,
                      required=True, help="The name of the input VCF file.")

  PARSER.add_argument('-o', '--output-file', dest='outfile', type=str,
                      required=True, help="The name of the output VCF file.")

  GROUP = PARSER.add_mutually_exclusive_group(required=True)

  GROUP.add_argument('-f', '--fasta', dest='fasta', type=str,
                     help='The FASTA file containing the whole genome'
                     + ' sequence used to build the Vcf2Maf object')

  GROUP.add_argument('-p', '--pickle', dest='pickle', type=str,
                     help='The name of a pickle file containing the saved'
                     + ' Vcf2Maf object to reload.')

  PARSER.add_argument('-g', '--gtf', dest='gtf', type=str, required=False,
                      help='The name of the GTF file containing gene model'
                      + ' (CDS, exon) info. Only required if rebuilding a'
                      + ' Vcf2Maf object from scratch.')

  PARSER.add_argument('-s', '--savefile', dest='savefile', type=str,
                      required=False,
                      help='The name of the file to save a pickled Vcf2Maf'
                      + ' object to (for later reloading; this is much faster'
                      + ' than rebuilding from scratch every time).')

  PARSER.add_argument('-r', '--random', dest='random', action='store_true',
                      help='Indicates that the script should generate'
                      + ' randomised SNVs in the output. This may be used'
                      + ' for generating null distributions.')

  PARSER.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                      help='Increase verbosity of the output (substantially).'
                      + ' This option will include the calculated peptide'
                      + ' sequences in the logging output.')

  ARGS = PARSER.parse_args()

  if ARGS.verbose:
    LOGGER.setLevel(logging.INFO)

  VCFCONV = Vcf2Maf(picklefile = ARGS.pickle,
                    fasta      = ARGS.fasta,
                    gtf        = ARGS.gtf,
                    savefile   = ARGS.savefile)

  VCFCONV.convert(infile  = ARGS.infile,
                  outfile = ARGS.outfile)

  sys.stderr.write("Vcf conversion complete.\n")
