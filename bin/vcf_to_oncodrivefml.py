#!/usr/bin/env python
#
# $Id: vcf_to_maf.py 4006 2017-03-31 09:47:13Z tfrayner $

'''
Script converts a VCF the tab-delimited format expected by
OncodriveFML. Currently requires genome annotation, although this
requirement may be dropped in subsequent versions.
'''

import sys
import logging

from vcf_walker.vcf2maf import Vcf2OncodriveFML, LOGGER

################################################################################

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(\
    description="Convert a VCF to the tab-delimited input format for OncodriveFML.")

  PARSER.add_argument('-i', '--input-file', dest='infile', type=str,
                      required=True, help="The name of the input VCF file.")

  PARSER.add_argument('-o', '--output-file', dest='outfile', type=str,
                      required=True, help="The name of the output VCF file.")

  PARSER.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                      help='Increase verbosity of the output (substantially).'
                      + ' This option will include the calculated peptide'
                      + ' sequences in the logging output.')

  ARGS = PARSER.parse_args()

  if ARGS.verbose:
    LOGGER.setLevel(logging.INFO)

  VCFCONV = Vcf2OncodriveFML(infile  = ARGS.infile)

  VCFCONV.convert(outfile = ARGS.outfile)

  sys.stderr.write("Vcf conversion complete.\n")
