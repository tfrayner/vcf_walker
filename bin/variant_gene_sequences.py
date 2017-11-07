#!/usr/bin/env python

'''
Script extracts the CDS sequence and translation covering each
variant in the supplied VCF.
'''

import os
from vcf_walker.gene_sequence import VcfGeneSequence


if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(\
    description="Report on the DNA and peptide sequence for CDS regions overlapping variants in a VCF file.")
  
  PARSER.add_argument('-i', '--input-file', dest='infile', type=str,
                      required=True, help="The name of the input VCF file.")

  PARSER.add_argument('-o', '--output-file', dest='outfile', type=str,
                      required=True, help="The name of the output sequence file.")

  PARSER.add_argument('-u', '--utr-length', dest='utrlen', type=int, default=0,
                      help="The length of the 3'- and 5'-UTR regions to dump."
                      + " Note that this assumes there are no introns within these UTRs.")

  PARSER.add_argument('-w', '--width', dest='width', type=int, default=50,
                      help="The printed width of sequence stanzas in the output file.")

  GROUP = PARSER.add_mutually_exclusive_group(required=True)

  GROUP.add_argument('-f', '--fasta', dest='fasta', type=str,
                     help='The FASTA file containing the whole genome'
                     + ' sequence used to build the VcfModel object')

  GROUP.add_argument('-p', '--pickle', dest='pickle', type=str,
                     help='The name of a pickle file containing the saved'
                     + ' VcfModel object to reload.')

  PARSER.add_argument('-g', '--gtf', dest='gtf', type=str, required=False,
                      help='The name of the GTF file containing gene model'
                      + ' (CDS, exon) info. Only required if rebuilding a'
                      + ' VcfModel object from scratch.')

  PARSER.add_argument('-s', '--savefile', dest='savefile', type=str,
                      required=False,
                      help='The name of the file to save a pickled VcfModel'
                      + ' object to (for later reloading; this is much faster'
                      + ' than rebuilding from scratch every time).')

  PARSER.add_argument('--preferred', dest='preferred', type=str,
                      required=False,
                      help='A comma-delimited string indicating a list of'
                      + ' preferred transcripts which will override the default'
                      + ' selection of "worst-affected transcript".')

  PARSER.add_argument('--with-dna-counter', dest='dnacounter', action='store_true',
                      help='Whether or not to include a DNA base count. Note that'
                      + ' this count will not correspond to genome coordinates and'
                      + ' is just a relative count at this point.')

  ARGS = PARSER.parse_args()

  if os.path.exists(ARGS.outfile):
    raise StandardError("Will not overwrite pre-existing output file %s"
                        % ARGS.outfile)

  PREF = ARGS.preferred.split(',') if ARGS.preferred is not None else None

  VCFOBJ = VcfGeneSequence(infile     = ARGS.infile,
                           picklefile = ARGS.pickle,
                           fasta      = ARGS.fasta,
                           gtf        = ARGS.gtf,
                           savefile   = ARGS.savefile,
                           preferred  = PREF)

  VCFOBJ.dump_gene_sequence(ARGS.outfile, ARGS.utrlen, ARGS.width, ARGS.dnacounter)
