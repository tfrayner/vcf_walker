#!/usr/bin/env python

'''
Script to filter out variants from a VEP-annotated VCF based on
user-supplied criteria.
'''

import re
import vcf
from vcf_walker.vcf_meta import VcfMeta, LOGGER

CSQ_DELIM = '|'

class VcfFilter(VcfMeta):

  def filter(self, outfile, csqtag='CSQ',
             filters=('synonymous_variant',
                      'intergenic_variant',
                      'intron_variant')):

    # Filter out all variants with either no annotated effects, or
    # annotated effects which are entirely within the filter set.
    filters = list(filters) + ['']

    # Check that the CSQ tag is defined in the header; abort if not.
    if csqtag in self.reader.infos:
      colstr = re.sub('.*Format: ', '', self.reader.infos[csqtag].desc)
      csqcol = colstr.split(CSQ_DELIM).index('Consequence')
    else:
      raise StandardError("VCF header does not contain %s tag." % csqtag)
    
    with open(outfile, 'w') as out_fh:
      vcf_writer = vcf.Writer(out_fh, self.reader)

      for record in self.reader:
        effects = [ csq
                    for csqlist in record.INFO[csqtag]
                    for csq in csqlist.split(CSQ_DELIM)[csqcol].split('&') ]

        if any([ item not in filters for item in effects ]):
          vcf_writer.write_record(record)

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(description=\
                          "Filter variants from a VEP-annotated VCF based on user-supplied criteria.")

  PARSER.add_argument('-i', '--input-file', dest='infile', type=str,
                      required=True, help="The name of the input VCF file.")

  PARSER.add_argument('-o', '--output-file', dest='outfile', type=str,
                      required=True, help="The name of the output VCF file.")

  PARSER.add_argument('-f', '--filtered-variants', dest='filters', type=str,

                      # A list of variant consequence types which do
                      # not explicitly code for nonsynonymous changes,
                      # and which we generally want to remove from the
                      # output.
                      default=('synonymous_variant',
                               'intergenic_variant',
                               'intron_variant',
                               'upstream_gene_variant',
                               'downstream_gene_variant',
                               'non_coding_transcript_variant',
                               'non_coding_transcript_exon_variant',
                               'NMD_transcript_variant',
                               'stop_retained_variant',
                               'incomplete_terminal_codon_variant',
                               'coding_sequence_variant',
                               '3_prime_UTR_variant',
                               '5_prime_UTR_variant'),
                      nargs='+', help="The variant effect SO terms to remove from the output."
                      + " Note that variants affecting multiple transcript isoforms in"
                      + " different ways may retain some of these terms in the output.")

  ARGS = PARSER.parse_args()

  FILTER = VcfFilter(infile=ARGS.infile)

  FILTER.filter(ARGS.outfile, filters=ARGS.filters)
