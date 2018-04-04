#!/usr/bin/env python

'''
Script to extract variant loci known to be mutated in one of the
indicated samples, and omit all other loci.
'''

import re
import vcf
from vcf_walker.vcf_meta import VcfMeta, LOGGER

class SampleExtractor(VcfMeta):

  def filter(self, outfile, samppatt):

    sampre = re.compile(samppatt)

    with open(outfile, 'w') as out_fh:
      vcf_writer = vcf.Writer(out_fh, self.reader)

      for record in self.reader:
        if any([ sampre.search(call.sample)
                 and call.data.GT not in ('.', None)
                 for call in record.samples ]):
          vcf_writer.write_record(record)

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(description=\
                          "Extract variants from a VCF based on sample names.")

  PARSER.add_argument('-i', '--input-file', dest='infile', type=str,
                      required=True, help="The name of the input VCF file.")

  PARSER.add_argument('-o', '--output-file', dest='outfile', type=str,
                      required=True, help="The name of the output VCF file.")

  PARSER.add_argument('-s', '--sample-pattern-file', dest='sampfile', type=str,
                      required=True, help="A file containing sample ID patterns,"
                      + " one per line (commented lines will be ignored).")

  ARGS = PARSER.parse_args()

  with open(ARGS.sampfile) as sfh:
    sample_ids = [ x.strip() for x in sfh if not re.match(r'#', x) ]

  # Build a general sample ID regex with field delimiters based on
  # what's in the supplied IDs.
  omit = ''
  if any([re.search(r'[a-zA-Z]', x) for x in sample_ids]):
    omit += 'a-zA-Z'
  if any([re.search(r'\d', x) for x in sample_ids]):
    omit += '0-9'
  for test in (r'_', r'-'):
    if any([re.search(test, x) for x in sample_ids]):
      omit += test

  # Wrap each sample ID in negative look-ahead and -behind assertions
  # to avoid things like do888 matching do8888.
  samppatt = r'(' + r'|'.join([ r'(?<![' + omit + r'])' + idstr + r'(?![' + omit + r'])'
                                for idstr in sample_ids ]) + r')'

  FILTER = SampleExtractor(infile=ARGS.infile)

  FILTER.filter(ARGS.outfile, samppatt=samppatt)
