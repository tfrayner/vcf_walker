#!/usr/bin/env python
#
# $Id: vcf_to_maf.py 4006 2017-03-31 09:47:13Z tfrayner $

'''
Script converts a VCF to MAF, using a pickled VcfGeneModel object
to infer variant effects and affected gene IDs.
'''

import sys
import logging

import vcf
from vcf_walker.vcf2maf import Vcf2Maf, LOGGER, detect_maf_var_type
from vcf_walker.constants import VEP_TAG

################################################################################
class Vcf2LCEMaf(Vcf2Maf):

  def __init__(self, *args, **kwargs):
    super(Vcf2LCEMaf, self).__init__(*args, **kwargs)

    # called with (record, altnums, call) as arguments.
    self._dispatch_table = [
      (self.gene_id_type,                lambda (x): ";".join(sorted(list(set(x[0].INFO['GENE_ID'])))) ),
      (self.gene_symbol_type,            lambda (x): ";".join(sorted(list(set(x[0].INFO['GENE_NAME'])))) ),
      # Will be set from Vcf header upon initial file open.
      (self.assembly_type,               None),
      ('Chromosome',                     lambda (x): x[0].CHROM ),
      ('Start_Position',                 lambda (x): x[0].POS ),
      ('End_Position',                   lambda (x): x[0].POS + len(x[0].REF) - 1 ),
      ('Strand',                         '+'),
      ('Variant_Context',                lambda (x): x[0].INFO['CONTEXT']),
      ('Variant_Allele_Frequency',       lambda (x): x[2].data.VAF),
      ('Variant_Filter',                 lambda (x): 'PASS' if x[0].FILTER is None else ','.join(x[0].FILTER)),
      ('Variant_Classification',         lambda (x): x[0].INFO[VEP_TAG][x[1][1]] ),
      ('Variant_Type',                   lambda (x): detect_maf_var_type(x[0]) ),
      ('Reference_Allele',               lambda (x): x[0].REF),
      # This will generally be the reference allele, except for 1/1, 1/2 etc. calls.
      ('Tumor_Seq_Allele1',              lambda (x): x[0].REF if x[1][0] < 0 else str(x[0].ALT[x[1][0]]) ),
      # Picks out the appropriate ALT allele.
      ('Tumor_Seq_Allele2',              lambda (x): x[0].REF if x[1][1] < 0 else str(x[0].ALT[x[1][1]]) ),
      # These currently give potentially misleading output since the sort order here will differ from that for GENE_ID etc.
#      ('Peptide_Position',               lambda (x): ";".join(sorted(list(set(x[0].INFO['PEP_POS'])))) ),
#      ('Peptide_Variant',                lambda (x): ";".join(sorted(list(set(x[0].INFO['PEP_VAR'])))) ),
      ('Tumor_Sample_Barcode',           lambda (x): x[2].sample ),
    ]

  def pre_record_hook(self, record):
    '''
    Add variant consequence (VEP) and affected TRANSCRIPT tags to the record.
    '''
    if any( [ x not in record.INFO for x in ('CSQ', 'CONTEXT') ] ):
      raise StandardError("CSQ and/or CONTEXT tags not found for record; cannot proceed.")

    # Pull out some of the more useful information from the VEP CSQ field.
    csqdat = [ csq.split('|') for csq in record.INFO['CSQ'] ]
    record.INFO[VEP_TAG]     = [ csq[1] for csq in csqdat ]
    record.INFO['GENE_ID']   = [ csq[4] for csq in csqdat ]
    record.INFO['GENE_NAME'] = [ csq[3] for csq in csqdat ]
    record.INFO['PEP_POS']   = [ csq[14] for csq in csqdat ]
    record.INFO['PEP_VAR']   = [ csq[15] for csq in csqdat ]

    return record

################################################################################

class FilteredVcf2Maf(Vcf2LCEMaf):

  def __init__(self, vaf_cutoff=None, vcf_filter=None, *args, **kwargs):
    super(FilteredVcf2Maf, self).__init__(*args, **kwargs)

    self.vaf_cutoff = vaf_cutoff
    if vcf_filter is None:
      self.vcf_filter = None  # Need to set the attribute to prevent delegation.
    else:
      self.vcf_filter = vcf_filter.split(',')

  def filtered_record(self, record):
    '''
    Method takes a VCF record object and returns a (possibly modified)
    record to be output. If the record is to be filtered from the
    eventual output the method will return False.
    '''
    # PyVCF sets PASS record.FILTER to None; we accept all such records.
    if record.FILTER is not None and self.vcf_filter is not None:
      if all([ term not in record.FILTER for term in self.vcf_filter ]):
        return None

    if self.vaf_cutoff is not None:
      if 'VAF' not in record.FORMAT.split(':'):
        raise StandardError("Attempted to use VAF filter on a VCF lacking VAF FORMAT fields.")

      passcount = 0
      for call in record.samples:
        if call.data.VAF is not None:
          if call.data.VAF < self.vaf_cutoff:
            fields = call.data._asdict().keys()
            CallData = vcf.model.collections.namedtuple('CallData', fields)
            data = dict(zip(fields, [None for n in fields]))
            call.data = CallData(**data)
          else:
            passcount += 1

      if passcount == 0: # No calls with passing VAF.
        return None

    return record

  def pre_record_hook(self, record):
    '''
    Hook method detects records to be filtered and omits them from the
    output.
    '''
    record = super(FilteredVcf2Maf, self).pre_record_hook(record)
    if record is None:
      return None

    return self.filtered_record(record)

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

  PARSER.add_argument('--vaf-cutoff', dest='vafcutoff', type=float, required=False,
                      help='The lower bound cutoff of the minimum VAF'
                      + ' allowed for output records.')

  PARSER.add_argument('--vcf-filter', dest='vcffilter', type=str, required=False,
                      help='VCF FILTER values which are allowed in the output records'
                      + ' ("PASS" is always accepted). Multiple values should be separated by commas.')

  PARSER.add_argument('--gene-symbol-type', dest='symboltype', type=str, default='Hugo_Symbol',
                      help='The output heading to use for the MAF gene symbol column.')

  PARSER.add_argument('--gene-id-type', dest='geneidtype', type=str, default='Entrez_Gene_Id',
                      help='The output heading to use for the MAF gene ID column.')

  PARSER.add_argument('--assembly-type', dest='assemblytype', type=str, default='NCBI_Build',
                      help='The output heading to use for the genome assembly column.')

  PARSER.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                      help='Increase verbosity of the output (substantially).'
                      + ' This option will include the calculated peptide'
                      + ' sequences in the logging output.')

  ARGS = PARSER.parse_args()

  if ARGS.verbose:
    LOGGER.setLevel(logging.INFO)

  VCFCONV = FilteredVcf2Maf(infile     = ARGS.infile,
                            picklefile = ARGS.pickle,
                            fasta      = ARGS.fasta,
                            gtf        = ARGS.gtf,
                            savefile   = ARGS.savefile,
                            vaf_cutoff = ARGS.vafcutoff,
                            vcf_filter = ARGS.vcffilter,
                            gene_symbol_type = ARGS.symboltype,
                            gene_id_type     = ARGS.geneidtype,
                            assembly_type    = ARGS.assemblytype)

  VCFCONV.convert(outfile = ARGS.outfile)

  sys.stderr.write("Vcf conversion complete.\n")
