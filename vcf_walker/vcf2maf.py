'''
Classes and functions used to add variant effects and mutational
context sequence tags to a VCF.
'''

import sys
import re
from warnings import warn, catch_warnings, showwarning

import vcf

from .annotator import VcfAnnotator, LOGGER
from .utils import AnnotationWarning, flexi_open
from .cds import get_cds_id

################################################################################

# Map SO ontology terms to those expected by the MAF specification
# (because reasons). This should be kept up to date with variant
# effect terms in the .constants module.
VEP_MAP = {
  'synonymous_variant'      : 'Silent',
  'missense_variant'        : 'Missense_Mutation',
  'frameshift_variant'      : 'Frame_Shift', # compromise (Frame_Shift_Del or Frame_Shift_Ins)
  'splice_acceptor_variant' : 'Splice_Site',
  'splice_donor_variant'    : 'Splice_Site',
  'stop_gained'             : 'Nonsense_Mutation',
  'stop_lost'               : 'Nonstop_Mutation',
  'start_lost'              : 'Translation_Start_Site',
  'inframe_insertion'       : 'In_Frame_Ins',
  'inframe_deletion'        : 'In_Frame_Del',
  'intron_variant'          : 'Intron',
  'intergenic_variant'      : 'IGR',
  'sequence_variant'        : '', # Nothing really fits as a catch-all.
}

def _detect_var_type(record):

  if len(record.REF) == 1:
    return 'SNP'
  elif len(record.REF) == 2:
    return 'DNP'
  elif len(record.REF) == 3:
    return 'TNP'
  elif len(record.REF) <= 4:
    return 'ONP'
  else:
    raise ValueError("Unexpected reference allele length.") # e.g. indels

class Vcf2Maf(VcfAnnotator):
  '''
  Class designed to convert a VCF to MAF, pulling in variant effect
  annotations as necessary.

  MAF format specification taken from::

  https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification
  '''

  def _dispatch_column(self, disp, record, altnum, sample):
    '''
    Simple dispatch table interpreting method.
    '''
    if type(disp) == str:
      return disp
    else:
      return disp((record, altnum, sample))

  def add_gene(self, record):
    '''
    Method adds GENE_ID and GENE_NAME INFO tags to the record
    according to its overlapping CDSs.
    '''
    if record.CHROM not in self.chromosomes:
      raise IndexError(\
        "VCF contains chromosome not found in fasta reference: %s"
        % record.CHROM)

    record = self.clean_record_ends(record)

    # Annotate variant effect.
    affected = self.intervaldict[ record.CHROM ].find(record.start,
                                                      record.end)

    # Annotate record.INFO['GENE_ID'] etc. This is far from ideal as
    # often we're only given an Ensembl transcript ID rather than
    # Entrez Gene ID.
    if len(affected) == 0:

      # No gene hit; probably intergenic. Values here as given by TCGA
      # spec for "no gene within 3kb".
      record.INFO['GENE_ID']   = '0'
      record.INFO['GENE_NAME'] = 'Unknown'

    else:
      
      ids = [ get_cds_id(x.value['cdsobj']) for x in affected ]
      uniqids = list(set(ids))
      if len(uniqids) == 1:

        # One gene hit; simplest case.
        gene_id = uniqids[0]
      else:

        # Multiple gene IDs hit; look for worst effect and pick that gene.
        LOGGER.info("Multiple genes for record %s:%s; picking worst-hit gene.", record.CHROM, record.POS)
        effect_ranks = self.effects_by_affected(affected, record)
        gene_id = ids[effect_ranks.index(min(effect_ranks))]

      record.INFO['GENE_ID']   = gene_id
      record.INFO['GENE_NAME'] = gene_id

    return record

  def convert(self, infile, outfile,
              center='cruk.cam.ac.uk', control='control', status='Somatic',
              source='WGS', sequencer='Illumina X10'):
    '''
    Read in a VCF and write out a MAF file.
    '''
    self._dispatch_table = [
      ('Hugo_Symbol',                    lambda (x): x[0].INFO['GENE_NAME'] ),
      # NOTE this is not actually Entrez.
      ('Entrez_Gene_Id',                 lambda (x): x[0].INFO['GENE_ID'] ),
      ('Center',                         center),
      # Will be set from Vcf header upon initial file open.
      ('NCBI_Build',                     None),
      ('Chromosome',                     lambda (x): x[0].CHROM ),
      ('Start_Position',                 lambda (x): x[0].POS ),
      ('End_Position',                   lambda (x): x[0].POS + len(x[0].REF) - 1 ),
      # The only MAF-supported strand type.
      ('Strand',                         '+'),
      # Here we limit the classification to just Tumor_Seq_Allele2
      ('Variant_Classification',         lambda (x): VEP_MAP[ x[0].INFO['VEP'][x[1][1]] ] ),
      ('Variant_Type',                   lambda (x): _detect_var_type(x[0]) ),
      ('Reference_Allele',               lambda (x): x[0].REF),
      # This will generally be the reference allele, except for 1/1, 1/2 etc. calls.
      ('Tumor_Seq_Allele1',              lambda (x): x[0].REF if x[1][0] < 0 else str(x[0].ALT[x[1][0]]) ),
      # Picks out the appropriate ALT allele.
      ('Tumor_Seq_Allele2',              lambda (x): x[0].REF if x[1][1] < 0 else str(x[0].ALT[x[1][1]]) ),
      ('dbSNP_RS',                       ''),
      ('dbSNP_Val_Status',               ''),
      ('Tumor_Sample_Barcode',           lambda (x): x[2] ),
      ('Matched_Norm_Sample_Barcode',    control),
      ('Match_Norm_Seq_Allele1',         ''),
      ('Match_Norm_Seq_Allele2',         ''),
      ('Tumor_Validation_Allele1',       ''),
      ('Tumor_Validation_Allele2',       ''),
      ('Match_Norm_Validation_Allele1',  ''),
      ('Match_Norm_Validation_Allele2',  ''),
      ('Verification_Status',            'Unknown'),
      ('Validation_Status',              'Untested'),
      ('Mutation_Status',                status),
      ('Sequencing_Phase',               ''),
      ('Sequence_Source',                source),
      ('Validation_Method',              'none'),
      ('Score',                          'NA'),
      ('BAM_File',                       'NA'),
      ('Sequencer',                      sequencer),
      ('Tumor_Sample_UUID',              lambda (x): x[2] ),
      ('Matched_Norm_Sample_UUID',       control),
    ]
    
    with flexi_open(infile, 'r') as in_fh:
      vcf_reader = vcf.Reader(in_fh)

      # Genome assembly information is stored in the VCF header.
      mdat = dict(vcf_reader.metadata)
      if 'contig' in mdat and 'assembly' in mdat['contig'][0]:
        assembs = list(set([ x['assembly'].strip('"') for x in mdat['contig'].values() ]))
        assembly = ";".join(assembs)
      else:
        assembly = 'unknown'

      # Modify the dispatch table on the fly.
      assert self._dispatch_table[3][0] == 'NCBI_Build'
      self._dispatch_table[3] = ('NCBI_Build', assembly)

      warningcount = 0

      with open(outfile, 'w') as out_fh:

        # Header line.
        out_fh.write("\t".join([ x[0] for x in self._dispatch_table ]) + "\n")
        
        for record in vcf_reader:
          with catch_warnings(record=True) as warnlist:
            if 'VEP' not in record.INFO:
              LOGGER.info("Inferring VEP tag for record %s:%s", record.CHROM, record.POS)
              record = self.add_vep(record)

            record = self.add_gene(record)

            # One row per sample:snv call. Multiallelic records need handling.
            for sampnum in range(len(record.samples)):
              call = record.samples[sampnum]
              if call.data.GT not in (None, '0/0'):

                # This should handle 0/1, 0/2, 1/1, 1/2 etc. (altnum
                # of -1 will indicate REF).
                altnums = [ int(x) - 1 for x in call.data.GT.split('/') ]
                rowvals = [ str(self._dispatch_column(x[1], record,
                                                      altnums, call.sample))
                            for x in self._dispatch_table ]
                out_fh.write("\t".join(rowvals) + "\n")

            # Count AnnotationWarnings, show all others.
            annwarns = 0
            for wrn in warnlist:
              if issubclass(wrn.category, AnnotationWarning):
                annwarns += 1
              else:
                showwarning(wrn.message, wrn.category,
                            wrn.filename, wrn.lineno,
                            wrn.file, wrn.line)

            if annwarns > 0:
              warningcount += 1

    # Finally, report on AnnotationWarnings (typically non-consecutive exons).
    if warningcount > 0:
      sys.stderr.write("Detected AnnotationWarnings for %d variants.\n"
                       % warningcount)
