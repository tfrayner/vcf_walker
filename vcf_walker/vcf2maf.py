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
from .constants import SEVERITY
from .vcf_meta import VcfMeta

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

# A quick check for self-consistency.
assert all( [ x in VEP_MAP for x in SEVERITY ] )

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

################################################################################

class Vcf2Tab(VcfMeta):
  '''
  Mixin class providing functionality to subclasses
  dedicated to converting VCF to various tabular formats.
  '''
  def get_dispatch_table(self):
    '''
    This method must be implemented in subclasses to return a dispatch
    table constructed as a list of 2-element tuples: (column_heading,
    function).
    '''
    raise NotImplementedError()

  def post_init_hook(self):
    '''
    This method is called prior to reading the first VCF record, to
    allow the record dispatch table to be modified to include information
    from the VCF header.
    '''
    pass

  def pre_record_hook(self, record):
    '''
    This method is called once for each VCF record, prior to any
    further calculations. It may be used to add annotation to a
    record prior to writing it out to disk.
    '''
    return record

  def _build_tabrow(self, record, altnums, sample):
    '''
    Build a list of table row values using the internal dispatch table.
    '''
    rowvals = [ str(self._dispatch_column(x[1], record,
                                          altnums, sample))
                for x in self.get_dispatch_table() ]
    return rowvals
  
  def _dispatch_column(self, disp, record, altnum, sample):
    '''
    Simple dispatch table interpreting method.
    '''
    if type(disp) == str:
      return disp
    else:
      return disp((record, altnum, sample))

  def convert(self, outfile):
    '''
    Read in a VCF and write out a MAF file.
    '''
    self.post_init_hook()

    warningcount = 0

    with open(outfile, 'w') as out_fh:

      # Header line.
      out_fh.write("\t".join([ x[0] for x in self.get_dispatch_table() ]) + "\n")
        
      for record in self.reader:
        with catch_warnings(record=True) as warnlist:
          record = self.pre_record_hook(record)

          # One row per sample:snv call. Multiallelic records need handling.
          for call in record.samples:
            altnums = self.call_snv_genotype(call, record)
            if altnums is not None and max(altnums) > 0:

              # This should handle 0/1, 0/2, 1/1, 1/2 etc. (altnum
              # of -1 will indicate REF).
              altnums = [ x - 1 for x in altnums ]
              rowvals = self._build_tabrow(record, altnums, call.sample)
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


################################################################################

class Vcf2Maf(VcfAnnotator, Vcf2Tab):
  '''
  Class designed to convert a VCF to MAF, pulling in variant effect
  annotations as necessary.

  MAF format specification taken from::

  https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification
  '''

  def __init__(self, center='cruk.cam.ac.uk', control='control', status='Somatic',
               source='WGS', sequencer='Illumina X10', *args, **kwargs):

    super(Vcf2Maf, self).__init__(*args, **kwargs)

    # TRANSCRIPT and VEP tags below are included by the call to
    # self.add_vep in the Vcf2Tab superclass.
    self._dispatch_table = [
      ('Hugo_Symbol',                    lambda (x): ";".join(x[0].INFO['TRANSCRIPT']) ),
      # NOTE this is not actually Entrez.
      ('Entrez_Gene_Id',                 lambda (x): ";".join(x[0].INFO['TRANSCRIPT']) ),
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

  def get_dispatch_table(self):
    return self._dispatch_table

  def post_init_hook(self):
    '''
    Retrieve genome assembly information from the VCF header. The
    exact location seems to vary with different PyVCF versions. We
    fall back to the 'reference' tag if contig information does not
    include the assembly.
    '''
    assembly = None

    if hasattr(self.reader, 'contigs'): # PyVCF 0.6.4 and above
      contigs = dict(self.reader.contigs).values()
      if hasattr(contigs[0], 'assembly'): # Not yet supported in PyVCF 0.6.8.
        assembs = list(set([ x.assembly.strip('"') for x in contigs ]))
        assembly = ";".join(assembs)

    if assembly is None: # Fall back to generic header tags.
      mdat = dict(self.reader.metadata)
      if 'contig' in mdat and 'assembly' in mdat['contig'][0]:
        assembs = list(set([ x['assembly'].strip('"') for x in mdat['contig'].values() ]))
        assembly = ";".join(assembs)

      elif 'reference' in mdat: # This is the most likely fallback at the moment.
        assembly = mdat['reference']
        
      else:
        assembly = 'unknown'

    # Modify the dispatch table on the fly.
    assert self._dispatch_table[3][0] == 'NCBI_Build'
    self._dispatch_table[3] = ('NCBI_Build', assembly)

    return

  def pre_record_hook(self, record):
    '''
    Add variant consequence (VEP) and affected TRANSCRIPT tags to the record.
    '''
    if any( [ x not in record.INFO for x in ('VEP','TRANSCRIPT') ] ):
      LOGGER.info("Inferring VEP tag for record %s:%s", record.CHROM, record.POS)
      record = self.add_vep(record)

    return record

################################################################################

class Vcf2OncodriveFML(Vcf2Tab):
  '''
  A somewhat overengineered class which simply converts VCF to the
  tab-delimited format expected by OncodriveFML.
  '''
  # Currently this class does not need to access the VcfGeneModel
  # annotation despite requiring it for initialisation. FIXME remove
  # this initialisation requirement, possibly by converting Vcf2Tab to
  # a mixin and moving annotation code into Vcf2Maf.

  def __init__(self, *args, **kwargs):

    super(Vcf2OncodriveFML, self).__init__(*args, **kwargs)

    self._dispatch_table = [
      ('CHROMOSOME', lambda (x): x[0].CHROM ),
      ('POSITION',   lambda (x): x[0].POS ),
      ('REF',        lambda (x): x[0].REF),
      # This line picks out the appropriate ALT allele. Note that this
      # is assuming Het throughout (ignores the first allele in the
      # genotype)
      ('ALT',        lambda (x): x[0].REF if x[1][1] < 0 else str(x[0].ALT[x[1][1]]) ),
      ('SAMPLE',     lambda (x): x[2] ),
    ]

  def get_dispatch_table(self):
    return self._dispatch_table
