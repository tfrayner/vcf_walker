'''
Top-level class providing basic metadata analysis for variants
described by a VCF.
'''

import sys
import gzip
import logging
from warnings import warn

import vcf

from . import __package_version__
from .utils import AnnotationWarning, order
from .constants import STRELKA_INDELS, STRELKA_ALLELES

LOGGER = logging.getLogger()
LOGGER.addHandler(logging.StreamHandler(sys.stdout))
LOGGER.handlers[0].setFormatter(\
  logging.Formatter("[%(asctime)s]VEPANN_%(levelname)s: %(message)s"))
LOGGER.setLevel(logging.WARNING)

################################################################################

class VcfMeta(object):
  '''
  Class used to handle basic metadata calculations for VCFs.
  '''
  def __init__(self, homozygous_normal=True):
    sys.stdout.write("# Initialising vcf_walker version %s #\n" % __package_version__)
    self.homozygous_normal = homozygous_normal

  def add_genotypes(self, record):
    '''
    Add in GT records on the fly; this is less memory-intensive than
    the addStrelkaGT code in the OdomLabLiverCancer R package.
    '''
    
    if 'GT' not in record.FORMAT.split(':'):
      record.FORMAT = 'GT:' + record.FORMAT 
    for call in record.samples:
      fields = call.data._asdict().keys()
      gt = self._infer_genotype(call, record)
      if 'GT' in fields:
        call.data = call.data._replace(GT=gt)
      else:
        CallData = vcf.model.collections.namedtuple('CallData', ['GT'] + fields)
        call.data = CallData(GT=gt, **call.data._asdict())
    
    return record

  def _infer_genotype(self, call, record):
    '''
    Infer our genotype call (currently assumes Strelka SNV input).
    '''

    # First, handle null calls in merged VCFs.
    if call.data.DP is None or call.data.DP == 0:
      return None

    if record.is_indel:

      # Now decipher the record SGT INFO field.  Assumes Strelka indel
      # coding. FIXME using SGT isn't ideal here as it won't
      # necessarily work for multi-sample VCFs.
      [ normal, tumor ] = record.INFO['SGT'].split('->')
      if call.sample == 'NORMAL':
        gt_all = normal
      elif call.sample == 'TUMOR':
        gt_all = tumor
      else:
        warn("Sample label is neither NORMAL nor TUMOR; assuming TUMOR.", AnnotationWarning)
        gt_all = tumor

      return STRELKA_INDELS[gt_all]
    
    elif record.is_snp:

      allindices = self.call_snv_genotype(call, record)

      return "/".join([ str(y) for y in sorted(allindices) ])

    else:
      raise ValueError("Record is neither SNV nor indel, and so is unsupported for genotype calling.")

  def call_snv_genotype(self, call, record):
    '''
    Code the genotype call for an SNV as a two-element list of
    integers: 0 for record.REF, 1 and above for record.ALT alleles.
    '''
    if not record.is_snp:
      raise StandardError("Non-SNV variant passed to method expecting SNVs only.")

    # First, handle null calls in merged VCFs.
    if call.data.DP is None or call.data.DP == 0:
      return None

    # Assumes this is a Strelka SNV.
    tier1  = [ getattr(call.data, att)[0] for att in STRELKA_ALLELES ]      # count tier1 reads
    ord    = order(tier1)

    allowable = record.alleles

    # Identify the top allele candidates in descending order;
    # requires at least two tier1 reads to be detected for any given
    # allele. This is not particularly stringent.
    gt_all    = [ STRELKA_ALLELES[n][0] for n in ord[::-1] if tier1[n] >= 2 ]

    # Filter out invalid alleles.
    gt_all    = [ base for base in gt_all if base in allowable ]

    # Assumes all normals are homozygous, if that's the desired behaviour.
    if len(gt_all) == 1 or (self.homozygous_normal and call.sample == 'NORMAL'):
      gt_all = [ gt_all[0], gt_all[0] ]

    # Trim off any spare alleles to take just the top two.
    gt_all    = gt_all[:2]

    # Convert to [ 0, 1 ] etc.
    allindices = [ allowable.index(x) for x in gt_all ]

    return allindices

  def calculate_vaf(self, call, record, altnum=None):
    '''
    Calculate the frequency for variant alleles at the given
    locus. Multiple ALT alleles at a given locus are summed together
    for each call.
    '''
    # First, handle null calls in merged VCFs.
    if call.data.DP is None or call.data.DP == 0:
      return None

    # Indels are not going to work here due to differences in Strelka coding.
    if not record.is_snp:
      raise StandardError("VAF calculation for non-SNV records is unsupported.")

    # Assumes this is a Strelka SNV.
    tier2 = [ getattr(call.data, att)[1] for att in STRELKA_ALLELES ] # count tier2 reads

    # Find read counts for the allowable allele calls; REF is first in list.
    ref = str(record.REF) + 'U'
    if altnum is None:

      # Identify the main ALT allele by tier2 read count.
      allowable = [ str(n) + 'U' for n in record.alleles ]
      allcounts = [ tier2[ STRELKA_ALLELES.index(n) ] for n in allowable ]
      counts    = [ allcounts[0], max(allcounts[1:]) ]

    else:

      # Use the specified ALT allele.
      allowable = [ ref, str(record.ALT[altnum]) + 'U' ]
      counts    = [ tier2[ STRELKA_ALLELES.index(n) ] for n in allowable ]

    # Assumes all normals are homozygous, if that's the desired behaviour.
    if len(counts) == 1 or (self.homozygous_normal and call.sample == 'NORMAL'):
      counts = [ counts[0], 0 ]

    # All selected ALT alleles summed, discarding highly speculative
    # alleles for which only a single read was detected.
    vaf = round(float(sum([ x for x in counts[1:] if x >= 2])) / sum(counts), 3)

    return vaf

  def add_vaf(self, record):
    '''
    Add Variant Allele Frequency (VAF) calculation to VCF, where each
    VAF value is the sum of all tier2 ALT allele reads divided by the
    total tier2 read depth.
    '''
    if 'VAF' not in record.FORMAT.split(':'):
      record.FORMAT = 'VAF:' + record.FORMAT

    for call in record.samples:

      # Calculate VAF for the ALT allele with the greatest number of
      # tier2 reads.
      vaf = self.calculate_vaf(call, record)

      fields = call.data._asdict().keys()
      if 'VAF' in fields:
        call.data = call.data._replace(VAF=vaf)
      else:
        CallData = vcf.model.collections.namedtuple('CallData', ['VAF'] + fields)
        call.data = CallData(VAF=vaf, **call.data._asdict())

    return record

