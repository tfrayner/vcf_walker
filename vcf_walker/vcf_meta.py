'''
Top-level class providing basic metadata analysis for variants
described by a VCF.
'''

import sys
import re
import gzip
import logging
from warnings import warn
from copy import deepcopy
from math import sqrt

import vcf

from . import __package_version__
from .utils import AnnotationWarning, order, mean_and_sstdev, is_zipped, \
  GenotypingError, ReadDepthError
from .constants import STRELKA_INDELS, STRELKA_ALLELES

LOGGER = logging.getLogger(__name__)
LOGGER.addHandler(logging.StreamHandler(sys.stdout))
LOGGER.handlers[0].setFormatter(\
  logging.Formatter("[%(asctime)s]VEPANN_%(levelname)s: %(message)s"))
LOGGER.setLevel(logging.WARNING)

GERMVAR_TAG = 'GermlineVariant'

################################################################################

class VcfMeta(object):
  '''
  Class used to handle basic metadata calculations for VCFs.
  '''
  def __init__(self, infile,
               homozygous_normal=True,
               source_regex=r'_(\d{5})_',
               lenient_germline=False):
    sys.stdout.write("# Initialising vcf_walker version %s #\n" % __package_version__)
    self.homozygous_normal = homozygous_normal
    self.source_regex      = re.compile(source_regex)
    self.lenient_germline  = lenient_germline
    self.reader            = self._initialise_vcf_reader(infile)

  def _initialise_vcf_reader(self, infile):
    # Gzip file support is deferred to PyVCF itself, but our gzip file
    # detection is a bit more robust.
    return vcf.Reader(filename=infile, compressed=is_zipped(infile))

  def add_genotypes(self, record):
    '''
    Add in GT records on the fly; this is less memory-intensive than
    the addStrelkaGT code in the OdomLabLiverCancer R package.
    '''
    
    if 'GT' not in record.FORMAT.split(':'):
      record.FORMAT = 'GT:' + record.FORMAT 
    for call in record.samples:
      fields = call.data._asdict().keys()
      gt = self.infer_genotype(call, record)
      if 'GT' in fields:
        call.data = call.data._replace(GT=gt)
      else:
        CallData = vcf.model.collections.namedtuple('CallData', ['GT'] + fields)
        call.data = CallData(GT=gt, **call.data._asdict())
    
    return record

  def infer_genotype(self, call, record):
    '''
    Infer our genotype call (currently assumes Strelka SNV or indel
    input). May be overridden in subclasses, e.g. to add support for
    other input types such as Manta SVs.
    '''
    if 'DP' not in call.data._asdict().keys():
      raise GenotypingError("Strelka DP read depth field not found in call tags.")

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
    # requires at least one tier1 read to be detected for any given
    # allele. This is not at all stringent.
    gt_all    = [ STRELKA_ALLELES[n][0] for n in ord[::-1] if tier1[n] >= 1 ]

    # Filter out invalid alleles.
    gt_all    = [ base for base in gt_all if base in allowable ]

    # Assumes all normals are homozygous, if that's the desired behaviour.
    if len(gt_all) == 1 or (self.homozygous_normal and call.sample == 'NORMAL'):
      gt_all = [ gt_all[0], gt_all[0] ]

    # Trim off any spare alleles to take just the top two based on read depth.
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
    numcounts = sum(counts)
    if numcounts == 0:
      LOGGER.warning("Tier2 read count of zero detected for record %s, call %s.", record, call)
      vaf = 0
    else:
      vaf = round(float(sum([ x for x in counts[1:] if x >= 2])) / numcounts, 3)

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

  def calc_germline_probability(self, record, hits, num_srcs, somatic_likelihood=0.98):
    '''
    Calculate the probability that the nth call in the given record is
    germline-derived (as opposed to somatic).

    The algorithm is briefly described as follows:

    For each locus, calculate the likely minimum read depth available
    for any uncalled samples, and use that to calculate the
    theoretical likelihood of missing a germline mutation
    (i.e. VAF=0.5) for every uncalled sample. This likelihood is then
    multiplied by (1 - somp), where somp=the likelihood of all calls
    being somatic rather than germline (based on the specified
    somatic_likelihood). This is a fudge factor, included to loosen
    the requirements for classifying variants as somatic in cases
    where only 2 or 3 samples are taken from a given source. In
    practice, somatic_likelihood=0.98 will only allow the 2-sample
    case; increase to 0.985 to include the 3-sample case. Reduce to
    zero to remove the fudge entirely.
    '''
    num_vars = sum(hits)
    assert num_vars > 1
    assert num_srcs > 1

    if 'DP' not in record.FORMAT.split(':'):
      raise ReadDepthError("Unable to perform germline probability calculation without DP tag.")

    # The minimum read depth that is likely for uncalled samples
    # (based on the depth at called samples; this assumes bulk
    # coverage is equivalent across all samples, which for the
    # current project is a reasonable assumption). The len(depths) must be >= 2.
    depths = [ record.samples[n].data.DP for n in range(len(hits)) if hits[n] ]
    (mu, sd) = mean_and_sstdev(depths)
    lower_bound = mu - (sd / sqrt(num_vars)) # Mean - SEM
    lower_bound -= 1 # Strelka would need more than one read to call a variant.
    LOGGER.debug("Lower bound: %.2f", lower_bound)

    # The probability that the uncalled samples are actually false
    # negatives, and are instead germline variants with VAF of
    # 0.5.
    fnp = min( 1, (0.5 ** lower_bound) ** (num_srcs - num_vars) )

    # The probability that none of the called variants are actually
    # somatic. This is a fudge to allow variants called in 2-3
    # samples from a single source to contribute to the outcome.
    somp = 1 - (somatic_likelihood ** num_vars)

    # The probability that this source has germline contamination
    # (fnp, with the somp fudge).
    prob_germline = somp * fnp

    return prob_germline

  def flag_germline_variants(self, record, threshold=0.05):
    '''This method identifies sample/record signs of germline
    contamination and returns a list of such records, split into
    multiple source groups if necessary, annotated with an appropriate
    FILTER flag.
    
    The reality is that most germline contaminants will be detected as
    an ALT call in every sample from that mouse, but for low coverage
    regions we can also accept thresholds a bit lower than 100%.

    '''
    source_matches = [ self.source_regex.search(call.sample) for call in record.samples ]
    if any([ x is None for x in source_matches ]):
      raise StandardError("Source regex fails to match at least one sample.")
    sources        = [ mobj.group(1) for mobj in source_matches ]
    srcset         = list(set(sources))

    records = []

    formtags = record.FORMAT.split(':')
    nullform = dict( ( x, None ) for x in formtags )

    for src in srcset:
      is_src = [ x == src for x in sources ]
      num_srcs = sum(is_src)
      if num_srcs > 1: # Shortcut the heavy processing for the vast majority of variants
        hits = [ is_src[n] and record.samples[n].data.GT is not None for n in range(len(sources)) ]

        num_vars = sum(hits)
        if num_vars <= 1: # Some sources will simply have one or no variants at this locus.
          continue        # We shortcut the calculation in such cases.

        prob_germline = self.calc_germline_probability(record, hits, num_srcs)

        LOGGER.debug("Germline probability: %.2f", prob_germline)

        # For threshold = 0.05, this would correspond to one uncalled
        # sample with minimum likely read depth 5, or two such samples
        # with min likely read depth 3. More uncalled samples would
        # usually be interpreted as evidence that this is a somatic
        # variant.
        if prob_germline > threshold:

          # Deal with problem SNVs here.
          LOGGER.warning("Probable germline contaminant variant identified in source %s: %s", src, record)

          # First determine whether we even need to split the
          # record. For conservative filtering, we just mark all SNVs
          # at a dubious locus as germline; otherwise, we check that
          # there are no variants from out-of-source.
          if (not self.lenient_germline) or sum([ not is_src[n] and record.samples[n].data.GT is not None
                                            for n in range(len(sources)) ]) == 0:

            # No non-src hits, or conservative filtering; set the flag and move on.
            if GERMVAR_TAG not in record.FILTER:
              record.FILTER += [ GERMVAR_TAG ]
          
          else:

            # Otherwise, deep copy the record and modify the FORMAT
            # fields for record.samples appropriately. Also append the
            # src string to the new record.ID. Check that the
            # delimiter is valid in VCF. Note that the order of
            # operations here is deliberate; for lenient filtering we
            # allow records at a GermlineVariant locus to PASS if the
            # probability is low for that source.
            newrec         = deepcopy(record)
            if GERMVAR_TAG not in record.FILTER:
              newrec.FILTER += [ GERMVAR_TAG ]
            newrec.ID     += '|%s' % src

            for n in range(len(record.samples)):
              if is_src[n]:

                # Remove germline SNV metadata from original record
                record.samples[n].data = record.samples[n].data._replace(**nullform)

              else:

                # Remove possibly-good SNV metadata from the germline record
                newrec.samples[n].data = newrec.samples[n].data._replace(**nullform)

            records += [ newrec ]

    # No need to include the stripped record if there's nothing of interest remaining.
    if any([ record.samples[n].data.GT is not None for n in range(len(record.samples)) ]):
      records = [ record ] + records

    return records

