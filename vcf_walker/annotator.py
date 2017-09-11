'''
Classes and functions used to add variant effects and mutational
context sequence tags to a VCF. Strelka format conventions are assumed
throughout.
'''

import sys
import re
import gzip
from copy import deepcopy
from warnings import warn, showwarning, catch_warnings
from random import uniform, choice

from Bio.Seq import Seq, MutableSeq
from Bio.SeqFeature import FeatureLocation

import vcf
from vcf.model import _Substitution as VcfSubstitution
from vcf.model import _Record       as VcfRecord
from vcf.model import _Call         as VcfCall

from math import sqrt

from .genome import VcfGenome, LOGGER
from .utils import AnnotationWarning, is_zipped, cumsum, order, mean_and_sstdev
from .cds import check_splicesites, decide_indel_or_nonsense, \
  get_cds_repr, trim_hanging_cds_end, spliced_cds_sequence, \
  generate_mutant_cds, get_subfeatures
from .constants import SYNM, MISS, FRSH, SPAC, SPDN, STGN, STLS, \
  SRLS, IFIS, IFDL, INTV, INTG, SQVR, SEVERITY, NULLPATT, STRELKA_INDELS, STRELKA_ALLELES

################################################################################

class VcfAnnotator(VcfGenome):
  '''
  Class designed to annotate variants provided in VCF with some
  information on the variant effects on protein and transcript structure.
  '''
  def __init__(self, context_bases=3,
               vaf=False, genotypes=False,
               effect=True, context=True,
               filter_germline=False, source_regex=r'_(\d{5})_',
               *args, **kwargs):

    super(VcfAnnotator, self).__init__(*args, **kwargs)

    self.context_bases     = context_bases

    # filter_germline uses genotype information, so we just activate
    # the genotypes option here.
    self.genotypes         = genotypes or filter_germline
    self.vaf               = vaf
    self.effect            = effect
    self.context           = context
    self.filter_germline   = filter_germline
    self.source_regex      = re.compile(source_regex)

  def generate_ref_objects(self, cds, record, utrlen=0):
    '''
    Generates a copy of the CDS and record objects alongside a
    reference sequence objects, with internal coordinates rewritten to
    minimise the in-memory size of the latter.
    '''
    cdscopy = deepcopy(cds)
    reccopy = deepcopy(record)

    # If the record lies entirely within the CDS region (as it will in
    # the vast majority of cases), we can shorten the sequence we
    # analyse to improve performance.
    cdsmin = min(cds.location.start-utrlen, cds.location.end-utrlen)
    cdsmax = max(cds.location.start+utrlen, cds.location.end+utrlen)
    if all([ recpos < cdsmax and recpos > cdsmin
             for recpos in (record.start, record.end) ]):

      # We want to work on just the chromosome sequence within the
      # bounds of the top-level CDS. First we remap the CDS FeatureLocations.
      genstart = cds.location.start-utrlen
      subfeats = get_subfeatures(cdscopy)
      if subfeats is not None:
        for subfeat in subfeats:
          subfeat.location = FeatureLocation(subfeat.location.start - genstart,
                                             subfeat.location.end   - genstart,
                                             subfeat.location.strand)
      cdscopy.location = FeatureLocation(cdscopy.location.start - genstart,
                                         cdscopy.location.end   - genstart,
                                         cdscopy.location.strand)

      # Remap the record location appropriately
      reccopy.start = reccopy.start - genstart
      reccopy.end   = reccopy.end   - genstart

      # Pull out the appropriate chromosome regions
      baseseq = self.chromosomes[reccopy.CHROM]\
          [(cds.location.start-utrlen):(cds.location.end+utrlen)].seq.tostring()
    else:

      # Overhanging variant (probably an indel). This treatment is
      # slower, but will handle overhanging indels etc.
      warn("Detected variant overhanging CDS end: %s:%d-%d"
           % (record.CHROM, record.start, record.end))
      baseseq = self.chromosomes[reccopy.CHROM].seq.tostring()

    return (cdscopy, reccopy, baseseq)

  def determine_variant_effect(self, cds, record, altnum=0):
    '''
    Given a top-level CDS object (inferred_parent) alongside a VCF
    record and the mutated chromosome (precomputed here for the sake
    of performance), figure out what the effect of the mutation
    is. In cases where multiple ALT alleles are described in the
    input VCF, the altnum argument is used to disambiguate.
    '''

    stub = get_cds_repr(cds)

    # First, check intron and splice-site variants because we don't
    # need much in the way of processing to figure those out.

    # FIXME also consider adding {3,5}_prime_UTR_variant (although
    # this is trickier as we do not store UTR ranges at the moment).

    # N.B. this also detects intronic variants INTV. I'm pretty sure
    # this should indeed be cds/cdscopy rather than mutcopy, as we're
    # not looking at downstream CDS features here.
    splicemut = check_splicesites(cds, record)
    if splicemut is not None:
      LOGGER.info('%s: %s', splicemut, stub)
      return splicemut

    # Heavy lifting begins here.
    (cdscopy, reccopy, baseseq) = self.generate_ref_objects(cds, record)

    origchrom = Seq(baseseq)
    mutchrom  = MutableSeq(baseseq)

    reference_allele = str(mutchrom[reccopy.start:reccopy.end])
    if reference_allele.lower() != reccopy.REF.lower():
      raise ValueError(("REF allele does not agree with reference: %s:%d-%d"
                        + " (REF %s, reference %s).")
                       % (record.CHROM, record.start, record.end,
                          reccopy.REF.upper(), reference_allele.upper()))

    mutchrom[reccopy.start:reccopy.end] = str(reccopy.ALT[altnum])

    # Indels need handling by rewriting the CDS model on the fly.
    mutcopy = generate_mutant_cds(cdscopy, reccopy, altnum)

    # First, establish whether there's any difference in the
    # output peptide. Then we'll figure out what that difference is.
    cdsseq  = spliced_cds_sequence(cdscopy, origchrom)
    peptide = cdsseq.translate(to_stop=True)

    mutseq  = spliced_cds_sequence(mutcopy, mutchrom)
    mutpep  = mutseq.translate(to_stop=True)

    if str(peptide) == str(mutpep):
      LOGGER.info('%s: %s', SYNM, stub)
      return SYNM # Synonymous

    base_loss = len(reccopy.REF) - len(reccopy.ALT[altnum])
    assert len(origchrom) - len(mutchrom) == base_loss
    pept_loss = len(peptide) - len(mutpep)

    # This is a quick check for start lost due to ATG disruption at
    # start of sequence (which results in an IndexError when accessing
    # mutpep[0]).
    if len(mutpep) <= 0:
      LOGGER.info('%s: %s %s -> %s', SRLS, stub, peptide, mutpep)
      return SRLS

    # This really shouldn't happen, and is probably indicative of gene
    # model inconsistencies.
    if len(peptide) <= 0:
      warn('Zero- or negative-length input peptide detected.')
      LOGGER.warning('Zero- or negative-length input peptide detected; returning %s.', SQVR)
      return SQVR

    vep = decide_indel_or_nonsense(base_loss, pept_loss, peptide[0], mutpep[0])

    LOGGER.info('%s: %s %s -> %s', vep, stub, peptide, mutpep)

    return vep

  def effects_by_affected(self, affected, record, altnum=0):
    '''
    Given a list of affected CDS records, and a VCF SNV record,
    identify the relative severity (0=worst) of the effect of the
    variant in each CDS record.  For VCF records with multiple ALT
    alleles the altnum argument can be used to disambiguate.
    '''
    if len(affected) > 0:
      effects = [ self.determine_variant_effect(cds.value['cdsobj'],
                                                record, altnum)
                  for cds in affected ]
    else:

      # Assume intergenic only.
      effects = [ INTG ]

    return [ SEVERITY.index(efct) for efct in effects ]

  def _find_worst_effects(self, affected, record, altnum=0):
    '''
    Given a list of affected CDS records, and a VCF SNV record,
    identify the worst effect that the SNV has on any of the CDS
    features. For VCF records with multiple ALT alleles the altnum
    argument can be used to disambiguate.
    '''
    effect_ranks = self.effects_by_affected(affected, record, altnum)

    worst = SEVERITY[ min(effect_ranks) ]

    # Identify transcript ID.
    if len(affected) > 0:
      widx = min(xrange(len(effect_ranks)), key=effect_ranks.__getitem__)
      txid = get_cds_repr(affected[ widx ].value['cdsobj'])
    else:
      txid = ''

    return (worst, txid)

  def clean_record_ends(self, record):
    '''
    Make sure our variant ends are not offset by null chars ('-' and
    '.' so far). Also rewrite the REF and ALT allele representations
    in our record copy.
    '''
    if record.CHROM not in self.chromosomes:
      raise IndexError(\
        "VCF contains chromosome not found in fasta reference: %s"
        % record.CHROM)
    record.end = record.end - len(NULLPATT.findall(record.REF))
    record.REF = NULLPATT.sub('', record.REF)
    record.ALT = [ VcfSubstitution( NULLPATT.sub('', str(allele)) )
                   for allele in record.ALT ]

    return record

  def add_vep(self, record):
    '''
    Method adds EFFECT, TRANSCRIPT and TXNSTR INFO tags to the record according to
    its effects on overlapping CDS regions.
    '''
    if record.CHROM not in self.chromosomes:
      raise IndexError(\
        "VCF contains chromosome not found in fasta reference: %s"
        % record.CHROM)

    record = self.clean_record_ends(record)

    # Annotate variant effect.
    affected = self.intervaldict[ record.CHROM ].find(record.start,
                                                      record.end)

    LOGGER.info("Examining %s:%s (%s -> %s)",
                record.CHROM, record.POS, record.REF,
                "".join([ str(x) for x in record.ALT]))
    worst_by_alt = [ self._find_worst_effects(affected, record, altnum)
                     for altnum in range(len(record.ALT)) ]

    record.INFO['EFFECT']     = [ x[0] for x in worst_by_alt ]
    record.INFO['TRANSCRIPT'] = [ x[1] for x in worst_by_alt ]

    # Annotate transcription strand.
    strands = list(set([ cds.value['cdsobj'].location.strand
                         for cds in affected ]))
    if len(strands) == 0:
      record.INFO['TXNSTR'] = 'N'
    elif len(strands) == 1:
      record.INFO['TXNSTR'] = 'F' if strands[0] == 1 else 'R'
    else:
      record.INFO['TXNSTR'] = 'B'

    return record

  def add_context(self, record):
    '''
    Method adds CONTEXT tag to the record, based on
    the underlying genome sequence.
    '''
    cxtbases = self.context_bases
    if cxtbases == 0:
      return record

    if record.CHROM not in self.chromosomes:
      raise IndexError(\
        "VCF contains chromosome not found in fasta reference: %s"
        % record.CHROM)

    # Rewrite indel records where necessary to include a base of context.
    chrom   = self.chromosomes[ record.CHROM ]
    alleles = [ str(record.REF) ] + [ str(allele) for allele in record.ALT ]
    if any([ len(allele) == 0 for allele in alleles ]):
      if record.start > 1:
        extrabase    = chrom[record.start - 1].upper()
        record.REF   = extrabase + record.REF
        record.ALT   = [ VcfSubstitution(extrabase + str(allele))
                         for allele in record.ALT ]
        record.start = record.start - 1
        record.POS   = record.POS   - 1
        # Consider rewriting record.ID as well? - maybe better not to?
      else: # beginning of chromosome as special case.
        extrabase    = chrom[record.end + 1].upper()
        record.REF   = record.REF + extrabase
        record.ALT   = [ VcfSubstitution(str(allele) + extrabase)
                         for allele in record.ALT ]
        record.end   = record.end + 1

    # Annotate base context.
    context = list(chrom[ record.start - cxtbases:record.end + cxtbases ]\
                     .seq.tostring().upper())
    contend = cxtbases + record.end - record.start

    # Quick sanity check.
    if "".join(context[ cxtbases:contend ]) != record.REF.upper():
      raise ValueError(("Extracted context sequence %s disagrees"
                        + " with reference allele %s.")
                       % ("".join(context), record.REF))

    # Replace the REF allele with x characters.
    context[ cxtbases:contend ] = ['x'] * (contend - cxtbases)
    record.INFO['CONTEXT'] = "".join(context)

    return record

  def flag_germline_variants(self, record, threshold=0.05, somatic_likelihood=0.98):
    '''This method identifies sample/record signs of germline
    contamination and returns a list of such records, split into
    multiple source groups if necessary, annotated with an appropriate
    FILTER flag.
    
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

    The reality is that most germline contaminants will be detected as
    an ALT call in every sample from that mouse, but for low coverage
    regions we can also accept thresholds a bit lower than 100%.

    '''
    source_matches = [ self.source_regex.search(call.sample) for call in record.samples ]
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

        # The minimum read depth that is likely for uncalled samples
        # (based on the depth at called samples; this assumes bulk
        # coverage is equivalent across all samples, which for the
        # current project is a reasonable assumption). The len(depths) must be >= 2.
        depths = [ record.samples[n].data.DP for n in range(len(sources)) if hits[n] ]
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

        LOGGER.debug("Germline probability: %.2f", prob_germline)

        # For threshold = 0.05, this would correspond to one uncalled
        # sample with minimum likely read depth 5, or two such samples
        # with min likely read depth 3. More uncalled samples would
        # usually be interpreted as evidence that this is a somatic
        # variant.
        if prob_germline > threshold:

          # Deal with problem SNVs here.
          LOGGER.warning("Probable germline contaminant variant identified in source %s: %s", src, record)

          # First determine whether we even need to split the record.
          if sum([ not is_src[n] and record.samples[n].data.GT is not None
                   for n in range(len(sources)) ]) == 0:

            # No non-src hits; set the flag and move on.
            record.FILTER += [ 'GermlineVariant' ]
          
          else:

            # Otherwise, deep copy the record and modify the FORMAT
            # fields for record.samples appropriately. Also append the
            # src string to the new record.ID. Check that the
            # delimiter is valid in VCF.
            newrec = deepcopy(record)
            newrec.FILTER += [ 'GermlineVariant' ]
            newrec.ID += '|%s' % src

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

  def _annotate_record(self, record):
    '''
    Method annotates a given record (variant effects, and a number of
    context positions). The return value is a list of records, to
    support cases where a record has to be split.
    '''

    # At some point we'll need actual record IDs, which are not always
    # present.
    if record.ID is None:
      record.ID = ("%s:%s_%s/%s"
                   % (record.CHROM,
                      record.POS,
                      record.REF,
                      ",".join([ str(x) for x in record.ALT ])))
      
    if self.effect:
      record = self.add_vep(record)

    if self.context:
      record = self.add_context(record)

    if self.vaf:
      record = self.add_vaf(record)

    if self.genotypes:
      record = self.add_genotypes(record)

    if self.filter_germline:
      return self.flag_germline_variants(record)

    else:
      return [ record ]

  def _randomise_record(self, record):
    '''
    Given a record, randomise the genome location and variant allele.
    '''
    if len(record.REF) != 1:
      raise ValueError("Variant randomisation only supported for SNVs.")

    vloc = int(uniform(0, self.chrlen_cumsums[-1]))
    which = [     vloc > self.chrlen_cumsums[x]
              and vloc <= self.chrlen_cumsums[x+1]
              for x in range(len(self.chrlen_cumsums)) ]

    record.CHROM = self.chridx[which.index(True)]
    record.POS   = vloc - self.chrlen_cumsums[which.index(True)]
    record.start = record.POS - 1
    record.end   = record.start + len(record.REF)

    chrom      = self.chromosomes[ record.CHROM ]
    record.REF = chrom[ record.start:record.end ].seq.tostring().upper()

    bases = ('A','C','G','T')
    newalt = []
    for _allele in range(len(record.ALT)):
      if len(newalt) == 3:
        raise ValueError("Already selected three variant alleles;"
                         + " there are only four bases!")
      choices = set(bases) - set([ str(x) for x in newalt ] + [ record.REF ])
      chosen  = choice(list(choices))
      newalt += [ VcfSubstitution(chosen) ]

    record.ALT = newalt

    # Also fix record.alleles and ID for consistency. Note that there
    # may be other tags (e.g. SGT for strelka) which contain allele
    # info which we aren't fixing.
    record.alleles = [ record.REF ] + record.ALT
    record.ID = "%s:%d_%s/%s" % (record.CHROM, record.POS,
                                 record.REF, "".join([ str(x)
                                                       for x in record.ALT ]))

    return record

  def _initialise_vcf_reader(self, infile):

    # Gzip file support is deferred to PyVCF itself, but our gzip file
    # detection is a bit more robust.
    vcf_reader = vcf.Reader(filename=infile, compressed=is_zipped(infile))

    # Not 100% sure this is a supported part of the PyVCF model.
    Info = vcf.model.collections.namedtuple('Info',
                                            ['id', 'num', 'type', 'desc'])
    vcf_reader.infos['EFFECT'] = Info(id='EFFECT', num='A', type='String',
                                   desc='The predicted effect'
                                   + ' of the ALT allele.')
    vcf_reader.infos['TRANSCRIPT'] = Info(id='TRANSCRIPT', num='A', type='String',
                                          desc='The affected transcript ID.')
    vcf_reader.infos['CONTEXT'] = Info(id='CONTEXT', num=1, type='String',
                                       desc='The base context of the locus.')
    vcf_reader.infos['TXNSTR'] = Info(id='TXNSTR', num=1, type='String',
                                      desc='The direction of transcription at'
                                      + ' the locus, relative to the forward'
                                      + ' strand. F=Forward, R=Reverse,'
                                      + ' B=Both, N=Neither.')

    Format = vcf.model.collections.namedtuple('Format',
                                              ['id', 'num', 'type', 'desc'])
    if self.genotypes and 'GT' not in vcf_reader.formats:
      vcf_reader.formats['GT'] = Format(id='GT', num=1, type='String', desc='Inferred genotype call.')

    if self.vaf and 'VAF' not in vcf_reader.formats:
      vcf_reader.formats['VAF'] = Format(id='VAF', num=1, type='Float', desc='Variant Allele Frequency'
                                         + ' for the primary ALT allele of each sample.')

    if self.filter_germline and 'GermlineVariant' not in vcf_reader.filters:
      Filter = vcf.model.collections.namedtuple('Filter',
                                                ['id', 'desc'])
      vcf_reader.filters['GermlineVariant'] = Filter(id='GermlineVariant',
                                                     desc='Variant is likely to be a germline contaminant'
                                                     + ' due to shared tumour origin across multiple samples.')

    return vcf_reader

  def annotate_vcf(self, infile, outfile, random=False, vaf_cutoff=None):
    '''
    Read in a VCF, add annotation and write it back out again.
    '''
    warningcount = 0

    vcf_reader = self._initialise_vcf_reader(infile)

    with open(outfile, 'w') as out_fh:
      vcf_writer = vcf.Writer(out_fh, vcf_reader)

      for record in vcf_reader:

        if vaf_cutoff is not None:
          vafs = [ self.calculate_vaf(call, record) for call in record.samples ]
          fvafs = [ x for x in vafs if x is not None ]
          if max(fvafs) < vaf_cutoff:
            continue

        if random:
          record = self._randomise_record(record)

        with catch_warnings(record=True) as warnlist:
          newrecs = self._annotate_record(record)
          for rec in newrecs:
            vcf_writer.write_record(rec)

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

    if warningcount > 0:
      sys.stderr.write("Detected AnnotationWarnings for %d variants.\n"
                       % warningcount)

################################################################################

