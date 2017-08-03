'''
Functions used to manipulate CDS records in the context of the VCF
being processed.
'''

from copy import deepcopy
from warnings import warn

from Bio.SeqFeature import FeatureLocation

from .utils import AnnotationWarning
from .constants import SYNM, MISS, FRSH, SPAC, SPDN, STGN, STLS, \
  SRLS, IFIS, IFDL, INTV, INTG, SQVR, SEVERITY, NULLPATT

def check_splicesites(cds, record):
  '''
  Quick check on mutations lying within 2bp of the internal exon
  boundaries. Will also return an appropriate value if the SNV record
  falls in the larger body of the intron.
  '''
  # Splice site mutations are classed as any mutation in the two
  # nucleotides downstream of the donor (5') site or upstream of
  # the acceptor (3') site (see SO for definitions).
  #
  # I.e., mutations in the N residues here, where ^ represents the
  # actual RNA cleavage sites:
  #  ...ATTTACGAATCG^NNxxxxxx......xxxxxxxxNN^TACGTACGATACGTA...

  subfeats = get_subfeatures(cds)
  if subfeats is None:
    return None

  # Classic Schwartzian transform territory here (sorting so we know
  # which CDS subfeatures are at the start and the end, so that we
  # treat them appropriately).
  feats = [ tup[1] for tup in sorted([ (subfeat.location.start, subfeat)
                                       for subfeat in subfeats ],
                                     key=lambda tup: tup[0]) ]

  # Internal splice junctions for CDS are equivalent to exons, which
  # differ only in including 5' and 3' UTRs.
  donordists    = [ record.start - exon.location.end for exon in feats[:-1] ]
  acceptordists = [ exon.location.start - record.end for exon in feats[1:]  ]

  # Exactly which site is designated acceptor and which is donor
  # depends on the strand.
  if any([ dist <= 2 and dist > 0 for dist in donordists ]):
    return SPDN if cds.location.strand == 1 else SPAC
  elif any([ dist <= 2 and dist > 0 for dist in acceptordists ]):
    return SPAC if cds.location.strand == 1 else SPDN

  # If both distances are positive then the SNV falls in the intron
  # somewhere.
  elif any([ donordists[n] > 0 and acceptordists[n] > 0
             for n in range(len(donordists)) ]):
    return INTV

  return None

def decide_indel_or_nonsense(base_loss, pept_loss,
                             peptfirst, mutfirst):
  '''
  Simply figure out the likely variant effect given the observed
  loss of nucleotides and amino acid residues. Only applies to
  indels and nonsense stop_lost/stop_gained mutations. The starting
  amino acids are also compared to see if there is a start_lost
  mutation.
  '''
  if base_loss % 3 != 0: # Out-of-frame gain or loss of bases.

    # Frameshift. We assume the edge case where the ALT allele
    # introduces an immediate stop codon should still be categorised
    # as FRSH rather than STGN.
    return FRSH

  if pept_loss == 0:
    if peptfirst != mutfirst:
      return SRLS # Initial start (M) is no longer M. Start lost.
    else:
      return MISS # Must be otherwise missense in some way.

  assert pept_loss != 0

  if base_loss > 0: # In-frame loss of bases.

    if pept_loss < 0: # with gain of peptide; lost a stop codon.
      return STLS # stop_lost

    elif base_loss/3 == pept_loss: # base_loss matches pept_loss
      return IFDL # inframe_deletion

    else: # pept_loss > 0 but unexpectedly so.
      return STGN # stop_gained

  elif base_loss < 0: # In-frame gain of bases.

    if pept_loss > 0: # with loss of peptide; gained a stop codon.
      return STGN # stop_gained

    elif base_loss/3 == pept_loss: # base_loss matches pept_loss
      return IFIS # inframe_insertion

    else: # pept_loss < 0 but unexpectedly so.
      return STLS # stop_lost

  # base_loss == 0. Fallback is a simple nonsense mutation.
  if pept_loss > 0:
    return STGN # Gained a stop codon
  else:
    return STLS # Lost a stop codon

def get_cds_repr(cds):
  '''
  Return a string representation of the CDS, typically for logging
  purposes.
  '''
  if 'transcript_id' in cds.qualifiers:
    stub = ";".join(cds.qualifiers['transcript_id'])
  else:
    stub = get_cds_id(cds)

  return stub

def get_cds_id(cds):
  '''
  Attempt to pull out a meaningful ID from the CDS.
  '''
  if hasattr(cds, 'id'):
    stub = getattr(cds, 'id')
  else:
    stub = str(cds)

  if stub == '':
    if hasattr(cds, 'qualifiers'):
      qual = cds.qualifiers
      for item in ('gene_id', 'transcript_id', 'protein_id', 'gene_name', 'transcript_name', 'ccds_id'):
        if item in qual:
          stub = qual[item]
          if type(stub) is list:
            stub = ";".join(stub)
          break

  return stub

def trim_hanging_cds_end(seq):
  '''
  Add trailing Ns so that len(fullseq) is divisible by 3 (a
  biopython recommendation prior to translation).
  '''
  # Not currently used because if we fill in with Ns then there's a
  # tendency for biopython to somehow translate them out;
  # single-base frameshifts become missense.
  codon_length  = 3
  hanging = len(seq) % codon_length
  if hanging != 0:
    seq = seq[:-hanging]

  return seq

def spliced_cds_sequence(cds, chrom):
  '''
  We manually splice together exons because the standard SeqFeature
  translate method does not seem to do this reliably.
  '''
  subfeats = get_subfeatures(cds)
  if subfeats is None or len(subfeats) == 0:
    return trim_hanging_cds_end(cds.extract(chrom))

  subfeats = sort_subfeatures(cds)

  seqs = [ cds.extract(chrom) for cds in subfeats ]
  fullseq = seqs[0]
  if len(seqs) > 1:
    for nextseq in seqs[1:]:
      fullseq += nextseq

  return trim_hanging_cds_end(fullseq)

def sort_subfeatures(cds):
  '''
  Sort subfeatures by exon_number (if available) to ensure they are
  processed in the correct order.
  '''
  subfeats = get_subfeatures(cds)
  if 'exon_number' in subfeats[0].qualifiers:

    # Sort by exon number, test for consecutive numbers
    sortedtuples = sorted([ (int(cds.qualifiers['exon_number'][0]), cds)
                            for cds in subfeats ],
                          key=lambda tup: tup[0])
    exnums = [ tup[0] for tup in sortedtuples ]
    previous = exnums[0]
    if len(exnums) > 1:
      for num in range(1, len(exnums)):
        if exnums[num] - previous is not 1:
          # This is more common than I would like.
          stub = get_cds_repr(cds)
          warn("Nonconsecutive exon numbers found for %s." % stub,
               AnnotationWarning)
        previous = exnums[num]

    # All is presumably well
    subfeats = [ tup[1] for tup in sortedtuples ]

  else:
    warn("No exon numbers available; assuming provided"
         + " subfeature CDS order is consistent.",
         AnnotationWarning)

  return subfeats

def get_subfeatures(cds):
  '''
  Simple wrapper function to handle the alternative APIs provided by
  biopython <1.68 and >=1.68 to handle sub_features.
  '''
  # FIXME this is only really necessary because the BCBio.GFF parser
  # has not yet been properly fixed to use CompoundLocation with
  # biopython 1.68 and above. At some point the _sub_features
  # attribute is likely to vanish, at which point we will need to
  # recode this again.
  try:
    subfeats = cds.sub_features    # up to biopython 1.67
  except AttributeError:
    subfeats = cds._sub_features   # 1.68 and later (for now).

  return subfeats

def generate_mutant_cds(cds, record, altnum):
  '''
  Generate a mutated copy of the supplied CDS region based on indel
  info in the supplied variant record and alt allele number.
  '''
  base_loss = len(record.REF) - len(record.ALT[altnum])
  mutcopy = deepcopy(cds)
  if base_loss != 0:
    subfeats = get_subfeatures(mutcopy)
    if subfeats is not None:
      for subfeat in subfeats:
        substart = subfeat.location.start
        subend   = subfeat.location.end
        if subfeat.location.start > record.start:
          substart = substart - base_loss
        if subfeat.location.end   > record.start:
          subend   = subend   - base_loss
        subfeat.location = FeatureLocation(substart, subend,
                                           subfeat.location.strand)
    mutcopy.location = FeatureLocation(mutcopy.location.start,
                                       mutcopy.location.end - base_loss,
                                       mutcopy.location.strand)
  return mutcopy

