'''
Core class used to link VCF annotation to the genome.
'''

import sys

from bx.intervals.intersection import Intersecter, Interval

from Bio import SeqIO

from BCBio import GFF

from . import __format_version__
from .utils import cumsum, flexi_open
from .cds import get_subfeatures

class VcfGeneModel(object):
  '''
  Class used to link VCF records to gene- and sequence-level annotation.
  '''
  __slots__ = ('chromosomes', 'intervaldict', 'fasta', 'gtf',
               'chrlen_cumsums', 'chridx')

  def __init__(self, fasta, gtf):

    self.fasta = fasta
    self.gtf   = gtf

    sys.stderr.write("Reading FASTA file...\n")
    with flexi_open(fasta, 'rU') as handle:
      chromosomes = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))

    self.intervaldict = dict()

    sys.stderr.write("Reading GTF file (this will take some time)...\n")
    limit_info = dict(gff_type = ('CDS',))
    with flexi_open(gtf, 'r') as handle:
      for rec in GFF.parse(handle, limit_info=limit_info,
                           base_dict=chromosomes):

        # Fix strand info.
        for feature in rec.features: # Each top-level CDS
          if hasattr(feature, 'strand'):
            if feature.strand is None:
              # An unfortunate effect of the GFF parser.
              # Check subfeatures, take the first defined strand.
              # Note: may be bad assumption in weird species (ciliates??).
              for subfeat in get_subfeatures(feature):
                if subfeat.strand is not None:
                  feature.strand = subfeat.strand
                  break

        # Now, chromosomes['X'] is a SeqRecord with features each of
        # which has a feature.extract method which can be used to
        # access the underlying DNA sequence. CDS features are nested
        # automatically. See also feature.qualifiers for a list of IDs
        # (gene name, ID etc) associated with it.
        chromosomes[rec.id] = rec

        # We need to create some interval trees to identify
        # affected CDS features etc.
        self._index_record_in_intervaldict(rec)

    self.chromosomes = chromosomes

    self._precompute_chrlens()

    sys.stderr.write("Object initialisation complete.\n")

  def _index_record_in_intervaldict(self, rec):
    '''
    Insert CDS intervals from a chromosome SeqRecord into our cached
    interval tree dict.
    '''
    intersector = self.intervaldict.setdefault(rec.id, Intersecter())
    for feature in rec.features: # Each top-level CDS
      intvl = dict(chrom=rec.id,
                   value={'cdsobj': feature})
      if hasattr(feature, 'strand'):
        intvl['strand'] = feature.strand # could still be None.
      locus = Interval(int(feature.location.start),
                       int(feature.location.end),
                       **intvl)
      intersector.add_interval(locus)

  def _precompute_chrlens(self):
    '''
    Method precomputes chromosome lengths and stores them for later
    use (typically when random-sampling variant loci from the genome).
    '''
    # Order unimportant; just keep it consistent
    self.chridx = self.chromosomes.keys()

    lens = [ 0 ] + [ len(self.chromosomes[key]) for key in self.chridx ]

    self.chrlen_cumsums = list(cumsum(lens))

  def __getstate__(self):
    # Included to support pickling, since the initial object creation
    # is quite expensive.
    return dict(chromosomes=self.chromosomes, version=__format_version__,
                fasta=self.fasta, gtf=self.gtf)

  def __setstate__(self, state):
    # Included to support pickling. Note that intervaldict needs to be
    # rebuilt on unpickling.
    if state['version'] != __format_version__:
      raise StandardError("Unable to reset state: version mismatch"
                          + " (expected: %s; found: %s)."
                          % (__format_version__, state['version']))
    self.fasta        = state['fasta']
    self.gtf          = state['gtf']
    self.chromosomes  = state['chromosomes']
    sys.stderr.write("Rebuilding genome interval index...\n")
    self.intervaldict = dict()
    for rec in self.chromosomes.values():
      self._index_record_in_intervaldict(rec)
    self._precompute_chrlens()
    sys.stderr.write("Object state restored.\n")

