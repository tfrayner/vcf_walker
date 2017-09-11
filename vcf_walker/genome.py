'''
Abstract base class which handles composition and delegation to the VcfGeneModel class.
'''

import sys
import gzip
import logging

from .vcf_meta import VcfMeta
from .gene_model import VcfGeneModel
from .utils import flexi_open

LOGGER = logging.getLogger()
LOGGER.addHandler(logging.StreamHandler(sys.stdout))
LOGGER.handlers[0].setFormatter(\
  logging.Formatter("[%(asctime)s]VEPANN_%(levelname)s: %(message)s"))
LOGGER.setLevel(logging.WARNING)

################################################################################

class NonDelegatableItem(AttributeError):
  pass

################################################################################

class VcfGenome(VcfMeta):
  '''
  Class used to manage VcfGeneModel instances.
  '''
  # These attributes are currently all delegated to self._gene_model
  _delegated_attrs = ('chromosomes', 'intervaldict',
                      'chrlen_cumsums', 'chridx')

  def __init__(self, picklefile, fasta, gtf, savefile, *args, **kwargs):

    super(VcfGenome, self).__init__(*args, **kwargs)

    setattr(self, '_gene_model', None)
    setattr(self, '_picklefile', picklefile)
    setattr(self, '_fasta',      fasta)
    setattr(self, '_gtf',        gtf)
    setattr(self, '_savefile',   savefile)

  def _read_gene_model(self):
    '''
    Build a new VcfGeneModel and pickle it to disk, or restore a
    previously pickled object.
    '''
    model = None

    # Restore pickled object.
    if self._picklefile is not None:
      try:
        from cPickle import load
      except ImportError:
        from pickle import load
      sys.stderr.write("Loading VcfGeneModel object from file...\n")
      with flexi_open(self._picklefile, 'rb') as pkfh:
        model = load(pkfh)

    # FASTA + GTF object initialisation.
    if self._fasta is not None:
      if self._gtf is None:
        raise StandardError("Must provide GTF file containing gene model (-g).")
      model = VcfGeneModel(fasta = self._fasta, gtf = self._gtf)
      if self._savefile is not None:
        try:
          from cPickle import dump
        except ImportError:
          from pickle import dump
        sys.stderr.write("Dumping VcfGeneModel object to file...\n")
        with gzip.open(self._savefile, 'wb') as pkfh:
          dump(model, pkfh, -1) # fh must be opened in binary mode.

    if model is None:
      raise StandardError("No FASTA+GTF or saved genome metadata files provided."
                          + " Unable to build in-memory gene model")

    setattr(self, '_gene_model', model)

  # The following adapted from
  # http://stackoverflow.com/a/21966266/5555544; see that post for
  # how to delegate methods rather than attributes.
  def __getattr__(self, attr_name):
    return self._delegate_attribute(attr_name)

  def _check_attr_name_is_delegated(self, attr_name):
    if attr_name not in self._delegated_attrs:
      raise NonDelegatableItem('{} can not be delegated'.format(attr_name))

  def _delegate_attribute(self, attr_name):
    '''
    Method caches a delegated attribute by linking the attribute
    to this class instance.
    '''
    self._check_attr_name_is_delegated(attr_name)

    # This is where the gene model is actually loaded, upon first
    # access to a delegated attribute.
    if self._gene_model is None:
      self._read_gene_model()

    attribute = getattr(self._gene_model, attr_name)
    setattr(self, attr_name, attribute)

    return attribute
