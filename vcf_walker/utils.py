'''
Basic utility functions used in VCF processing.
'''
import gzip
from contextlib import contextmanager
from itertools import izip

################################################################################

class AnnotationWarning(UserWarning):
  '''
  Custom warning class used to tag a fairly verbose set of warnings
  which we can do little about.
  '''
  pass

class GenotypingError(Exception):
  '''
  Error encountered in genotyping samples. May indicate that the
  record does not contain the expected FORMAT tags.
  '''
  pass

class ReadDepthError(Exception):
  '''
  Operation attempted which requires of DP tags, which are absent.
  '''
  pass

################################################################################

@contextmanager
def flexi_open(filename, *args, **kwargs):
  '''
  Simple context manager function to seamlessly handle gzipped and
  uncompressed files.
  '''
  if is_zipped(filename):
    handle = gzip.open(filename, *args, **kwargs)
  else:
    handle = open(filename, *args, **kwargs)

  yield handle

  handle.close()

def is_zipped(filename):
  '''
  Test whether a file is zipped or not, based on the file magic
  number. Simplified version of the osqpipe.pipeline.utilities function.
  '''
  return bool(open(filename).read(2) == '\x1f\x8b') # gzipped file magic number.

def cumsum(itr):
  '''
  Simple cumsum implementation (see http://stackoverflow.com/a/9258634).
  '''
  total = 0
  for val in itr:
    total += val
    yield total

def order(x):
  '''
  Returns the order of each element in x as a list.
  '''
  # See http://code.activestate.com/recipes/491268-ordering-and-ranking-for-lists/
  L = len(x)
  rangeL = range(L)
  z = izip(x, rangeL)
  z = izip(z, rangeL) #avoid problems with duplicates.
  D = sorted(z)
  return [d[1] for d in D]
  
def rank(x):
  '''
  Returns the rankings of elements in x as a list. Does not handle ties or None values.
  '''
  # See http://code.activestate.com/recipes/491268-ordering-and-ranking-for-lists/
  L = len(x)
  ordering = order(x)
  ranks = [0] * len(x)
  for i in range(L):
    ranks[ordering[i]] = i
  return ranks

# The following from https://stackoverflow.com/a/27758326 and edited
# to switch to sample STDEV and to only calculate the mean once.
def mean(data):
    """
    Return the sample arithmetic mean of data.
    """
    n = len(data)
    if n < 1:
        raise ValueError('mean requires at least one data point')
    return sum(data)/float(n)

def mean_and_sstdev(data):
    """
    Calculates the mean and sample standard deviation.
    """
    n = len(data)
    if n < 2:
        raise ValueError('variance requires at least two data points')
    mu = mean(data)
    ss = sum((x-mu)**2 for x in data)
    pvar = ss/(n-1) # the sample variance
    return (mu, pvar**0.5)
