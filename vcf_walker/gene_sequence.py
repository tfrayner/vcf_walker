'''
Classes and functions used to extract sequence information for
genes underlying the variants in a VCF.
'''

import sys
import re
from warnings import warn, showwarning, catch_warnings
from collections import OrderedDict

import vcf

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .annotator import VcfAnnotator, LOGGER
from .utils import flexi_open, AnnotationWarning
from .cds import get_cds_id, sort_subfeatures

# This is designed simply to pull out the numeric part of an Ensembl
# ID; the idea is that smaller, earlier-assigned numbers will
# typically be of more interest/use to the biologist who are the
# target market for this module.
IDREGEX=re.compile('^\w+?(\d+)$')

def _numeric_interval_id(interval):
  return int(IDREGEX.sub('\\1', get_cds_id(interval.value['cdsobj'])))

################################################################################

class VcfGeneSequence(VcfAnnotator):
  '''
  Class designed to extract the underlying CDS DNA sequence and
  print it out in a form that's useful to bench biologists (e.g. for
  primer design when verifying variant calls).
  '''
  def __init__(self, preferred, *args, **kwargs):
    if preferred is None:
      preferred = []
    self.preferred_transcripts = preferred
    super(VcfAnnotator, self).__init__(*args, **kwargs)

  def retrieve_sequence(self, record, utrlen=0):
    '''
    Method retrieves a SeqRecord object for the gene worst affected by
    the given variant record.
    '''
    if record.CHROM not in self.chromosomes:
      raise IndexError(\
        "VCF contains chromosome not found in fasta reference: %s"
        % record.CHROM)

    # Annotate variant effect.
    affected = self.intervaldict[ record.CHROM ].find(record.start,
                                                      record.end)

    LOGGER.info("Examining %s:%s (%s -> %s)",
                record.CHROM, record.POS, record.REF,
                "".join([ str(x) for x in record.ALT]))

    # Find the transcripts affected by the variant.
    affected = self.intervaldict[ record.CHROM ].find(record.start,
                                                      record.end)

    if len(affected) == 0:
      return None
    else:

      # Unless manually overridden, identify the lowest-numbered,
      # worst-affected transcript.
      worst = None
      affected.sort(key=lambda x: _numeric_interval_id(x))

      for interval in affected:
        for pref in self.preferred_transcripts:
          if get_cds_id(interval.value['cdsobj']) == pref:
            worst = interval
            break

      if worst is None:
        effect_ranks = self.effects_by_affected(affected, record)
        worst = affected[effect_ranks.index(min(effect_ranks))]

      # These objects have correct relative coordinates but will
      # differ from the reference.
      (cdscopy, reccopy, baseseq) = self.generate_ref_objects(worst.value['cdsobj'], record,
                                                              utrlen=utrlen)
      origchrom = Seq(baseseq)
      cdsseqrec = SeqRecord(seq=Seq(baseseq),
                            features=[cdscopy],
                            id=get_cds_id(cdscopy))

      return cdsseqrec

  def _pad_utrs(self, seqstr, utrlen=0):
    '''
    Pads a given string or Seq object with spaces representing 3'- and 5'-UTRs.
    '''
    return (" " * utrlen) + seqstr + (" " * utrlen)

  def format_output_seqrec(self, seqrec, utrlen=0, width=50, with_dna_counter=False):
    '''
    Create a string representation of the given SeqRecord object.
    '''
    printstr = ""
    
    for feature in seqrec.features:
      featdna = feature.extract(seqrec) # seqrecord; .seq for seq only.
      if feature.type == 'CDS':
        pep = featdna.seq.translate(to_stop=True)

        # The exon ends will make slightly more intuitive sense if the
        # peptide residues are printed at the end of each codon.
        pep = Seq("".join([ "  " + x for x in pep ]))
        cdsmark = "-" * len(featdna.seq)

      elif feature.type == 'inferred_parent':
        (pep, cdsmark) = self._format_inferred_parent_feature(feature, seqrec)
        
      else:
        raise ValueError("Unexpected feature type: %s" % feature.type)

      gendna = seqrec.seq
      if feature.strand == -1:
        gendna = gendna.reverse_complement()

      pepmark = self._format_sequence_counter(pep, step=10)

      seqhash = OrderedDict()
      seqhash['Residue' ] = self._pad_utrs(pepmark, utrlen)
      seqhash['Peptide']  = self._pad_utrs(pep,     utrlen)
      seqhash['Exon']     = self._pad_utrs(cdsmark, utrlen)
      seqhash['DNA']      = gendna
      if with_dna_counter:
        seqhash['Base'] = self._format_sequence_counter(gendna, step=width)
      printstr += self.wrap_sequences(seqhash, width=width)

    return printstr

  def _format_sequence_counter(self, seq, step=10, start=1):
    '''
    Add a counter track to a given sequence (skipping blank residues).
    '''
    # This is perhaps a little inefficient, but very reusable code
    # nonetheless. Works on both DNA and peptide.
    marks = [" "] * len(seq)
    resnum = start - 1
    for posnum in range(len(seq)):
      if seq[posnum] != ' ':
        resnum += 1
        if resnum % step == 0:
          marker = str(resnum)
          offset = posnum - len(marker)
          marks[(offset+1):(posnum+1)] = marker

    return "".join(marks)
  
  def _format_inferred_parent_feature(self, feature, seqrec):

    # Peptide translation; includes spaces
    pep = [" "] * len(feature)

    # CDS/exon marker
    cdsmark = [" "] * len(feature)

    # Carry over untranslated bases from one exon to the next.
    prev_overhang = Seq("")

    for subfeat in sort_subfeatures(feature):

      if subfeat.type != 'CDS':
        raise ValueError("Unexpected subfeature type: %s" % subfeat.type)

      cdsseq = prev_overhang + subfeat.extract(seqrec)

      overbases = len(cdsseq) % 3
      overhang  = cdsseq[-overbases:].seq if overbases > 0 else Seq("")

      start    = subfeat.location.start - feature.location.start
      end      = subfeat.location.end   - feature.location.start

      if feature.strand == -1:
        tmp = start
        start = len(feature) - end
        end   = len(feature) - tmp

      pepstart = start - len(prev_overhang)
      pepend   = end - overbases

      # The exon ends will make slightly more intuitive sense if the
      # peptide residues are printed at the end of each codon. (x + 2, below).
      positions = [ x + 2 for x in range(pepstart, pepend, 3) ]
      cdspep    = str(cdsseq.seq.translate()) # silently discards overhang.
      for resnum in range(len(positions)):
        basenum = positions[resnum]
        pep[ basenum ] = cdspep[ resnum ]

      cdsmark[ start:end ] = ["-"] * ( end - start )

      prev_overhang = overhang

    pep     = Seq("".join(pep))
    cdsmark = "".join(cdsmark)
    return (pep, cdsmark)

  def wrap_sequences(self, seqhash, width=50):
    '''
    Given a hash of labelled strings or sequences (all having the same
    length), create a string with the hash contents interleaved, with
    stanzas of the given printed width.
    '''
    seqlens = list(set([ len(x) for x in seqhash.values() ]))
    if len(seqlens) > 1:
      raise ValueError("wrap_sequences seqhash argument must have values all the same length")

    sidebar = "%%%ds  " % max([ len(x) for x in seqhash.keys() ])

    printstr = ''
    maxlen = seqlens[0]
    for start in range(0, maxlen, width):
      end = start + width
      if end > maxlen:
        end = maxlen
      for (key, seq) in seqhash.iteritems():
        printstr += sidebar % key
        printstr += str(seq[start:end]) + "\n"
      printstr += "\n"

    return printstr

  def dump_gene_sequence(self, infile, outfile, utrlen=0, width=50, with_dna_counter=False):
    '''
    Read in a VCF and print out gene sequence annotation for the
    variants it contains. Typically this will be used on a tiny subset
    of variants that one desires to validate.
    '''
    warningcount = 0

    with flexi_open(infile, 'r') as in_fh:
      vcf_reader = vcf.Reader(in_fh)

      with open(outfile, 'w') as out_fh:

        for record in vcf_reader:

          with catch_warnings(record=True) as warnlist:
            cdsseqrec = self.retrieve_sequence(record, utrlen)

            if cdsseqrec is None:
              continue

            out_fh.write("############### %s ###############\n\n" % cdsseqrec.id)
            out_fh.write(self.format_output_seqrec(cdsseqrec, utrlen,
                                                   width=width, with_dna_counter=with_dna_counter))
            out_fh.write("\n\n")

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

