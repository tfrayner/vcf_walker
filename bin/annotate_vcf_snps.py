#!/usr/bin/env python
#
# $Id: annotate_vcf_snps.py 4060 2017-07-14 15:04:32Z tfrayner $

'''
Script takes a VCF, a GTF gene model file containing CDS and exon
features, and the reference genome assembly fasta file, and uses these
to annotate SNPs in the former as to whether they are synonymous,
missense, nonsense, frameshift or splice site altering.
'''

import sys
import logging

from vcf_walker.annotator import VcfAnnotator, LOGGER

################################################################################

if __name__ == '__main__':

  from argparse import ArgumentParser

  PARSER = ArgumentParser(\
    description="Annotate variants in a VCF as missense, nonsense etc.")

  PARSER.add_argument('-i', '--input-file', dest='infile', type=str,
                      required=True, help="The name of the input VCF file.")

  PARSER.add_argument('-o', '--output-file', dest='outfile', type=str,
                      required=True, help="The name of the output VCF file.")

  GROUP = PARSER.add_mutually_exclusive_group(required=True)

  GROUP.add_argument('-f', '--fasta', dest='fasta', type=str,
                     help='The FASTA file containing the whole genome'
                     + ' sequence used to build the VcfAnnotator object')

  GROUP.add_argument('-p', '--pickle', dest='pickle', type=str,
                     help='The name of a pickle file containing the saved'
                     + ' VcfAnnotator object to reload.')

  PARSER.add_argument('-g', '--gtf', dest='gtf', type=str, required=False,
                      help='The name of the GTF file containing gene model'
                      + ' (CDS, exon) info. Only required if rebuilding a'
                      + ' VcfAnnotator object from scratch.')

  PARSER.add_argument('-s', '--savefile', dest='savefile', type=str,
                      required=False,
                      help='The name of the file to save a pickled VcfAnnotator'
                      + ' object to (for later reloading; this is much faster'
                      + ' than rebuilding from scratch every time).')

  PARSER.add_argument('-r', '--random', dest='random', action='store_true',
                      help='Indicates that the script should generate'
                      + ' randomised SNVs in the output. This may be used'
                      + ' for generating null distributions.')

  PARSER.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                      help='Increase verbosity of the output (substantially).'
                      + ' This option will include the calculated peptide'
                      + ' sequences in the logging output.')

  PARSER.add_argument('--with-gt', dest='genotypes', action='store_true',
                      help='Flag indicating that GT records should be inferred'
                      + ' (Strelka VCFs only, for the moment).')

  PARSER.add_argument('--with-vaf', dest='vaf', action='store_true',
                      help='Flag indicating that the variant allele frequency'
                      + ' should be calculated and stored for each locus'
                      + ' in each sample.')

  PARSER.add_argument('--without-effect', dest='noeffect', action='store_true',
                      help='Flag indicating that the inference of variant'
                      + ' effects on coding sequences should be skipped'
                      + ' (e.g. for performance reasons).')

  PARSER.add_argument('--without-context', dest='nocontext', action='store_true',
                      help='Flag indicating that the nucleotide context'
                      + ' need not be included in the output'
                      + ' (e.g. for performance reasons).')

  PARSER.add_argument('--filter-germline', dest='germline', action='store_true',
                      help='Flag indicating that germline contamination due to'
                      + ' samples sharing tumour sources should be detected and flagged.')

  PARSER.add_argument('--source-regex', dest='src_regex', type=str, default=r'_(\d{5})_',
                      help='Regex used to pull out source ID from the sample IDs'
                      + ' stored in the VCF. The first capture group should'
                      + ' correspond to the source ID. Used in conjunction'
                      + ' with the --filter-germline option.')

  PARSER.add_argument('--homozygous-normal', dest='homonorm', action='store_true',
                      help='Flag indicating that all samples labelled as NORMAL'
                      + ' should be assumed to be homozygous (used in the'
                      + ' context of the --withgt option above).')

  PARSER.add_argument('--context-bases', dest='cxtbases', type=int, default=3,
                      help='The number of context bases to include on'
                      + ' either side of the variant for the CONTEXT tag.')

  PARSER.add_argument('--vaf-cutoff', dest='vaf_cutoff', type=float, required=False,
                      help='The VAF threshold that each locus must achieve'
                      + ' before being included in the output. Note that'
                      + ' for loci mutated in multiple samples some calls'
                      + ' may still fall below this cutoff.')

  ARGS = PARSER.parse_args()

  if ARGS.verbose:
    LOGGER.setLevel(logging.INFO)

  VCFANN = VcfAnnotator(picklefile = ARGS.pickle,
                        fasta      = ARGS.fasta,
                        gtf        = ARGS.gtf,
                        savefile   = ARGS.savefile,
                        context_bases     = ARGS.cxtbases,
                        genotypes         = ARGS.genotypes,
                        vaf               = ARGS.vaf,
                        effect            = not ARGS.noeffect,
                        context           = not ARGS.nocontext,
                        filter_germline   = ARGS.germline,
                        source_regex      = ARGS.src_regex,
                        homozygous_normal = ARGS.homonorm)

  VCFANN.annotate_vcf(infile  = ARGS.infile,
                      outfile = ARGS.outfile,
                      random  = ARGS.random,
                      vaf_cutoff = ARGS.vaf_cutoff)

  sys.stderr.write("Vcf annotation complete.\n")
