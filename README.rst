==========
vcf_walker
==========

This package provides various functions to extract and/or modify the
information stored in VCF files. Scripts are provided which can
extract and embed genome annotation for the variants in a VCF, or
convert between variant file formats. Some functionality assumes that
the VCFs conform to Strelka output conventions.

Installation
------------

As usual::

    python setup.py install

Python package prerequisites
----------------------------

   - PyVCF
   - BioPython
   - bx-python
   - bcbio-gff

Scripts
-------

annotate_vcf_snps.py
  Annotate variants (SNVs or Indels) using a reference genome fasta file alongside a GTF gene model file. The output 
  contains the affected transcript ID, predicted variant effect, immediate sequence context, and relationship to transcribed 
  strand for each variant. It is also possible to include the inferred variant allele frequency (VAF) and genotype (GT) 
  information for each variant in the output.
vcf_to_maf.py
  Converts variants in VCF format to MAF output.
vcf_to_oncodrivefml.py
  Converts variants in VCF format to the tab-delimited output format expected by OncodriveFML.

The utils directory also contains the following scripts not installed by default:

variant_gene_sequences.py
  Extract the CDS sequence and translation covering each variant in the supplied VCF, formatting it for readability. This
  script was originally designed to help with the design of validation PCR primers.
extract_cds_metadata.py
  This is a relatively simple utility which calculates the GC content and overall CDS length for each feature in a 
  VcfAnnotator object.
filter_vep_csq.py
  A very simple script to filter a VCF previously annotated by VEP using user-supplied criteria (e.g., a set of variant sequence
  ontology tags to ignore).

Credits
-------

© Tim Rayner, University of Cambridge, 2017
