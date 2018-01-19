'''
Constants used in VCF parsing and annotation.
'''

import re

################################################################################

# Variant effects as constants; these terms are from the Sequence Ontology.
SYNM = 'synonymous_variant'      # SO:0001819
MISS = 'missense_variant'        # SO:0001583
FRSH = 'frameshift_variant'      # SO:0001589
SPAC = 'splice_acceptor_variant' # SO:0001574
SPDN = 'splice_donor_variant'    # SO:0001575
STGN = 'stop_gained'             # SO:0001587
STLS = 'stop_lost'               # SO:0001578
SRLS = 'start_lost'              # SO:0002012
IFIS = 'inframe_insertion'       # SO:0001821
IFDL = 'inframe_deletion'        # SO:0001822
INTV = 'intron_variant'          # SO:0001627
INTG = 'intergenic_variant'      # SO:0001628
SQVR = 'sequence_variant'        # SO:0001060 (fallback)

VEP_TAG = 'EFFECT'

# I guess the effects should probably be ordered in increasing
# severity thusly: synonymous (LOW), missense (MODERATE), inframe
# indel (MODERATE), splice site (HIGH), frameshift (HIGH), nonsense
# (HIGH)?
SEVERITY = (STGN, STLS, FRSH, SRLS, SPAC, SPDN,
            IFIS, IFDL, MISS, SYNM, INTV, INTG, SQVR)

# A regex matching null allele representations.
NULLPATT = re.compile('[-]')

################################################################################

# Strelka coding conventions.
STRELKA_INDELS = {'ref' : '0/0',
                  'het' : '0/1',
                  'hom' : '1/1'}

STRELKA_ALLELES = ('AU','CU','GU','TU')

