#!/usr/bin/env python3

# compute the hydrophobic moment for each peptide sequence in the input file
# input file must have one sequence per line

import sys
from modlamp.descriptors import PeptideDescriptor

input_fn = sys.argv[1]
for seq in open(input_fn).readlines():
    seq_clean = seq.strip()
    #Eisenberg hydrophobic moment
    desc = PeptideDescriptor(seq_clean,'eisenberg')
    desc.calculate_moment()
    y = desc.descriptor[0][0]
    print('%.3f' % y)
