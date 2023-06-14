#!/usr/bin/env python3

# Representation of amino acids as five-bit or three-bit patterns
# for filtering protein databases
# Avril Coghlan, Dónall A. Mac Dónaill, Nigel H. Buttimore
# Bioinformatics, Volume 17, Issue 8, August 2001, Pages 676–685,
# https://doi.org/10.1093/bioinformatics/17.8.676 *)

import sys

bits_of_aa = {}
# this encoding scheme only supports the 20 natural AAs
bits_of_aa['V'] = "11111"
bits_of_aa['W'] = "11100"
bits_of_aa['M'] = "11011"
bits_of_aa['G'] = "11010"
bits_of_aa['A'] = "11001"
bits_of_aa['P'] = "11000"
bits_of_aa['Q'] = "10111"
bits_of_aa['E'] = "10110"
bits_of_aa['D'] = "10011"
bits_of_aa['N'] = "10010"
bits_of_aa['T'] = "10001"
bits_of_aa['H'] = "10000"
bits_of_aa['F'] = "01101"
bits_of_aa['I'] = "01011"
bits_of_aa['C'] = "01010"
bits_of_aa['L'] = "01001"
bits_of_aa['S'] = "01000"
bits_of_aa['R'] = "00100"
bits_of_aa['Y'] = "00010"
bits_of_aa['K'] = "00001"
end_of_sequence = "00000" # right padding pattern

aa_of_bits = {}
aa_of_bits["11111"] = 'V'
aa_of_bits["11100"] = 'W'
aa_of_bits["11011"] = 'M'
aa_of_bits["11010"] = 'G'
aa_of_bits["11001"] = 'A'
aa_of_bits["11000"] = 'P'
aa_of_bits["10111"] = 'Q'
aa_of_bits["10110"] = 'E'
aa_of_bits["10011"] = 'D'
aa_of_bits["10010"] = 'N'
aa_of_bits["10001"] = 'T'
aa_of_bits["10000"] = 'H'
aa_of_bits["01101"] = 'F'
aa_of_bits["01011"] = 'I'
aa_of_bits["01010"] = 'C'
aa_of_bits["01001"] = 'L'
aa_of_bits["01000"] = 'S'
aa_of_bits["00100"] = 'R'
aa_of_bits["00010"] = 'Y'
aa_of_bits["00001"] = 'K'

def aa_to_bits(aa):
    try:
        return bits_of_aa[aa]
    except KeyError:
        print('pep_bits.py: aa_to_bits: unsupported AA: %s' % aa,
              file=sys.stderr)
        exit(1)

def bits_to_aa(bits):
    try:
        return aa_of_bits[bits]
    except KeyError:
        print('pep_bits.py: bits_to_aa: unsupported bit pattern: %s' % bits,
              file=sys.stderr)
        exit(1)

def bits_of_sequence(max_len, seq):
    n = len(seq)
    if n <= max_len:
        res = ""
        for aa in seq:
            # print(aa)
            bits = aa_to_bits(aa)
            # print(bits)
            res += bits
        rem = max_len - n
        for _i in range(rem):
            res += end_of_sequence
        assert(len(res) == max_len * 5)
        return res
    else:
        print('bits_of_sequence: ignored; too long: %s' % seq,
              file=sys.stderr)
        return None

def consume_5_chars(s):
    left = s[0:5]
    right = s[5:]
    return (left, right)

def sequence_of_bits(bitstring):
    n = len(bitstring)
    assert(n % 5 == 0)
    m = int(n / 5)
    rem = bitstring
    res = ""
    for i in range(m):
        aa_bits, rest = consume_5_chars(rem)
        rem = rest
        if aa_bits == end_of_sequence:
            break
        aa = bits_to_aa(aa_bits)
        res += aa
    return res

input_fn = ""
max_len = 24 # recommended
try:
    max_len = int(sys.argv[1])
    input_fn = sys.argv[2]
except (IndexError, ValueError):
    print('usage:\n./pep_bits.py <max_len:int> <sequences:filename>',
          file=sys.stderr)
    exit(1)

# regression tests
assert(bits_of_sequence(24, 'QVFTLIKGATQLIRKTLGEQ') == '101111111101101100010100101011000011101011001100011011101001010110010000001100010100111010101101011100000000000000000000')
assert(sequence_of_bits('101111111101101100010100101011000011101011001100011011101001010110010000001100010100111010101101011100000000000000000000') == 'QVFTLIKGATQLIRKTLGEQ')

# main
for line in open(input_fn).readlines():
    stripped = line.strip()
    split = stripped.split('\t')
    seq = split[0]
    mic = split[1]
    bitstring = bits_of_sequence(max_len, seq)
    if bitstring != None:
        print('%s,%s,%s' % (seq,mic,bitstring))
