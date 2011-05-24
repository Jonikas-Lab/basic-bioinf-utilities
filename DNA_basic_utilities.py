#! /usr/bin/env python

"""
Various basic DNA utilities I wrote (parse_fasta, complement, reverse_complement, find_seq_between, generate_seq_slices, ...) - most of them just imported from their own files, since I need them as separate programs too.
Weronika Patena, 2008-2010
"""

# help_functions.
from parse_fasta import parse_fasta,print_seq
from read_input import read_input
from transform_sequence_input import transform_sequence_input
# complement a sequence:
from complement import complement
# reverse-complement a sequence:
from reverse_complement import reverse_complement
# reverse a sequence:
from reverse import reverse

### some help functions I'm just going to write here... ###

def find_seq_between(seq, flankA, flankB, exclusive=True):
    """ Return the fragment of seq between flankA and flankB (first occurence). """
    posA = seq.find(flankA)
    posB = seq[posA:].find(flankB)  # only look for flankB in the region AFTER flankA
    if not posB==-1:    posB += posA        # if found, add the length of the omitted region to result 
    lenA = len(flankA)
    lenB = len(flankB)
    if posA==-1 or posB==-1:    return ''
    elif exclusive:             return seq[posA+lenA:posB]
    else:                       return seq[posA:posB+lenB]

def generate_seq_slices(gene_seq,slice_len,slice_step):
    """ Yield (start_pos,slice_seq) pairs with slices covering gene_seq with the given length and step size. """
    # if the gene will fit in one slice, just return that
    if slice_len>=len(gene_seq): 
        yield 0, gene_seq
        return
    for pos in range(0,len(gene_seq)-slice_len,slice_step):
        yield pos, gene_seq[pos:pos+slice_len]
    # if the sequence isn't evenly divided into correct-length slices, make one last slice covering the last slice_len 
    #  of the gene (the difference between this slice and the previous one will be under slice_step)
    if pos<len(gene_seq):
        yield len(gene_seq)-slice_len, gene_seq[-slice_len:]

