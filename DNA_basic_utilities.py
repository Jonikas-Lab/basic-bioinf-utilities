#! /usr/bin/env python

"""
Various basic DNA utilities I wrote (parse_fasta, complement, reverse_complement, find_seq_between, generate_seq_slices, ...) - most of them just imported from their own files, since I need them as separate programs too.
Weronika Patena, 2008-2010
"""

### basic library
import unittest

### my modules
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

### some constants 

SEQ_ENDS = ['5prime','3prime']
SEQ_STRANDS = ['+','-']
SEQ_DIRECTIONS = ['forward','reverse']
SEQ_ORIENTATIONS = ['sense','antisense']


### testing whether two sequences contain one another, or overlap (based purely on position, no sequence involved)

def position_test_contains(seq1_start, seq1_end, seq2_start, seq2_end):
    """ Return True if seq2 is completely contained within seq1, False otherwise. All positions should be numbers."""
    # make sure the arguments make sense
    if seq1_start > seq1_end or seq2_start > seq2_end:
        raise ValueError("sequence start positions must not be higher than end positions!")
    # the test is simple
    return bool(seq1_start <= seq2_start and seq2_end <= seq1_end)

def position_test_overlap(seq1_start, seq1_end, seq2_start, seq2_end, return_overlap=False):
    """ Returns True if seq1/seq2 overlap, else False; or the number of bases of overlap if return_overlap is True (min 0).
    Assumes the positions are inclusive (i.e. length = end-start+1). """
    # make sure the arguments make sense
    if seq1_start > seq1_end or seq2_start > seq2_end:
        raise ValueError("sequence start positions must not be higher than end positions!")
    ### if we don't want to measure the overlap, the test is simple
    if return_overlap is False:     return bool(seq2_start <= seq1_end and seq2_end >= seq1_start)
    ### if we do want to measure the overlap, need to do more work:
    # if seq1 contains seq2, or the opposite, return the length of the shorter sequence (see docstring for why +1)
    if position_test_contains(seq1_start, seq1_end, seq2_start, seq2_end):  return seq2_end - seq2_start + 1
    if position_test_contains(seq2_start, seq2_end, seq1_start, seq1_end):  return seq1_end - seq1_start + 1
    # if seq2 starts before seq1, or ends after seq1, the overlap is the distance between the "inner" start and end
    #   (there's a +1 in both cases, because if seq1_end==seq2_start, there's 1bp overlap)
    if seq2_start < seq1_start:                                             return max(seq2_end - seq1_start + 1, 0)  
    if seq2_end > seq1_end:                                                 return max(seq1_end - seq2_start + 1, 0)  


### other position-based functions

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



class Testing_everything(unittest.TestCase):
    """ Testing all functions/classes.etc. """

    def test__position_test_contains(self):
        # raises an error if start>end for either sequence
        self.assertRaises(ValueError, position_test_contains, 11,10,1,1)
        self.assertRaises(ValueError, position_test_contains, 1,1,11,10)
        # seq1 contains seq2 - True
        assert position_test_contains(1,10,2,3)==True
        assert position_test_contains(1,10,9,10)==True
        assert position_test_contains(1,10,10,10)==True
        assert position_test_contains(10,10,10,10)==True
        # seq2 contains seq1 - False
        assert position_test_contains(2,3,1,10)==False
        assert position_test_contains(9,10,1,10)==False
        assert position_test_contains(10,10,1,10)==False
        # seq2 and seq1 overlap - False
        assert position_test_contains(2,3,3,10)==False
        assert position_test_contains(2,4,3,10)==False
        assert position_test_contains(3,10,2,3)==False
        assert position_test_contains(3,10,2,4)==False
        # seq2 and seq1 don't overlap - False
        assert position_test_contains(2,3,4,10)==False
        assert position_test_contains(2,3,6,10)==False
        assert position_test_contains(4,10,2,3)==False
        assert position_test_contains(6,10,2,3)==False


    def test__position_test_overlap__returns_bool(self):
        # raises an error if start>end for either sequence
        self.assertRaises(ValueError, position_test_overlap, 11,10,1,1, False)
        self.assertRaises(ValueError, position_test_overlap, 1,1,11,10, False)
        # seq1 contains seq2 - True
        assert position_test_overlap(1,10,2,3,False)==True
        assert position_test_overlap(1,10,9,10,False)==True
        assert position_test_overlap(1,10,10,10,False)==True
        assert position_test_overlap(10,10,10,10,False)==True
        # seq2 contains seq1 - True
        assert position_test_overlap(2,3,1,10,False)==True
        assert position_test_overlap(9,10,1,10,False)==True
        assert position_test_overlap(10,10,1,10,False)==True
        # seq2 and seq1 overlap - True
        assert position_test_overlap(2,3,3,10,False)==True
        assert position_test_overlap(2,4,3,10,False)==True
        assert position_test_overlap(3,10,2,3,False)==True
        assert position_test_overlap(3,10,2,4,False)==True
        # seq2 and seq1 don't overlap - False
        assert position_test_overlap(2,3,4,10,False)==False
        assert position_test_overlap(2,3,6,10,False)==False
        assert position_test_overlap(4,10,2,3,False)==False
        assert position_test_overlap(6,10,2,3,False)==False

    def test__position_test_overlap__returns_number(self):
        # raises an error if start>end for either sequence
        self.assertRaises(ValueError, position_test_overlap, 11,10,1,1, True)
        self.assertRaises(ValueError, position_test_overlap, 1,1,11,10, True)
        # seq1 contains seq2 - True
        assert position_test_overlap(1,10,2,3,True)==2
        assert position_test_overlap(1,10,9,10,True)==2
        assert position_test_overlap(1,10,10,10,True)==1
        assert position_test_overlap(10,10,10,10,True)==1
        # seq2 contains seq1 - True
        assert position_test_overlap(2,3,1,10,True)==2
        assert position_test_overlap(9,10,1,10,True)==2
        assert position_test_overlap(10,10,1,10,True)==1
        # seq2 and seq1 overlap - True
        assert position_test_overlap(2,3,3,10,True)==1
        assert position_test_overlap(2,4,3,10,True)==2
        assert position_test_overlap(3,10,2,3,True)==1
        assert position_test_overlap(3,10,2,4,True)==2
        # seq2 and seq1 don't overlap - False
        assert position_test_overlap(2,3,4,10,True)==0
        assert position_test_overlap(2,3,6,10,True)==0
        assert position_test_overlap(4,10,2,3,True)==0
        assert position_test_overlap(6,10,2,3,True)==0


    def test__find_seq_between(self):
        pass
    # TODO implement!


    def test__generate_seq_slices(self):
        pass
    # TODO implement!


if __name__=='__main__':
    """ If module is ran directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()

