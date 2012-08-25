#! /usr/bin/env python
"""
Given any number of fasta or fastq files, print the N most common sequences in it (either full sequences, or Xbp prefixes) along with their counts as raw numbers and as percentages of the total count.  Optionally only print them if the percentage of the total count exceeds a minimum.
NOTE: fastq file qualities will always be written using Sanger (+33) encoding! (HTSeq only implements that)
 --Weronika Patena, Feb 2012
"""

# basic libraries
import sys, os
from collections import defaultdict
# other libraries
from Bio import SeqIO
# my modules
from deepseq_utilities import get_seq_count_from_collapsed_header, check_fasta_fastq_format


def subsequence_counts(infile_reader, seq_length=None, input_collapsed_to_unique=False):
    """ Given an iterator of Biopython seq objects and desired subsequence length/end, return subsequence:count dict.

    seq_length: if None, take whole seq; if N>0, take first N bases, if N<0, take last -N bases.
    If input_collapsed_to_unique is True, consider each sequence to be X reads, determined from seq.name, 
     using the format used by fastx_collapser from FastX Toolkit.
    """
    seq_counter = defaultdict(lambda: 0)    # a counter with a default value of 0
    for sequence in infile_reader: 
        N_seqs = get_seq_count_from_collapsed_header(sequence.name) if input_collapsed_to_unique else 1
        # using seq.data to convert Biopython Seq objects to plain strings - Seq objects aren't hashable correctly
        if seq_length > 0:  subsequence = sequence.seq[0:seq_length].data
        else:               subsequence = sequence.seq[seq_length:].data
        seq_counter[subsequence] += N_seqs
    return dict(seq_counter)


def seq_top_sequence_check(infiles, seq_length=None, n_to_print=3, min_percent_to_print=None, input_collapsed_to_unique=False):
    """ See module docstring and optparse option help messages, to avoid repeating the info. """

    if seq_length is None:  seqlen_info = ''
    elif seq_length>0:      seqlen_info = ' first %sbp'%seq_length
    elif seq_length<0:      seqlen_info = ' last %sbp'%(-seq_length)

    for infile in infiles:
        seq_format = check_fasta_fastq_format(infile)

        with open(infile) as INFILE:
            seq_counter = subsequence_counts(SeqIO.parse(INFILE, seq_format), seq_length, input_collapsed_to_unique)

        seq_list_by_count = sorted(seq_counter.items(), key=lambda (s,c): c, reverse=True)

        total_seqs = sum(seq_counter.values())

        # if not using the min_percent_to_print option, just print the top N sequences from each file
        if min_percent_to_print is None:
            seq_data_list = []
            for i in range(min(n_to_print,len(seq_list_by_count))):
                seq,count = seq_list_by_count[i]
                percent = count*100.0/total_seqs
                # "%.2g" is significant-digit-based formatting of floats!!  So 92.12345 is 92%, but 0.00045 is 0.00045%.
                percent_2_sig_digits = str(float("%.2g"%percent))
                if percent_2_sig_digits.endswith(".0"):     percent_2_sig_digits = percent_2_sig_digits[:-2]
                seq_data_list.append("%s%% %s (%d)"%(percent_2_sig_digits, seq, count))
            print " * %s (%s seqs, %s unique%s):"%(infile, total_seqs, len(seq_list_by_count), seqlen_info) 
            print ', '.join(seq_data_list)

        # if using the min_percent_to_print option, just print the top N sequences from each file
        else:
            print "min_percent_to_print NOT IMPLEMENTED!"
            # TODO implement!!


    # MAYBE-TODO could also make it output stuff in a table for all the files?  Especially short prefixes...

# MAYBE-TODO should this have unit-tests?


if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(__doc__)
    parser.add_option('-l','--seq_length', type='int', default=None, metavar='X',
                      help="Print most common Xbp prefixes (if X is positive) or suffixes (if X is negative) "
                          +"(None - full sequences) (default %default).")
    parser.add_option('-n','--n_to_print', type='int', default=3, metavar='N',
                      help="Print the N most common sequences (default %default).")
    parser.add_option('-m','--min_percent_to_print', type='int', default=None, metavar='M',
                      help="Only print the sequences that constitute over M% of total reads,"
                          +"INSTEAD OF using -n (default %default).")
    parser.add_option('-c','--input_collapsed_to_unique', action='store_true', default=False, 
                      help="Use this to get correct total counts if the infile was collapsed to unique sequences using fastx_collapser, with original sequence counts encoded in the headers (a '>2-572' header means there were 572 identical sequences); default %default).")
    (options, infiles) = parser.parse_args()

    seq_top_sequence_check(infiles, options.seq_length, options.n_to_print, options.min_percent_to_print, options.input_collapsed_to_unique)

    # MAYBE-TODO add an option to do multiple lengths per run, so that they're together for each file, instead of separate?  Like this:
    #  * cutadapt_RZ-2_prev_trimmed/10bp.fa (7148 seqs, 77 unique, 21 unique first 4bp, 69 unique first 8bp):
    #   88% ACTACCTACT (6262), 2.2% ACTAGTCCCA (159), 0.95% ACTACGATCT (68)
    #   98% ACTA (6981), 0.7% GGGT (50), 0.34% ACCT (24), 0.24% ACTG (17), 0.22% AACA (16)
    #   88% ACTACCTA (6263), 2.6% ACTAGTCC (186), 0.97% ACTACGAT (69), 0.8% ACTAGTAT (57), 0.7% GGGTCCAA (50)

