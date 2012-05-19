#! /usr/bin/env python
"""
Given any number of fasta or fastq files, print a list of seq lengths and counts of seqs with that length (or just the total seq count).  The file type is detected automatically, using the HTSeq library FastaReader/FastqReader methods.
 --Weronika Patena, June 2011
"""

# basic libraries
import sys
from collections import defaultdict
# other libraries
from HTSeq import FastaReader, FastqReader
# my modules
from deepseq_utilities import get_seq_count_from_collapsed_header


def _get_count_and_lengths(file_iterator, seqlen_counter=None, input_collapsed_to_unique=False):
    """ Given an iterator over sequences, add seq counts for each seq length to the seqlen_counter dictionary. """
    total_count = 0
    for seq in file_iterator: 
        N_seqs = get_seq_count_from_collapsed_header(seq.name) if input_collapsed_to_unique else 1
        if not seqlen_counter is None:
            seqlen_counter[len(seq)] += N_seqs
        total_count += N_seqs
        # MAYBE-TODO could reimplement it using file_iterator.get_sequence_lengths() function! Gives a header:len dict... But I'm not sure I WANT to keep that whole dict in memory even for a short time, if it's a large file, so maybe best not.
    return total_count


def _format_lengths(seqlen_counter, include_zeros=False, verbosity=1):
    """ Given a length:seqcount dictionary, format it for printing, return list of lines. """
    output_lines = []
    lengths = seqlen_counter.keys()
    if include_zeros:   lengths = range(min(lengths),max(lengths)+1)
    else:               lengths.sort()
    if verbosity>0:     
        output_lines.append("length\tseq count\n")
    for l in lengths:
        output_lines.append("%s\t%s\n"%(l,seqlen_counter[l]))
    if verbosity>0:     
        output_lines.append("Total %s seqs\n"%sum(seqlen_counter.values()))
    return output_lines


def seq_count_and_lengths(infiles, total_seq_number_only=False, input_collapsed_to_unique=False, 
                          include_zeros=False, verbosity=1, OUTPUT=sys.stdout):
    """ Given a list of fastq/fasta files, return and print a list of seq lengths and counts and the total seq count.
    
    If total_seq_number_only is True, only return/print total seq count.
    If input_collapsed_to_unique is True, program assumes infile was preprocessed with fastx_collapser, 
     and attempts to give original pre-collapsing seq_counts (based on headers).
    If include_zeros is False (default), only print non-zero seq counts; if True, print seq counts for 
     all lengths between min and max length, even if they're 0.
    Verbosity: if >1, print filetype and seqcount for each input file; if 0, don't print header or summary.
    Prints to stdout by default; to print to file, pass open file object as OUTPUT; to suppress printing, pass None."""

    # a counter with a default value of 0
    if total_seq_number_only:   seqlen_counter = None
    else:                       seqlen_counter = defaultdict(lambda: 0)
    total_seqcount = 0
    formatted_output = []
    # add the numbers from each file to seqlen_counter (using HTSeq, which auto-detects filetype etc)
    for infile in infiles:
        # FILETYPE RECOGNITION: first try Fasta: if it works, good. If it fails with some random error, raise it. 
        #   If it fails with a "file isn't fasta" exception, try fastq.  If fastq failse with some random error, raise it.
        #   If fastq fauls with a "file isn't fastq" exception, print the info and BOTH exceptions and exit.
        # Note that the HTSeq FastaReader/FastqReader fails when you try to do something, not when you initiate it!
        #   (so I have to actually run _get_count_and_lengths inside this try/except, instead of just defining the reader)
        # TODO actually this doesn't count empty sequences correctly! It only counts 1 0-length sequence even when there are many (see test_inputs/test.fa).  Fix that?
        try:    
            file_seqcount = _get_count_and_lengths(FastaReader(infile), seqlen_counter, input_collapsed_to_unique)
            if verbosity>1: 
                formatted_output.append("File %s recognized as fasta, with %s total seqs.\n"%(infile,file_seqcount))
        except (ValueError,AssertionError) as fasta_error_msg:
            if not str(fasta_error_msg).count("FASTA file does not start with"):  raise 
            try:    
                file_seqcount = _get_count_and_lengths(FastqReader(infile,qual_scale="solexa"), 
                                                       seqlen_counter, input_collapsed_to_unique)
                if verbosity>1: 
                    formatted_output.append("File %s recognized as fastq, with %s total seqs.\n"%(infile,file_seqcount))
            except (ValueError,AssertionError) as fastq_error_msg:
                if not str(fastq_error_msg).count("this is not FASTQ data"):   raise
                sys.exit("Error: input file %s is not recognized as either fasta or fastq. \
                         \n\tFasta error message: %s \n\tFastq error message: %s"\
                         %(infile,fasta_error_msg,fastq_error_msg))
        total_seqcount += file_seqcount
    # format and print (optionally) and return the output
    if total_seq_number_only:
        formatted_output.append("Total %s seqs\n"%total_seqcount)
    else:
        formatted_output += _format_lengths(seqlen_counter, include_zeros, verbosity)
    if not OUTPUT is None:
        for line in formatted_output:   OUTPUT.write(line)
    return formatted_output
            

# MAYBE-TODO should this have unit-tests?


if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(__doc__)
    parser.add_option('-n','--total_seq_number_only', action='store_true', default=False, 
                      help="Output only the total seq number, without the seq length distribution; default %default).")
    parser.add_option('-c','--input_collapsed_to_unique', action='store_true', default=False, 
                      help="Use to get correct total counts if the infile was collapsed to unique sequences using fastx_collapser, with original sequence counts encoded in the headers (a '>2-572' header means there were 572 identical sequences); default %default).")
    parser.add_option('-z','--include_zeros', action='store_true', default=False, 
                      help="Include lengths with 0 counts in the output (output the whole range from lowest to highest length, instead of just the lengths present; default %default).")
    # MAYBE-TODO add -u option to keep track of only unique sequences?
    # -v and -q modify the same variable (verbosity) - default 1, -v makes it 2, -q makes it 0.
    parser.add_option('-v','--verbose', action='store_const', const=2, dest="verbosity", default=1, 
                      help="Print file type and seq count info about each processed file (default off).")
    parser.add_option('-q','--quiet', action='store_const', const=0, dest="verbosity", 
                      help="Don't print the header or the total number of counts (default off).")
    (options, args) = parser.parse_args()
    if not args:
        parser.print_help()
        sys.exit("\nError: At least one argument file required.")
    infiles = args

    seq_count_and_lengths(infiles, options.total_seq_number_only, options.input_collapsed_to_unique, options.include_zeros, options.verbosity)
