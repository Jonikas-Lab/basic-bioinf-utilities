#! /usr/bin/env python
"""
Given any number of fasta or fastq files, print a list of read lengths and counts of reads with that length.  The file type is detected automatically, using the HTSeq library FastaReader/FastqReader methods.
 --Weronika Patena, June 2011
"""

from HTSeq import FastaReader, FastqReader

import sys
from collections import defaultdict

def get_lengths(file_iterator,readlen_counter):
    """ given an iterator over read sequences, add read counts for each read length to the readlen_counter dictionary. """
    total_count = 0
    for read in file_iterator: 
        readlen_counter[len(read)] += 1
        total_count += 1
    return total_count

def print_lengths(readlen_counter,include_zeros=False,verbose=1):
    """ Given a length:readcount dictionary, print it prettily to stdout. """
    lengths = readlen_counter.keys()
    if include_zeros:   lengths = range(min(lengths),max(lengths)+1)
    else:               lengths.sort()
    if verbose>0:     print "length\tread count"
    for l in lengths:
        print "%s\t%s"%(l,readlen_counter[l])
    if verbose>0:     print "Total %s reads"%sum(readlen_counter.values())

def read_length_distribution(infiles,include_zeros=False,verbose=1):
    """ Main program: Given a list of fastq/fasta files, print a list of read lengths and counts."""
    # a counter with a default value of 0
    readlen_counter = defaultdict(lambda: 0)
    for infile in infiles:
        # FILETYPE RECOGNITION: first try Fasta: if it works, good. If it fails with some random error, raise it. 
        #   If it fails with a "file isn't fasta" exception, try fastq.  If fastq failse with some random error, raise it.
        #   If fastq fauls with a "file isn't fastq" exception, print the info and BOTH exceptions and exit.
        # Note that the HTSeq FastaReader/FastqReader fails when you try to do something, not when you initiate it!
        #   (so I have to actually run get_lengths inside this try/except, instead of just defining the reader)
        try:    
            total_count = get_lengths(FastaReader(infile),readlen_counter)
            if verbose>1: print "File %s recognized as fasta, with %s total reads."%(infile,total_count)
        except (ValueError,AssertionError) as fasta_error_msg:
            if not str(fasta_error_msg).count("FASTA file does not start with"):  raise 
            try:    
                total_count = get_lengths(FastqReader(infile,qual_scale="solexa"),readlen_counter)
                if verbose>1: print "File %s recognized as fastq, with %s total reads."%(infile,total_count)
            except (ValueError,AssertionError) as fastq_error_msg:
                if not str(fastq_error_msg).count("this is not FASTQ data"):   raise
                sys.exit("Error: input file %s is not recognized as either fasta or fastq. \
                         \n\tFasta error message: %s \n\tFastq error message: %s"\
                         %(infile,fasta_error_msg,fastq_error_msg))
    print_lengths(readlen_counter,include_zeros)
            

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(__doc__)
    parser.add_option('-z','--include_zeros', action='store_true', default=False, help="Include lengths with 0 counts in the output (output the whole range from lowest to highest length, instead of just the lengths present; default %default).")
    parser.add_option('-v','--verbose', action='store_const', const=2, dest="verbose", default=1, help="Print file type and read count info about each processed file (default off).")
    parser.add_option('-q','--quiet', action='store_const', const=0, dest="verbose", help="Don't print the header or the total number of counts (default off).")  # dest=verbose means -v and -q modify the same variable, with default 1
    (options, args) = parser.parse_args()
    if not args:
        parser.print_help()
        sys.exit("\nError: At least one argument file required.")
    infiles = args

    read_length_distribution(infiles,options.include_zeros,options.verbose)
