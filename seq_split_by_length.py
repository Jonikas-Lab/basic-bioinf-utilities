#! /usr/bin/env python
"""
Given one fasta or fastq file, split the sequences by length, and output a folder matching the input filename (no extension), containing files named <original_filename>_20bp.fa etc.
The file type is detected automatically, using the HTSeq library FastaReader/FastqReader methods.
NOTE: fastq file qualities will always be written using Sanger (+33) encoding! (HTSeq only implements that)
 --Weronika Patena, Feb 2012
"""

# basic libraries
import sys, os
from collections import defaultdict
# other libraries
from HTSeq import FastaReader, FastqReader
# my modules
from seq_count_and_lengths import _format_lengths


def seq_split_by_length(infile, min_length=None, max_length=None, include_empty_files=False, 
                        input_collapsed_to_unique=False, verbosity=1):
    """ Given one fasta or fastq file, split the sequences by length. 
    Output a folder matching the input filename (no extension), containing files named <original_filename>_20bp.fa etc.
    The file type is detected automatically, using the HTSeq library FastaReader/FastqReader methods.
    Arguments:
     * min_length - output all sequences shorter than X to a shorter.fa file.
     * max_length - output all sequences longer than Y to a longer.fa file.
     * include_empty_files - also output (empty) files for lengths with no sequences with lengths between min and max in order not to have files missing from the range. If min_length/max_length are set, these will be taken as min/max; otherwise the actual minimum and maximum encountered sequence lengths will be used.
     * verbosity - verbosity level; print nothing to stdout if 0.
    """
    # file format recognition (I could do it by trying to use FastaReader/FastqReader on it, but it's annoying)
    fasta_extensions = ['fa','fasta']
    fastq_extensions = ['fq','fastq']
    extension = os.path.splitext(infile)[1].lower()[1:]
    if extension in fasta_extensions:   
        infile_reader = FastaReader(infile)
    elif extension in fastq_extensions: 
        infile_reader = FastqReader(infile,qual_scale="solexa")
    else:       sys.exit("Error: input file %s (extension %s) needs to have a %s extension to be recognized!"%(infile, 
                            extension, '/'.join(fasta_extensions+fastq_extensions)))

    ### make the output folder, and outfiles
    infile_base = os.path.splitext(infile)[0]
    outfolder = infile_base
    os.mkdir(outfolder)
    # a (length: open file object) dictionary, so I can keep them all open and close them at the end. 
    # Yes, I know I should really be using with/as, but I don't think you can do multiples of that at once, and I can't have a level of indent for every possible sequence length!
    len_to_outfile_dict = {}

    # a counter with a default value of 0
    seq_counter = defaultdict(lambda: 0)

    for seq in infile_reader: 
        seqlen = len(seq)
        # add the N_seqs to the seq counter
        N_seqs = get_seq_count_from_collapsed_header(seq.name) if input_collapsed_to_unique else 1
        seq_counter[seqlen] += N_seqs
        # special length cases for when min/max length is set
        if min_length is not None and seqlen<min_length:    seqlen = 'lower'
        elif max_length is not None and seqlen>max_length:  seqlen = 'higher'
        # make sure the output file for that length is open
        if seqlen not in len_to_outfile_dict.keys():
            len_to_outfile_dict[seqlen] = open(os.path.join(outfolder, "%s_%s.%s"%(infile_base,seqlen,extension)),'w')
        # write the sequence (fasta or fastq!) to the outfile!
        if extension in fasta_extensions:
            seq.write_to_fasta_file(len_to_outfile_dict[seqlen])
        else:
            seq.write_to_fastq_file(len_to_outfile_dict[seqlen])

    # optionally add the empty files that had no sequences of that length
    if include_empty_files:
        if min_length is None:  min_length = min(len_to_outfile_dict.keys())
        if max_length is None:  max_length = max(len_to_outfile_dict.keys())
        for seqlen in range(min_length+1,max_length):
            if seqlen not in len_to_outfile_dict.keys():
                len_to_outfile_dict[seqlen] = open(os.path.join(outfilder, "%s_%s.%s"%(infile_base,seqlen,extension)),'w')

    # close all the files
    for FILE in len_to_outfile_dict.values():
        FILE.close()

    # format and print the seq counts by length
    if verbosity:
        for line in _format_lengths(seqlen_dict, include_zeros, 1):     print(line)
            

# MAYBE-TODO should this have unit-tests?


if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(__doc__)
    parser.add_option('-m','--min_length', type='int', default=None, metavar='X',
                      help="Output all sequences shorter than X to a shorter.fa file (default %default).")
    parser.add_option('-M','--max_length', type='int', default=None, metavar='Y',
                      help="Output all sequences longer than Y to a longer.fa file (default %default).")
    parser.add_option('-e','--include_empty_files', action='store_true', default=False, 
                      help="Also output (empty) files for lengths with no sequences with lengths between min and max in order not to have files missing from the range. If -m/-M were set, these will be taken as min/max; otherwise the actual minimum and maximum encountered sequence lengths will be used.")
    parser.add_option('-c','--input_collapsed_to_unique', action='store_true', default=False, 
                      help="For seq totals printed to stdout only - use this to get correct total counts if the infile was collapsed to unique sequences using fastx_collapser, with original sequence counts encoded in the headers (a '>2-572' header means there were 572 identical sequences); default %default).")
    # -v and -q modify the same variable (verbosity) - default 1, -v makes it 2, -q makes it 0.
    parser.add_option('-q','--quiet', action='store_const', const=0, dest="verbosity", 
                      help="Don't print anything to stdout (default off).")
    (options, args) = parser.parse_args()
    try:
        [infile] = args
    except ValueError:
        parser.print_help()
        sys.exit("\nError: Exactly one input file required.")

    seq_split_by_length(infile, options.min_length, options.max_length, options.include_empty_files, options.input_collapsed_to_unique, options.verbosity)
