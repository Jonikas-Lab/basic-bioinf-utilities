#!/usr/bin/env python
""" 
Take a FASTQ file, remove a specified number of bases from the 5' and/or 3' end, output as fastq or fasta.
Author: Weronika Patena, 2010-2011 (using some code by Hunter Richards)

USAGE:  strip_seq_from_fastq.py [options] infile
"""

import sys

class FastQ_Generator:
    """ An iterator that, given a filename, will give successive FastQ records from the file (by Hunter)."""
    def __init__(self, filename):
        self.file = open(filename, 'rU')
    def __iter__(self):
        return self
    def next(self):
        # if this fails (end of file), that just means the iteration's over, it'll raise StopIteration for the caller.
        header1 = self.file.next().strip()
        # if any of the next three lines aren't in the infile or don't match the header right, raise an error
        try:
            sequence = self.file.next().strip()
            header2 = self.file.next().strip()
            quality = self.file.next().strip()
            if not header1[0]=='@' or not header2[0]=='+':  raise ValueError
            header1,header2 = header1[1:],header2[1:]
            if not (len(header2)==0 or header1==header2):   raise ValueError
            if not len(sequence)==len(quality):             raise ValueError
        except (StopIteration,ValueError):
            raise ValueError('Input FastQ file is malformed! At record %s'%header1)
        return (header1, sequence, quality)


if __name__ == "__main__":

    ### option parser
    from optparse import OptionParser
    parser = OptionParser(__doc__,add_help_option=False)
    parser.add_option('-5','--remove_5prime_bases',type='int', default=0, metavar='N', 
                      help="Trim the first N bases from the 5' end (default %default)")
    parser.add_option('-3','--remove_3prime_bases',type='int', default=0, metavar='N', 
                      help="Trim the first N bases from the 3' end (default %default)")
    parser.add_option('-f','--output_fasta', action='store_true', default=False, 
                      help="Output fasta instead of the default fastq")
    parser.add_option('-k','--keep_quality_header', action='store_true', default=False, 
                      help="Print the full header instead of just a + before the quality line (default %default)")
    global options
    (options, args) = parser.parse_args()
    try:
        [infile] = args
    except ValueError:
        parser.print_help()
        sys.exit("\nError: exactly one input file required")

    # convert options.remove_3prime_bases to useful list-slice format:
    #   if we're removing last 3 bases, we want to do seq[:-3]; if we're removing none, seq[:None]
    if options.remove_3prime_bases: options.remove_3prime_bases *= -1
    else:                           options.remove_3prime_bases = None
    # options.remove_5prime_bases can stay as is, since seq[1:] removes the first element of seq, etc

    # go over all sequences, remove bases as required, print all the sequences
    fastq_iterator = FastQ_Generator(infile)
    for (ID,sequence,quality) in fastq_iterator:
        trimmed_sequence = sequence[options.remove_5prime_bases:options.remove_3prime_bases]
        trimmed_quality = quality[options.remove_5prime_bases:options.remove_3prime_bases]
        if options.output_fasta:
            print ">%s"%ID
            print trimmed_sequence
        else:          
            print "@%s"%ID
            print trimmed_sequence
            if options.keep_quality_header: print "+%s"%ID
            else:                           print "+"
            print trimmed_quality
