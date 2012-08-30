#!/usr/bin/env python
""" 
Remove duplicate information from FASTQ file headers (read names):
1) completely remove the header of the quality lines, since it's identical to the header of the sequence lines
2) optionally remove some fields from the header
    Assumes header structure "@instrumentID-runID:lane:tile:x:y#index/which_PE_read"
    - option --remove_first_header_field removes the instrumentID-runID field.
    - option --remove_last_header_field removes the "#index/which_PE_read" section.
   In all cases the program saves the first X (probably 10) encountered values for those fields and prints them.

-- Weronika Patena, 2010-2012

USAGE:  minimize_fastq_headers.py [options] infile outfile
"""
# standard libraries
import sys
# my modules
from seq_basic_utilities import parse_fastq

### option parser
from optparse import OptionParser
parser = OptionParser(__doc__)
parser.add_option('-f','--remove_first_header_field',action='store_true', default=False)
parser.add_option('-l','--remove_last_header_field',action='store_true', default=False)
# MAYBE-TODO more specific options, like -i for instrument ID, -l for lane, -b for barcode, -p for pair, etc
# MAYBE-TODO also implement for new header format? (see http://en.wikipedia.org/wiki/FASTQ_format)

(options, args) = parser.parse_args()
try:
    [infile,outfile] = args
except ValueError:
    parser.print_help()
    sys.exit("\nError: exactly one input and output file required!")

max_values_kept = 10

if options.remove_first_header_field:   first_field_values = set()
if options.remove_last_header_field:    last_field_values = set()

# go over all sequences, remove bases as required, print all the sequences
with open(outfile,'w') as OUTFILE:
    for (header,sequence,quality) in parse_fastq(infile):
        if options.remove_first_header_field:
            first_field,header = header.split(':',1)
            if len(first_field_values)<max_values_kept:
                first_field_values.add(first_field)
        if options.remove_last_header_field:
            header,last_field = header.rsplit('#',1)
            if len(last_field_values)<max_values_kept:
                last_field_values.add(last_field)
        OUTFILE.write( "@%s\n"%header)
        OUTFILE.write(sequence + '\n')
        OUTFILE.write( "+\n")
        OUTFILE.write(quality + '\n')

def print_info_line(field_name, field_values):
    global max_values_kept
    if len(field_values)==max_values_kept:    
        N_values = '%s or more'%max_values_kept
        extra_info = ' (showing %s)'%max_values_kept
    else:                                           
        N_values = '%s'%(len(field_values))
        extra_info = ''
    print "Removed the %s field from all reads. Found %s different values%s:"%(field_name,N_values,extra_info)
    print "    %s"%(', '.join(field_values))

if options.remove_first_header_field:   
    print_info_line('first',first_field_values)
if options.remove_last_header_field:    
    print_info_line('last',last_field_values)

# TODO add unit-tests?
