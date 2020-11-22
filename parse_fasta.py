#! /usr/bin/env python3
"""
Parse/reformat a fasta file
 --Weronika Patena, nov2008

Read standard fasta files:
The input is a a generator - call as :parse_fasta(open(filename)):
The output is a generator - creates on the fly an iterator over 
all the (header, seq) pairs of a fasta file, to use like this:
>> for (header,seq) in parse_fasta(open(infile)): <do stuff>
"""

import sys

debug = 0

def print_seq(line):
    """ printing help function - takes (header,seq) tuples from fasta, or line strings from text. """
    # if the line is a string instead of a (header,seq) tuple, just print it
    if isinstance(line,str):
        print(line)
    # otherwise you got a (header,seq) tuple
    else:
        (header,seq) = line
        print(">%s\n%s"%(header,seq))

def parse_fasta(INPUT, not_standard_nucleotides=False):
    """ Usage: for (header,seq) in parse_fasta(input): <do stuff>. Input can be a filename or generator. """
    header,seq = '',''
    # If the input is a filename, open it first
    if type(INPUT)==str: 
        isfile = True
        INPUT = open(INPUT)
    else:
        isfile = False
    for line in INPUT:
        # DON'T skip empty lines - there may be empty sequences in the file!
        #if not line.strip():    continue
        # if you find a header line, start a new entry, first yielding the previous entry
        # (unless there is no previous entry, i.e no header, since a sequence can be empty)
        if line and line[0]=='>':
            if header:   
                yield (header,seq)
            header = line[1:].strip()      # line[1:] to get rid of the initial >
            seq = ''
        # otherwise it's a seq line - add it to the current sequence (strip spaces/newlines)
        else: 
            # get rid of all whitespace characters by splitting and joining (there's probably a better way...)
            seq += ''.join(line.strip().split())
            # if there are illegal seq characters, exit with an error 
            #  (checking for illegal characters by using translate to delete all legal characters and seeing if there's anything left)
            if not not_standard_nucleotides: 
                # maketrans makes a translation table - translate arg1 characters to arg2 characters and delete arg3 ones
                #  (I don't know why it's a method instead of a function, it's stupid, since it doesn't use the parent string!)
                wrong_characters = seq.upper().translate(''.maketrans('', '', 'ACTGURYKMSWBDHVN.-*'))
                # in python 2 this looked like this instead (empty translation table and a separate list of chars to delete)
                # from string import maketrans
                # wrong_characters = seq.upper().translate(maketrans('',''), 'ACTGURYKMSWBDHVN.-*'})
                if wrong_characters:
                    raise Exception("Error: invalid sequence line! %s - wrong characters \"%s\""%(seq, wrong_characters))
                    # TODO shouldn't really hard-code the allowed bases...
                    # MAYBE-TODO include option for proteins to check those?  And maybe an option for just ACTG?
            if not header: 
                raise Exception("Error: Found a sequence line without a header first! %s"%seq)
    if isfile:    INPUT.close()
    # also return the last entry, if it's non-empty!
    if header:
        yield (header,seq)

### If called directly, test everything
if __name__ == '__main__':
    # TODO unfortunately if the file's wrong read_input will still read it whole before noticing!  But if the parse_fasta function is called directly, it can be given a filename and read it line-by-line.
    import read_input
    input = read_input.read_input()
    for (header,seq) in parse_fasta(input):     print_seq((header,seq))        # nice printing
