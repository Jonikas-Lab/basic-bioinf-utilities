#! /usr/bin/env python2.7

"""
Convert the odd text output format from cutadapt -r option to fastq
 -- Weronika Patena, 2014
USAGE: program infile outfile
"""

# standard library
from __future__ import division
import sys
# other packages
# my modules

if __name__=='__main__':
    try:
        infile, outfile = sys.argv[1:]
    except ValueError:
        raise Exception("Exactly one infile and one outfile required!")

    with open(outfile, 'w') as OUTFILE:
        for line in open(infile):
            try:
                seq, ID = line.strip().split(' ',1)
            except ValueError:
                raise Exception("Can't parse line! %s"%line)
            OUTFILE.write('>%s\n%s\n'%(ID,seq))
