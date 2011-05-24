#! /usr/bin/env python

"""
Utilities to parse output files from different RNA folding programs, etc.
"""

import sys

error_message = "Problem parsing file %s - unexpected position %s instead of %s in column 1."


def parse_unafold_file(foldfile,n_bases=1):
    """ parse a unafold *.ext output file, return a list of single-stranded probabilities at each position. """
    ss_prob_by_pos = []
    for line in open(foldfile):
        line = line.strip()
        # ignore header line
        if line[0]=='i':    continue
        line_fields = line.split('\t')
        try:    
            pos = int(line_fields[0])
            if n_bases==1:      ss_prob = float(line_fields[1])
            elif n_bases==2:    ss_prob = float(line_fields[2])
        except (ValueError, IndexError):  
            print "Warning: couldn't parse line \"%s\" in file %s, skipping line."%(line,foldfile)
            continue
        ss_prob_by_pos.append(ss_prob)
        # I'm just appending the lines from the file in order instead of by position - check that the position matches
        if not pos==len(ss_prob_by_pos):    sys.exit(error_message%(foldfile, pos,len(ss_prob_by_pos)))
    return ss_prob_by_pos


def parse_sfold_file(foldfile,n_bases=1):
    """ parse an sfold sstrand.out output file, return a list of single-stranded probabilities at each position. """
    ss_prob_by_pos = []
    import re
    for line in open(foldfile):
        line = line.strip()
        # remove all extra spaces first!  The numbers of spaces separating fields may not be consistent.
        line = re.sub('\s\s+',' ',line)
        line_fields = line.split(' ')
        try:    
            pos = int(line_fields[0])
            if n_bases==1:      ss_prob = float(line_fields[3])
            elif n_bases==2:    ss_prob = float(line_fields[4])
            elif n_bases==4:    ss_prob = float(line_fields[5])
        except (ValueError, IndexError):  
            print "Warning: couldn't parse line \"%s\" in file %s, skipping line."%(line,foldfile)
            continue
        ss_prob_by_pos.append(ss_prob)
        # I'm just appending the lines from the file in order instead of by position - check that the position matches
        if not pos==len(ss_prob_by_pos):    sys.exit(error_message%(foldfile, pos,len(ss_prob_by_pos)))
    return ss_prob_by_pos


def parse_mfold_file(foldfile,n_bases=1):
    """ parse an mfold ss-count output file, return a list of single-stranded probabilities at each position. """
    ss_prob_by_pos = []
    total_folds = 0
    for line in open(foldfile):
        line = line.strip()
        line_fields = line.split(' ')
        try:    
            # read the number of total folds from the header line if it hasn't been read yet
            if not total_folds and len(line_fields)==1:     
                total_folds = int(line)
                continue
            pos = int(line_fields[0])
            ss_prob = float(line_fields[1])/total_folds
        except (ValueError, IndexError):  
            print "Warning: couldn't parse line \"%s\" in file %s, skipping line."%(line,foldfile)
            continue
        ss_prob_by_pos.append(ss_prob)
        # I'm just appending the lines from the file in order instead of by position - check that the position matches
        if not pos==len(ss_prob_by_pos):    sys.exit(error_message%(foldfile, pos,len(ss_prob_by_pos)))
    return ss_prob_by_pos
