#! /usr/bin/env python2.7

"""
______
 -- Weronika Patena, 2019
USAGE: _____
"""

# standard library
from __future__ import division
import sys
import collections
import math
# other packages
import matplotlib.pyplot as mplt
# my modules


def parse_blast(infile):
    """ Parses a tabularized blast result file (generated with -m9 option). 
    
    Returns query:list_of_results dict, where each item in the list of results is a namedtuple containing all the blast data:  
        subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score.
    """
    # I can skip the header lines, because each result line gives me the query + subject. 
    H = 'subject_id identity_percent aln_length mismatches gaps query_start query_end subject_start subject_end e_value bit_score'
    Blast_aln = collections.namedtuple('Blast_aln', H.split())
    result_dict = collections.defaultdict(list)
    for line in open(infile):
        if line.startswith('# '):   continue
        f = line.split()
        fields = f[:2] + [float(f[2])] + [int(x) for x in f[3:10]] + [float(x) for x in f[10:]]
        result_dict[fields[0]].append(Blast_aln(*fields[1:]))
    return result_dict


def plot_blast(infile=None, input_data=None, 
               max_evalue=None, max_evalue_by_subject=None, subject_order=None, skip_cassette_overlaps=False, 
               subject_colors=None, cassette_gradient_cmap=None, query_lengths=None, 
               ncols=1, figsize=(12,8), titles=True, labels=True):
    """ Plot all the blast (or similar) alignments for each read.

    Infile should be a tabularized blast result file (generated with -m9 option).
    Alternatively, provide input_data, which should be a name:list_of_alignments result, where each alignment is a namedtuple
     or object with the following fields: subject_id query_start query_end subject_start subject_end e_value aln_length.
     (subject = chromosome the alignment is to; query = read).
    """
    if max_evalue is not None and max_evalue_by_subject is not None:
        raise Exception("Only specify max_evalue or max_evalue_by_subject, not both!")
    if (infile is None) + (input_data is None) != 1:
        raise Exception("Need one infile or input_data, but not both!")
    if infile is not None:
        input_data = parse_blast(infile)
    if cassette_gradient_cmap is not None:
        cassette_length = 2223          # TODO at some point I should probably provide an infile with all the subject lengths/etc?
    N = len(input_data)
    nrows = math.ceil(N/ncols)
    mplt.figure(figsize=figsize)
    for n, (query, alns) in enumerate(input_data.items()):
        if skip_cassette_overlaps:
            cassette_positions = [(aln.query_start, aln.query_end) for aln in alns if 'insertion_cassette' in aln.subject_id]
            filtered_alns = []
            for aln in alns:
                if not any(x in aln.subject_id for x in 'cassette LEAPseq adapter MinION'.split()):
                    if any(((start - aln.query_start) < 10 and (aln.query_end - end) < 10) for (start,end) in cassette_positions):
                        continue
                filtered_alns.append(aln)
        else:
            filtered_alns = alns
        if max_evalue is not None:
            filtered_alns = [aln for aln in filtered_alns if aln.e_value <= max_evalue]
        if max_evalue_by_subject is not None:
            filtered_alns = [aln for aln in filtered_alns if aln.e_value <= max_evalue_by_subject[aln.subject_id]]
        all_subjects = set(aln.subject_id for aln in filtered_alns)
        if subject_order is not None:
            subjects_sorted = sorted(all_subjects, key = lambda s: (subject_order(s), s))
        else:
            subjects_sorted = sorted(all_subjects)
        mplt.subplot(nrows, ncols, n+1)
        if query_lengths is not None:   query_len = query_lengths[query]
        else:                           query_len = max(aln.query_end for aln in alns) * 1.1
        mplt.barh(0.7, query_len, 0.3, 0, color='black')
        for i, subject in enumerate(subjects_sorted):
            curr_alns = [aln for aln in filtered_alns if aln.subject_id == subject]
            if ('insertion_cassette' in subject) and (cassette_gradient_cmap is not None):
                for aln in curr_alns:
                    # colormaps are 0-255, so scale the values appropriately, and use vmin/vmax 0/255 to make sure they stay there
                    color1 = aln.subject_start/cassette_length*256
                    color2 = aln.subject_end/cassette_length*256
                    mplt.imshow([[color1, color2], [color1, color2]], interpolation='bicubic', cmap=cassette_gradient_cmap, 
                              extent=(aln.query_start, aln.query_end, -i-0.5, -i+0.3), alpha=1, aspect='auto', vmin=0, vmax=255)
            else:
                kwargs = dict(edgecolor='None')
                if subject_colors is not None:  kwargs['color'] = subject_colors[subject]
                mplt.barh([-i-0.5 for aln in curr_alns], 
                          [aln.aln_length for aln in curr_alns], 0.8, 
                          [aln.query_start for aln in curr_alns], **kwargs)
            # put end markers on the cassette (triangle down for 5', up for 3');  
            #   markers on both ends of MinION/LEAPseq adapters (square for outer end, + for inner)
            #   * markers around the spots where the IB should be
            for (name, pos1, pos2, marker1, marker2, dist) in [('insertion_cassette_CIB1',      0, 2222, 'v', '^',    30),
                                                               ('MinION_adapter',               0, 60,   'o', '_',    10),
                                                               ('LEAPseq_genome_side',          0, 42,   '_', 'o',    10),
                                                               ('LEAPseq_cassette_side_5prime', 0, 189,  'o', '_',    10),
                                                               ('LEAPseq_cassette_side_3prime', 0, 156,  'o', '_',    10),
                                                               ('insertion_cassette_CIB1',      44, 67,     '*', '*', 5),
                                                               ('insertion_cassette_CIB1',      2155, 2178, '*', '*', 5),
                                                               ('LEAPseq_cassette_side_5prime', 122, 145,   '*', '*', 5),
                                                               ('LEAPseq_cassette_side_3prime', 89, 112,    '*', '*', 5)]:
                if subject == name:
                  for aln in curr_alns:
                    for subject_pos, query_pos in [(aln.subject_start, aln.query_start), (aln.subject_end, aln.query_end)]:
                      if abs(subject_pos-pos1) <= dist:    mplt.plot(query_pos, -i-0.1, marker=marker1, color='k', markersize=5)
                      if abs(subject_pos-pos2) <= dist:    mplt.plot(query_pos, -i-0.1, marker=marker2, color='k', markersize=5)
        # plot the query again to keep the xrange from getting screwed up by imshow
        mplt.barh(0.7, query_len, 0.3, 0, color='black')
        mplt.xticks(range(0, int(query_len), 100), [])
        mplt.xlim(query_len * -0.02, query_len * 1.02)
        if titles:  mplt.title(query)
        if labels:  mplt.yticks(range(-len(all_subjects), 1), list(reversed(subjects_sorted)) + ['query'])
        else:       mplt.yticks(range(-len(all_subjects), 1), [])   
        if not titles and not labels:   mplt.ylabel(n+1, rotation=0)
