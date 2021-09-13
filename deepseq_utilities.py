#! /usr/bin/env python3

""" Various basic deepseq-related utilities I wrote - see separate function docstrings.  For importing only.
 --Weronika Patena, 2011-2021
"""

# basic libraries
from __future__ import division
import unittest
# other packages
import HTSeq    # no longer available for python2
# my modules
import basic_seq_utilities


## NOTE: for fasta/fastq (raw data) utilities see basic_seq_utilities.py


################## aligned data utilities ################## 

### NOTE: I'm currently using HTSeq for parsing SAM-format data, but I could try something else - pysam seems good.

# SAM and GFF are 1-based, with the end being the last base; HTSeq is 0-based, with the end being the base AFTER last.

######### NOTES ON THE SAM FORMAT
### Header:
# MAYBE-TODO do something useful with the SAM header?  Or at least copy it to outfile?
### Alignment line fields:
# * query template name
# * bitwise flag (relevant bits: 4 = unmapped, 16 = reverse-complement, 512 = failed quality control)
# * reference sequence name (* = unmapped)
# * leftmost mapping position (1-based) (0 = unmapped)
# * mapping quality ("-10 * log10(probability that position is wrong)"; 255 = unknown)
# * CIGAR string - descriptions of alignment matches/mismatches/etc (M/= match, I/D ins/del, X mismatch, S/H clipping)
# * (PE only - reference name of the mate fragment)
# * (PE only - position of the mate fragment)
# * template length
# * fragment sequence
# * ASCII of Phred-scaled base quality + 33   (original deepseq read quality)
# * OPTIONAL FIELDS, lots of different possibilities:   MD is mismatch info string, NM is edit distance to reference
#       (for info on MD field format see SAM manual footnote, and 
#        sam_MD_field_examples_*.txt files in experiments/reference_data/aligner_format_info)
#########

CIGAR_TYPES_MATCH = ['=']
CIGAR_TYPES_NOOP = ['S','H','P']
CIGAR_TYPES_MUTATION = ['X','I','D']
CIGAR_TYPES_INTRON = ['N']     # 'N' is for introns, but we shouldn't be paying attention to those for genomic DNA seq
CIGAR_TYPES_UNKNOWN = ['M']
# MAYBE-TODO HTSeq doesn't appear aware of the = and X operations...  http://www-huber.embl.de/users/anders/HTSeq/doc/alignments.html#HTSeq.CigarOperation  - I emailed the author about it, no response

# SAM format flag meanings, from https://samtools.github.io/hts-specs/SAMv1.pdf
FLAG_MEANINGS = {1:    'the read is part of a pair', 
                 2:    'the read is part of a pair that properly aligned as a pair', 
                 4:    'unmapped', 
                 8:    'the read is part of a pair and the other mate in the pair is unmapped', 
                 16:   'alignment is reverse-complement', 
                 32:   'the other mate in the pair aligned reverse-complement', 
                 64:   'the read is mate 1 in a pair', 
                 128:  'the read is mate 2 in a pair', 
                 256:  'secondary alignment', 
                 512:  'not passing filters, such as platform/vendor quality controls', 
                 1024: 'PCR or optical duplicate', 
                 2048: 'supplementary alignment (whatever that is)'
                }


class DeepseqError(Exception):
    """ Exception in this module; no special behavior."""
    pass

### Parsing two fastq files in parallel (for paired-end deepseq data)

def parse_2fastq_parallel(file1, file2):
    """ Parse two fastq files in parallel - generator yielding (name, seq1, seq2, qual1, qual2) tuples.

    Doesn't check that the readnames match.
    """
    from Bio.SeqIO.QualityIO import FastqGeneralIterator    # Bio is the biopython package
    with open(file1) as INFILE1:
        with open(file2) as INFILE2:
            generator1 = FastqGeneralIterator(INFILE1)
            generator2 = FastqGeneralIterator(INFILE2)
            if_finished_1, if_finished_2 = False, False
            while True:
                try:                    name1, seq1, qual1 = next(generator1)
                except StopIteration:   if_finished_1 = True
                try:                    name2, seq2, qual2 = next(generator2)
                except StopIteration:   if_finished_2 = True
                name = name1.split()[0]
                if not if_finished_1 and not if_finished_2:
                    yield (name, seq1, seq2, qual1, qual2)
                elif if_finished_1 and if_finished_2:
                    return
                else:
                    raise DeepseqError("One file finished but the other one didn't! Read name %s"%(
                                                                        name if if_finished_2 else name2.split()[0]))
    # TODO unit-test!


def parse_fastx_sam_parallel(fastx_infile, sam_infile):
    """ Parse fastx and resulting sam file in parallel - generator yielding (name, seq, alignment_list) tuples.

    The sam file may contain multiple alignments per read.  Program checks that the readnames match.
    """
    fastx_generator = basic_seq_utilities.name_seq_generator_from_fasta_fastq(fastx_infile)
    sam_generator = iter(HTSeq.bundle_multiple_alignments(HTSeq.SAM_Reader(sam_infile)))
    if_finished_fastx, if_finished_sam = False, False
    while True:
        try:                    name, seq = next(fastx_generator)   # used to be generator.next() in python2
        except StopIteration:   if_finished_fastx = True
        try:                    alns = next(sam_generator)
        except StopIteration:   if_finished_sam = True
        # if both finished, good, we're done
        if if_finished_fastx and if_finished_sam:
            return
        # if one file was finished but the other wasn't, error!
        elif if_finished_fastx or if_finished_sam:
            raise DeepseqError("Parsing seq/aln files in parallel - inconsistent finished states! "
                              +"(If finished: %s %s, %s %s)"%(fastx_infile, if_finished_fastx, sam_infile, if_finished_sam))
        # if all the files still contained data, yield it
        else:
            name = name.split()[0]
            name2 = alns[0].read.name.split()[0]
            if not name2 == name:
                raise DeepseqError("Non-matching readnames between files! %s in %s, %s in %s"%(fastx_infile, name, 
                                                                                               sam_infile, name2))
            yield (name, seq, alns)


### Getting mutation counts from various SAM alignment format fields, as read by HTSeq

def _get_HTSeq_optional_field_either_version(val_or_tuple):
    """ Different HTSeq versions return either val or (name,val) from aln.optional_field(name) - convert either to val. """
    if isinstance(val_or_tuple, tuple): return val_or_tuple[1]
    else:                               return val_or_tuple


def get_HTSeq_optional_field(HTSeq_alignment, field_name):
    """ Return value of optional field (like NM, XM, etc). """
    return _get_HTSeq_optional_field_either_version(HTSeq_alignment.optional_field(field_name))


def check_mutation_count_by_CIGAR_string(HTSeq_alignment, treat_unknown_as='unknown', ignore_introns=False):
    """ Return number of mutations in HTSeq_alignment, based on CIGAR string; -1 if unknown ('M') by default.
    If treat_unknown_as is 'unknown', return -1 whenever an unknown (M, may be match or mismatch) operation is found; 
     if treat_unknown_as is 'mutation' or 'match', count unknowns accordingly.  Return -1 if read is unaligned.
    If ignore_introns is False, count introns (N) as mutations; otherwise don't."""
    global CIGAR_TYPES_MUTATION, CIGAR_TYPES_INTRON, CIGAR_TYPES_UNKNOWN
    # just return -1 for unaligned reads
    if HTSeq_alignment.cigar is None:
        return -1
    # figure out whether to consider intron-skipping ('N') as a mutation or not, based on argument
    if ignore_introns:
        # (need this []+ here so it's a copy, not a reference, and modifying it later doesn't modify the original)
        cigar_types_mutation = [] + CIGAR_TYPES_MUTATION
    else:
        cigar_types_mutation = CIGAR_TYPES_MUTATION + CIGAR_TYPES_INTRON   
    # figure out how to treat unknown matches ('M'), based on argument
    if treat_unknown_as=='unknown':
        # (need this []+ here so it's a copy, not a reference, and modifying it later doesn't modify the original)
        cigar_types_unknown = [] + CIGAR_TYPES_UNKNOWN
    elif treat_unknown_as=='mutation':
        cigar_types_mutation += CIGAR_TYPES_UNKNOWN
        cigar_types_unknown = []
    elif treat_unknown_as=='match':
        cigar_types_unknown = []
    else:
        raise ValueError("treat_unknown_as argument value must be 'mutation', 'match' or 'unknown'")
    # count the mutations, return total count (or instantly return -1 on finding an unknonw)
    mutations = 0
    for cigar_op in HTSeq_alignment.cigar:
        if cigar_op.type in cigar_types_mutation:
            mutations += cigar_op.size
        # if there's an unknown, just return -1, no need to count
        elif cigar_op.type in cigar_types_unknown:
            return -1
    return mutations


def check_mutation_count_by_optional_NM_field(HTSeq_alignment, negative_if_absent=True):
    """ Return #errors in HTSeq_alignment, based on optional NM field; -1 or exception if field missing."""
    # for unalign reads NM field is missing - returns -1
    try:                        return get_HTSeq_optional_field(HTSeq_alignment, 'NM')
    except KeyError:    
        if negative_if_absent:  return -1
        else:                   raise DeepseqError("Optional NM field missing in read %s - can't determine #errors!"%HTSeq_alignment.read.name)


def check_mutation_count_by_optional_MD_field(HTSeq_alignment):
    """ Return number of mutations in HTSeq_alignment, based on optional MD field; -1 if unknown (MD field missing)."""
    # for info on MD field format see SAM manual footnote, 
    #   and sam_MD_field_examples_*.txt files in experiments/reference_data/aligner_format_info
    #       basically a number means matches, a letter means a mismatch to reference (or insertion? is that different?), 
    #       letters preceded by ^ mean deletion from the reference
    try:                mutation_string = get_HTSeq_optional_field(HTSeq_alignment, 'MD')
    except KeyError:    return -1
    # for unalign reads MD field is missing - returns -1
    mutation_letters = [c for c in mutation_string if not (c.isdigit() or c=='^')]
    #   (^ is used in describing a mutation but it shouldn't be counted as a separate mutation - only letters count.)
    return len(mutation_letters)


def check_mutation_count_try_all_methods(HTSeq_alignment, treat_unknown_as='unknown', ignore_introns=False):
    """ Return number of mutations in HTSeq_alignment (look at CIGAR string and NM and MD optional fields); -1 if unknown.
    First check the CIGAR string but only accept the answer if there are no unknown ('M') characters; 
     then check the NM and MD fields and return the result if those fields exist.
    If the CIGAR string is ambiguous and neither of the optional fields exist:
     - if treat_unknown_as is 'unknown', return -1
     - if treat_unknown_as is 'mutation' or 'match', return the CIGAR string result with unknowns counted accordingly.
    If ignore_introns is False, count introns (N) in CIGAR string as mutations; otherwise don't .
    Does NOT guarantee returning a sensible value if the CIGAR, NM and MD fields contain inconsistent information.
    """
    mutation_count = check_mutation_count_by_CIGAR_string(HTSeq_alignment, treat_unknown_as='unknown', 
                                                          ignore_introns=ignore_introns)
    if not mutation_count==-1:  
        return mutation_count
    mutation_count = check_mutation_count_by_optional_NM_field(HTSeq_alignment)
    if not mutation_count==-1:  
        return mutation_count
    mutation_count = check_mutation_count_by_optional_MD_field(HTSeq_alignment)
    if not mutation_count==-1:  
        return mutation_count
    if treat_unknown_as=='unknown':     
        return -1
    return check_mutation_count_by_CIGAR_string(HTSeq_alignment, treat_unknown_as=treat_unknown_as, 
                                                          ignore_introns=ignore_introns)


### Other SAM alignment utilities

def aln_read_coverage(HTSeq_alignment):
    """ Given an HTSeq alignment, return (start,end) tuple describing the part of the read it covers. """
    cigar = HTSeq_alignment.cigar
    if cigar[0].type == 'S':        start = cigar[1].query_from
    else:                           start = 0
    if cigar[-1].type == 'S':       end = cigar[-2].query_to
    else:                           end = cigar[-1].query_to
    return start, end


def aln_read_coverage_fraction(HTSeq_alignment, percent_string=False):
    """ Given an HTSeq alignment, return fraction of the read it covers (or percentage string). """
    s, e = aln_read_coverage(HTSeq_alignment)
    fraction = (e-s)/len(HTSeq_alignment.read.seq)
    if percent_string:  return "%.0f%%"%(fraction*100)
    else:               return fraction


def read_coverage_all_alns(alns):
    """ Given a list of HTSeq alignment objects, return (start,end) tuple describing the part of the read covered by any of them. """
    starts, ends = zip(*[aln_read_coverage(a) for a in alns])
    return min(starts), max(ends)


def print_aln_list_info(alns, sort=True):
    """ Return string of useful information about a list of alignments, optionally sorted. """
    if not alns:    return "no alignments!"
    aln_data = [('-' if a.not_primary_alignment else '+', '1' if a.pe_which=='first' else '2',   
                aln_read_coverage_fraction(a, True), a.optional_field('NM'), 
                a.iv.chrom, a.iv.start, a.iv.end) 
               for a in alns]
    if sort: aln_data.sort(key=lambda x: (x[0], x[1], -int(x[2][:-1]), x[3]))
    # TODO make this a nicer string instead of just str() of a list!
    return alns[0].read.name + '\n' + '\n'.join(str(x) for x in aln_data)


def find_best_aln(alns, min_coverage=0.8, bad_coverage=.5, max_errors=.02, bad_errors=0.05):
    """ If there's exactly one alignment that meets min_coverage and max_errors 
        while all other alignments are bad_coverage and bad_errors or worse, return that; otherwise return None.
    """
    acceptable_alns = [a for a in alns if aln_read_coverage_fraction(a) > bad_coverage 
                       and check_mutation_count_try_all_methods(a)/len(a.read.seq) < bad_errors]
    if len(acceptable_alns) != 1:   
        return None
    a = acceptable_alns[0]
    if aln_read_coverage_fraction(a) > min_coverage and check_mutation_count_try_all_methods(a)/len(a.read.seq) < max_errors:
        return a
    else:
        return None


def primary_or_best_aln(alns, min_coverage=0.8, bad_coverage=.5, max_errors=.02, bad_errors=0.05):
    """ Return primary alignment if there is one, or the only alignment if there's only one, otherwise same as find_best_aln. """
    if len(alns) == 1:      return alns[0]
    primary = [a for a in alns if not a.not_primary_alignment]
    if len(primary) == 1:   return primary[0]
    elif len(primary) > 1:  print("Multiple primary alignments for %s - shouldn't happen!"%a.read.name)
    return find_best_aln(alns, min_coverage, bad_coverage, max_errors, bad_errors)


def interpret_flags(flag_number, details=False):
    """ Given an overall flag number (like 35), split it into bits (32, 2, 1), optionally with verbal interpretations. """
    flags = []
    for i, bit in enumerate(reversed(bin(flag_number))):
        if bit=='1':    flags.append(2**i)
        elif bit=='0':  pass
        elif bit=='b':  break
        else:           print('This binary interpretation has weird bits! %s -> %s'%(flag_number, bin(flag_number)))
    if details:         
        for f in flags:
            print('%s - %s'%(f, FLAG_MEANINGS[f]))
    return flags



################## unit tests ################## 

class Fake_deepseq_objects:
    """ Fake deepseq data objects for testing. """
    # NOTE: not all of those are used in the unit-tests for this module, but they're also imported elsewhere!

    class Fake_HTSeq_cigar_op:
        """ Fake CIGAR operation, mimicking HTSeq cigar object."""
        size = 1
        def __init__(self,string):  
            self.type = string

    class Fake_HTSeq_genomic_pos:
        """ Fake HTSeq.GenomicPosition. """
        def __init__(self, chrom, strand, start, end):
            self.chrom = chrom
            self.strand = strand
            self.start = start
            self.end = end

    class Fake_HTSeq_read:
        """ Fake read, as in HTSeq_alignment.read. """
        def __init__(self,seq='AAA',name='test'):
            self.seq = seq
            self.name = name

    class Fake_HTSeq_alignment:
        """ Fake HTSeq.Alignment object."""

        def __init__(self, seq='AAA', readname='test', unaligned=False, pos=('chr_','+',0,0), 
                     cigar_string=None, optional_field_data={}):    
            self.read = Fake_deepseq_objects.Fake_HTSeq_read(seq,readname)
            if unaligned:
                self.aligned = False
                self.iv = None
            else:
                self.aligned = True
                self.iv = Fake_deepseq_objects.Fake_HTSeq_genomic_pos(*pos)
            self.optional_field_data = optional_field_data
            if cigar_string is None:    self.cigar = None
            else:                       self.cigar = [Fake_deepseq_objects.Fake_HTSeq_cigar_op(c) for c in cigar_string]

        def optional_field(self,field):             
            return self.optional_field_data[field]


class Testing(unittest.TestCase):
    """ Unit-tests for all the functions/classes in this module. """

    def test__parse_fastx_sam_parallel(self):
        output = list(parse_fastx_sam_parallel('_test_inputs/test_parallel2.fq', '_test_inputs/test_parallel2.sam'))
        assert len(output) == 3
        assert [len(x) for x in output] == [3, 3, 3]
        assert [len(x[2]) for x in output] == [2, 1, 1]
        assert output[0][2][0].read.name == output[0][0] == 'ROCKFORD:4:1:1680:975#0/1'
        # this .decode() is needed because the first one is a bytes type: b'ACTAATACGCGGCCTGGAGCTGGACGTTGGAACCAA'
        assert output[0][2][0].read.seq.decode() == output[0][1] ==  'ACTAATACGCGGCCTGGAGCTGGACGTTGGAACCAA'
        # the generator isn't really run until you ask for its results, so I have to run list on it to get the error
        self.assertRaises(DeepseqError, list, parse_fastx_sam_parallel('_test_inputs/test.fq', '_test_inputs/test_parallel2.sam'))

    def test__check_mutation_count_by_CIGAR_string(self):
        # no alignment (CIGAR is None)
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment()
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == -1
        # CIGAR is unambiguous, no MD or NM given (or needed)
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='==')
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == 0
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='XX')
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == 2
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='DD')
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == 2
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='S=')
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == 0
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='SX')
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == 1
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='N=')
        assert check_mutation_count_by_CIGAR_string(fake_alignment, ignore_introns=True) == 0
        assert check_mutation_count_by_CIGAR_string(fake_alignment, ignore_introns=False) == 1
        # CIGAR is ambiguous (contains M's) - return -1, 2 or 0 depending on what treat_unknown_as is set to
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='MM')
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='unknown') == -1
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='mutation') == 2
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='match') == 0
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='M=')
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='unknown') == -1
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='mutation') == 1
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='match') == 0
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='MX')
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='unknown') == -1
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='mutation') == 2
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='match') == 1

    def test__check_mutation_count_by_optional_NM_field(self):
        """ the tested function should return -1 if no NM field, otherwise return value of NM field. """
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment()
        assert check_mutation_count_by_optional_NM_field(fake_alignment, negative_if_absent=True) == -1
        self.assertRaises(DeepseqError, check_mutation_count_by_optional_NM_field, fake_alignment, negative_if_absent=False)
        for x in range(10):
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(optional_field_data={'NM':x})
            assert check_mutation_count_by_optional_NM_field(fake_alignment) == x

    def test__check_mutation_count_by_optional_MD_field(self):
        """ see ~/experiments/reference_data/aligner_format_info/* files for MD field examples."""
        fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment({})
        assert check_mutation_count_by_optional_MD_field(fake_alignment) == -1
        for s in [str(x) for x in range(30)]:
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(optional_field_data={'MD': s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 0
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(optional_field_data={'MD': s+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 0
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(optional_field_data={'MD': s+'A'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 1
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(optional_field_data={'MD': s+'A0G'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 2
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(optional_field_data={'MD': s+'A2G'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 2
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(optional_field_data={'MD': s+'A2G2T2C2N'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 5
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(optional_field_data={'MD': s+'^A'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 1
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(optional_field_data={'MD': s+'^AGC'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 3

    def test__check_mutation_count_try_all_methods(self):
        """ The order of check is CIGAR, NM, MD; CIGAR is skipped if ambiguous; NM and MD skipped if inexistent. 
        Not attempting to deal with inconsistent states sensibly."""
        # all measures agree there are no mutations (with 0-2 of NM/MD fields present)
        for opt_data in [{'NM':0, 'MD':'10'}, {'NM':0}, {'MD':'10'}, {}]:
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='='*10, optional_field_data=opt_data)
            assert check_mutation_count_try_all_methods(fake_alignment) == 0
        # all measures agree there is a mutation (with 0-2 of NM/MD fields present)
        for opt_data in [{'NM':1, 'MD':'A9'}, {'NM':1}, {'MD':'A9'}, {}]:
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='X'+'='*9,optional_field_data=opt_data)
            assert check_mutation_count_try_all_methods(fake_alignment) == 1
        # CIGAR is ambiguous, there are no mutations according to NM/MD (NM, MD or both are present)
        for opt_data in [{'NM':0, 'MD':'10'}, {'NM':0}, {'MD':'10'}]:
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='M'*10, optional_field_data=opt_data)
            assert check_mutation_count_try_all_methods(fake_alignment) == 0
        # CIGAR is ambiguous, there is a  mutation according to NM/MD (NM, MD or both are present)
        for opt_data in [{'NM':1, 'MD':'A9'}, {'NM':1}, {'MD':'A9'}]:
            fake_alignment = Fake_deepseq_objects.Fake_HTSeq_alignment(cigar_string='M'*10, optional_field_data=opt_data)
            assert check_mutation_count_try_all_methods(fake_alignment) == 1

    def test__interpret_flags(self):
        for i in range(1,12):
            self.assertEqual(interpret_flags(2**i), [2**i])
            self.assertEqual(interpret_flags(2**i + 1), [1,2**i])
            self.assertEqual(interpret_flags(2**i - 1), [2**j for j in range(i)])
        self.assertEqual(interpret_flags(1), [1])
        self.assertEqual(interpret_flags(3), [1,2])
        self.assertEqual(interpret_flags(4), [4])
        # and the details option I just tested by hand


if __name__ == "__main__":
    """ Allows both running and importing of this file. """
    print("*** This is a module to be imported to other files - running the built-in test suite. ***")
    unittest.main()

