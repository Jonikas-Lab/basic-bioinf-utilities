#! /usr/bin/env python

""" Various basic deepseq-related utilities I wrote - see separate function docstrings.  For importing only.
 --Weronika Patena, 2011
"""

# basic libraries
import unittest


################## fasta/fastq (raw data) utilities ################## 

def parse_fastq(infile):
    """ Given a fastq file, yield successive (header,sequence,quality) tuples. """
    with open(infile) as INFILE:
        while True:
            header = INFILE.next().strip()
            try:
                seq = INFILE.next().strip()
                header2 = INFILE.next().strip()
                qual = INFILE.next().strip()
            except (StopIteration):
                raise Exception("Input FastQ file is malformed! Last record didn't have all four lines!")

            if not header[0]=='@':  
                raise Exception("Malformed input fastq file! Expected seq-header line (@ start), found %s"%header)
            if not header2[0]=='+':  
                raise Exception("Malformed input fastq file! Expected qual-header line (+ start), found %s"%header2)
            header,header2 = header[1:],header2[1:]
            if not (len(header2)==0 or header==header2):   
                raise Exception("Malformed input fastq file! Qual-header %s doesn't match seq-header %s"%(header2,header))
            if not len(seq)==len(qual):             
                raise Exception("Malformed input fastq file! Seq length doesn't match qual length (%s,%s)"%(seq, qual))

            yield (header, seq, qual)

    # TODO add unit-test?


def get_seq_count_from_collapsed_header(header, return_1_on_failure=False):
    """ Given a sequence header from fastx_collapser, return the original sequence count ('>1-243' means 243 sequences).
    If cannot parse the header, exits with an error message, unless return_1_on_failure is True (then returns 1). """
    try:                        header_fields = header.split('-')
    except AttributeError:      header_fields = [header]
    if len(header_fields) > 1:
        try:                    return int(header_fields[-1])
        except ValueError:      pass
    # if didn't find a '-' to split on, or couldn't get an int from the string, either return 1 or fail
    if return_1_on_failure: 
        return 1
    else:                   
        raise ValueError("Can't parse header %s to get original pre-fasxt_collapser sequence count!"%header)
    

################## aligned data utilities ################## 

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
# MAYBE-TODO HTSeq doesn't appear aware of the = and X operations...  http://www-huber.embl.de/users/anders/HTSeq/doc/alignments.html#HTSeq.CigarOperation  - I emailed the author about it

### Getting mutation counts from various SAM alignment format fields, as read by HTSeq

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


def check_mutation_count_by_optional_NM_field(HTSeq_alignment):
    """ Return number of mutations in HTSeq_alignment, based on optional NM field; -1 if unknown (NM field missing)."""
    # for unalign reads NM field is missing - returns -1
    try:                return HTSeq_alignment.optional_field('NM')
    except KeyError:    return -1


def check_mutation_count_by_optional_MD_field(HTSeq_alignment):
    """ Return number of mutations in HTSeq_alignment, based on optional MD field; -1 if unknown (MD field missing)."""
    # for info on MD field format see SAM manual footnote, 
    #   and sam_MD_field_examples_*.txt files in experiments/reference_data/aligner_format_info
    #       basically a number means matches, a letter means a mismatch to reference (or insertion? is that different?), 
    #       letters preceded by ^ mean deletion from the reference
    try:                mutation_string = HTSeq_alignment.optional_field('MD')
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


################## unit tests ################## 

class Testing(unittest.TestCase):
    """ Unit-tests for all the functions/classes in this module. """

    class Fake_cigar_op:
        """ Fake CIGAR operation, mimicking HTSeq cigar object."""
        def __init__(self,string):  self.type = string
        size = 1
    class Fake_alignment:
        """ Fake alignment with optional_field lookup function, mimicking HTSeq alignment object."""
        def __init__(self,data):        self.data = data
        def optional_field(self,x):     return self.data[x]

    def test__check_mutation_count_by_CIGAR_string(self):
        # no alignment (CIGAR is None)
        class fake_alignment:
            cigar = None
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == -1
        # CIGAR is unambiguous, no MD or NM given (or needed)
        class fake_alignment:
            cigar = [self.Fake_cigar_op('='),self.Fake_cigar_op('=')]
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == 0
        class fake_alignment:
            cigar = [self.Fake_cigar_op('X'),self.Fake_cigar_op('X')]
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == 2
        class fake_alignment:
            cigar = [self.Fake_cigar_op('D'),self.Fake_cigar_op('D')]
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == 2
        class fake_alignment:
            cigar = [self.Fake_cigar_op('S'),self.Fake_cigar_op('=')]
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == 0
        class fake_alignment:
            cigar = [self.Fake_cigar_op('S'),self.Fake_cigar_op('X')]
        assert check_mutation_count_by_CIGAR_string(fake_alignment) == 1
        class fake_alignment:
            cigar = [self.Fake_cigar_op('N'),self.Fake_cigar_op('=')]
        assert check_mutation_count_by_CIGAR_string(fake_alignment, ignore_introns=True) == 0
        assert check_mutation_count_by_CIGAR_string(fake_alignment, ignore_introns=False) == 1
        # CIGAR is ambiguous (contains M's) - return -1, 2 or 0 depending on what treat_unknown_as is set to
        class fake_alignment:
            cigar = [self.Fake_cigar_op('M'),self.Fake_cigar_op('M')]
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='unknown') == -1
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='mutation') == 2
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='match') == 0
        class fake_alignment:
            cigar = [self.Fake_cigar_op('M'),self.Fake_cigar_op('=')]
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='unknown') == -1
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='mutation') == 1
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='match') == 0
        class fake_alignment:
            cigar = [self.Fake_cigar_op('M'),self.Fake_cigar_op('X')]
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='unknown') == -1
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='mutation') == 2
        assert check_mutation_count_by_CIGAR_string(fake_alignment, treat_unknown_as='match') == 1

    def test__check_mutation_count_by_optional_NM_field(self):
        """ the tested function should return -1 if no NM field, otherwise return value of NM field. """
        fake_alignment = self.Fake_alignment({})
        assert check_mutation_count_by_optional_NM_field(fake_alignment) == -1
        for x in range(10):
            fake_alignment = self.Fake_alignment({'NM':x})
            assert check_mutation_count_by_optional_NM_field(fake_alignment) == x

    def test__check_mutation_count_by_optional_MD_field(self):
        """ see ~/experiments/reference_data/aligner_format_info/* files for MD field examples."""
        fake_alignment = self.Fake_alignment({})
        assert check_mutation_count_by_optional_MD_field(fake_alignment) == -1
        for s in [str(x) for x in range(30)]:
            fake_alignment = self.Fake_alignment({'MD': s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 0
            fake_alignment = self.Fake_alignment({'MD': s+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 0
            fake_alignment = self.Fake_alignment({'MD': s+'A'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 1
            fake_alignment = self.Fake_alignment({'MD': s+'A0G'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 2
            fake_alignment = self.Fake_alignment({'MD': s+'A2G'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 2
            fake_alignment = self.Fake_alignment({'MD': s+'A2G2T2C2N'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 5
            fake_alignment = self.Fake_alignment({'MD': s+'^A'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 1
            fake_alignment = self.Fake_alignment({'MD': s+'^AGC'+s })
            assert check_mutation_count_by_optional_MD_field(fake_alignment) == 3

    def test__check_mutation_count_try_all_methods(self):
        """ The order of check is CIGAR, NM, MD; CIGAR is skipped if ambiguous; NM and MD skipped if inexistent. 
        Not attempting to deal with inconsistent states sensibly."""
        # all measures agree there are no mutations (with 0-2 of NM/MD fields present)
        for optional_data in [{'NM':0, 'MD':'10'}, {'NM':0}, {'MD':'10'}, {}]:
            fake_alignment = self.Fake_alignment(optional_data)
            fake_alignment.cigar = [self.Fake_cigar_op('=') for x in range(10)]
            assert check_mutation_count_try_all_methods(fake_alignment) == 0
        # all measures agree there is a mutation (with 0-2 of NM/MD fields present)
        for optional_data in [{'NM':1, 'MD':'A9'}, {'NM':1}, {'MD':'A9'}, {}]:
            fake_alignment = self.Fake_alignment(optional_data)
            fake_alignment.cigar = [self.Fake_cigar_op('X')] + [self.Fake_cigar_op('=') for x in range(9)]
            assert check_mutation_count_try_all_methods(fake_alignment) == 1
        # CIGAR is ambiguous, there are no mutations according to NM/MD (NM, MD or both are present)
        for optional_data in [{'NM':0, 'MD':'10'}, {'NM':0}, {'MD':'10'}]:
            fake_alignment = self.Fake_alignment(optional_data)
            fake_alignment.cigar = [self.Fake_cigar_op('M') for x in range(10)]
            assert check_mutation_count_try_all_methods(fake_alignment) == 0
        # CIGAR is ambiguous, there is a  mutation according to NM/MD (NM, MD or both are present)
        for optional_data in [{'NM':1, 'MD':'A9'}, {'NM':1}, {'MD':'A9'}]:
            fake_alignment = self.Fake_alignment(optional_data)
            fake_alignment.cigar = [self.Fake_cigar_op('M') for x in range(10)]
            assert check_mutation_count_try_all_methods(fake_alignment) == 1

    def test__get_seq_count_from_collapsed_header(self):
        for bad_header in ['aaa','aaa-aa', 'a-3-a', 'a-3a', '3-a','a+3','a:3','a 3','3',3,0,100,None,True,False,[],{}]:
            assert get_seq_count_from_collapsed_header(bad_header, return_1_on_failure=True) == 1
            self.assertRaises(ValueError, get_seq_count_from_collapsed_header, bad_header, return_1_on_failure=False) 
        for header_prefix in ['a','a ','>a','> a','a b c','a-b-c','a-3-100']:
            for count in [0,1,3,5,123214]:
                assert get_seq_count_from_collapsed_header(header_prefix+'-'+str(count)) == count


if __name__ == "__main__":
    """ Allows both running and importing of this file. """
    print("*** This is a module to be imported to other files - running the built-in test suite. ***")
    unittest.main()

