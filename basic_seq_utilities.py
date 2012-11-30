#! /usr/bin/env python

"""
Various basic DNA utilities I wrote (parse_fasta, complement, reverse_complement, find_seq_between, generate_seq_slices, ...) - most of them just imported from their own files, since I need them as separate programs too.
Weronika Patena, 2008-2010
"""

### basic library
import unittest
import sys
import re
import random
### my modules
# help functions
from parse_fasta import parse_fasta,print_seq
from read_input import read_input
from transform_sequence_input import transform_sequence_input
# independently useful functions
from complement import complement
from reverse_complement import reverse_complement
from reverse import reverse

################## global constants ################## 

SEQ_ENDS = ['5prime','3prime']
SEQ_STRANDS = ['+','-']
SEQ_DIRECTIONS = ['forward','reverse']
SEQ_ORIENTATIONS = ['sense','antisense']

# the standard extensions for fasta/fastq files (basic form, processed later)
FASTA_EXTENSIONS = "fa fas fasta"
FASTQ_EXTENSIONS = "fq fastq"

# Biopython fastq quality encodings (standard Phred+33 one, new Illumina Phred+64 one, old Illumina nonPhred+64)
# Can be used like this:  "for read in SeqIO.parse(INFILE, "fastq-illumina"):"
FASTQ_QUALITY_ENCODINGS = ["fastq-sanger", "fastq-illumina", "fastq-solexa"]

# colors to use for different chromosome types, by default
CHROMOSOME_TYPE_COLORS = { 'chromosome': 'black', 'scaffold': 'blue', 'chloroplast': 'green', 'mitochondrial': 'red', 
                          'cassette': '0.6'}

################## end of global constants ################## 

### processing global constants
FASTA_EXTENSIONS = (FASTA_EXTENSIONS+" "+FASTA_EXTENSIONS.upper()).split(' ')
FASTQ_EXTENSIONS = (FASTQ_EXTENSIONS+" "+FASTQ_EXTENSIONS.upper()).split(' ')
assert not (set(FASTA_EXTENSIONS) & set(FASTQ_EXTENSIONS))

### formatting of biological data

def format_base_distance(x):
    """ Format base distance as a string: 1bp for 1, 12345bp for 12345, 2kb for 2000, 2Mb for 2000000. """
    x = int(x)
    if x==0:                return "0bp"
    elif not x % 1000000:   return "%sMb"%(x/1000000)
    elif not x % 1000:      return "%skb"%(x/1000)
    else:                   return "%sbp"%x

### basic fasta/fastq functions

def write_fasta_line(seqname, seq, OUTFILE=sys.stdout):
    """ Given a name and sequence, print in one-line fasta format to OUTFILE (default STDOUT). """
    OUTFILE.write(">%s\n%s\n"%(seqname, seq))


def name_seq_generator_from_fasta(fasta_infile):
    """ Yield successive (name,seq) pairs read from fasta file. """
    from Bio import SeqIO   # Bio is the biopython package
    with open(fasta_infile) as INFILE:
        for record in SeqIO.parse(INFILE, 'fasta'):
            # record.seq is a biopython Seq object, record.seq.tostring() is a string, which is what we want
            yield record.name, record.seq.tostring()


def name_seq_generator_from_fastq(fastq_infile):
    """ Yield successive (name,seq) pairs read from fastq file (ignores qualities!). """
    from Bio.SeqIO.QualityIO import FastqGeneralIterator    # Bio is the biopython package
    with open(fastq_infile) as INFILE:
        for name,seq,_ in FastqGeneralIterator(INFILE):
            yield name, seq


def check_fasta_fastq_format(infilename, verbose=False):
    """ Check fasta/fastq format based on infilename extension; return "fasta" or "fastq", raise ValueError if neither. 

    Note: the result can be used as the format argument to Bio.SeqIO, but for fastq ONLY if you don't care 
     about the quality encoding (otherwise you should use one of FASTQ_QUALITY_ENCODINGS, not just "fastq").
    """

    import os
    extension = os.path.splitext(infilename)[1][1:]
    if extension in FASTA_EXTENSIONS:       seq_format = "fasta"
    elif extension in FASTQ_EXTENSIONS:     seq_format = "fastq"
    else:
        raise ValueError("File %s has an unknown extension %s! Allowed extensions are fasta (%s) and fastq (%s)."%(
                            infilename, extension, ', '.join(FASTA_EXTENSIONS), ', '.join(FASTQ_EXTENSIONS)))
    if verbose:     
        formatted_output.append("File %s recognized as %s.\n"%(infilename, seq_format))
    return seq_format


def name_seq_generator_from_fasta_fastq(infile, verbose_filetype_detection=False):
    """ Yield successive (name,seq) pairs read from fasta or fastq file (filetype detected by extension). """
    seq_format = check_fasta_fastq_format(infile, verbose=verbose_filetype_detection)
    if seq_format=='fasta':     return name_seq_generator_from_fasta(infile)
    elif seq_format=='fastq':   return name_seq_generator_from_fastq(infile)
    else:                       raise ValueError("Unknown input file format %s (file %s)."%(seq_format, infile))
    

# Note: normally it's probably better to parse fastq using biopython or HTSeq or such!  But it can be convenient to have a simple locally defined function with no dependencies requiring installation, so I'll keep this.
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
        raise ValueError("Can't parse header %s to get original pre-fastx_collapser sequence count!"%header)


### utilities to deal with standard genomic chromosome names (sort them correctly etc)

def chromosome_type(chromosome, separate_other=False):
    """ Given the chromosome name, return type: chromosome/scaffold/chloroplast/mitochondrial/cassette/other.

    Check for chromosome/scaffold/chloroplast/mitochondrial/cassette (by substrings or startswith);
     if a name didn't match any, or matched more than one, return 'other', 
     or if separate_other is True, return the first _-separated word of the original name.
    """
    # first see which words are found in the name (check for all of them, just in case more than one is found)
    types = set()
    if chromosome.startswith('chr'):    types.add('chromosome')
    if 'chromosome' in chromosome:      types.add('chromosome')
    if 'scaffold' in chromosome:        types.add('scaffold')
    if 'chloro' in chromosome:          types.add('chloroplast')
    if 'mito' in chromosome:            types.add('mitochondrial')
    if 'cassette' in chromosome:        types.add('cassette')
    # now if the name matched exactly one type, return that, otherwise return 'other' or the first part of the name.
    if len(types)==1:       return types.pop()
    elif separate_other:    return chromosome.split('_')[0]
    else:                   return 'other'


def chromosome_sort_key(chromosome_name):
    """ Sort key: 1) sort chromosomes, then other, then scaffolds; 2) inside each category sort numerically if possible. 

    Numerical sort sorts 'chr1' before 'chr12' correctly.  Only numbers at the end of the chromosome name are recognized.
    Underscores are stripped from 
    Names are the full 
    The full original chromosome_name is used as the last part of the key, in case of two names resulting in the same key
     (like chr1 and chr_1 and chr01).
    """
    chromosome_data = re.search('^(?P<name>.*[^\d])?(?P<number>\d*)', chromosome_name)
    # chromosome_base is the original chromosome name without the last digits; if there were no non-digit characters, it's None
    chromosome_base = chromosome_data.group('name').strip('_')
    # chromosome_number is based on the number at the end of the chromosome name, or 0 if none was present
    chromosome_number = int(chromosome_data.group('number')) if chromosome_data.group('number') else 0
    if chromosome_base in ('chromosome', 'chr'):    return (1, 'chromosome', chromosome_number, chromosome_name)
    elif chromosome_base=='scaffold':               return (3, 'scaffold', chromosome_number, chromosome_name)
    else:                                           return (2, chromosome_base, chromosome_number, chromosome_name) 


def chromosome_color(chromosome, color_for_other=None):
    """ Return a color by chromosome type: black=chromosome, blue=scaffold, green=chloroplast, red=mitochondrial, grey=cassette.

    Uses chromosome_type to determine type. 
    For other types, raise exception if color_for_other is None, otherwise return color_for_other.
    """
    chr_type = chromosome_type(chromosome)
    try:
        return CHROMOSOME_TYPE_COLORS[chr_type]
    except KeyError:
        if color_for_other is not None:     return color_for_other
        else:                               raise ValueError("Can't parse chromosome name %s!"%chromosome)


### testing whether two sequences contain one another, or overlap (based purely on position, no sequence involved)

def position_test_contains(seq1_start, seq1_end, seq2_start, seq2_end):
    """ Return True if seq2 is completely contained within seq1, False otherwise. All positions should be numbers."""
    # make sure the arguments make sense
    if seq1_start > seq1_end or seq2_start > seq2_end:
        raise ValueError("sequence start positions must not be higher than end positions!")
    # the test is simple
    return bool(seq1_start <= seq2_start and seq2_end <= seq1_end)

def position_test_overlap(seq1_start, seq1_end, seq2_start, seq2_end, return_overlap=False):
    """ Returns True if seq1/seq2 overlap, else False; or the number of bases of overlap if return_overlap is True (min 0).
    Assumes the positions are inclusive (i.e. length = end-start+1). """
    # make sure the arguments make sense
    if seq1_start > seq1_end or seq2_start > seq2_end:
        raise ValueError("sequence start positions must not be higher than end positions!")
    ### if we don't want to measure the overlap, the test is simple
    if return_overlap is False:     return bool(seq2_start <= seq1_end and seq2_end >= seq1_start)
    ### if we do want to measure the overlap, need to do more work:
    # if seq1 contains seq2, or the opposite, return the length of the shorter sequence (see docstring for why +1)
    if position_test_contains(seq1_start, seq1_end, seq2_start, seq2_end):  return seq2_end - seq2_start + 1
    if position_test_contains(seq2_start, seq2_end, seq1_start, seq1_end):  return seq1_end - seq1_start + 1
    # if seq2 starts before seq1, or ends after seq1, the overlap is the distance between the "inner" start and end
    #   (there's a +1 in both cases, because if seq1_end==seq2_start, there's 1bp overlap)
    if seq2_start < seq1_start:                                             return max(seq2_end - seq1_start + 1, 0)  
    if seq2_end > seq1_end:                                                 return max(seq1_end - seq2_start + 1, 0)  


### other position-based functions

def find_seq_between(seq, flankA, flankB, exclusive=True):
    """ Return the fragment of seq between flankA and flankB (first occurence). """
    posA = seq.find(flankA)
    posB = seq[posA:].find(flankB)  # only look for flankB in the region AFTER flankA
    if not posB==-1:    posB += posA        # if found, add the length of the omitted region to result 
    lenA = len(flankA)
    lenB = len(flankB)
    if posA==-1 or posB==-1:    return ''
    elif exclusive:             return seq[posA+lenA:posB]
    else:                       return seq[posA:posB+lenB]

def generate_seq_slices(seq, slice_len, step=1, extra_last_slice=True):
    """ Generate (start_pos,slice_seq) pairs with slices covering seq with the given slice length and step size. 

    Start_pos is zero-based.
    
    When seq doesn't cleanly divide into slices with given step/slice_len (i.e. len(seq)-slice_len isn't divisible by step)
     (for example for 'aaccttgg' with step=3 and slice_len=4):
     - if extra_last_slice is False, ignore the remainder of seq (at most slice_len-1) 
        (just generate 'aacc', 'cttg' - the last 'g' is ignored)
     - if it's True, yield an extra slice from the last slice_len bases of seq, that's less than step distant from the previous slice
        (after 'aacc' and 'cttg' yield 'ttgg' to cover the last 'g', even though the offset between 'cttg' and 'ttgg' is 1, not 3)
    """
    if slice_len <= 0:   raise ValueError("slice_len argument to generate_seq_slices must be a positive integer!")
    if step <= 0:        raise ValueError("step argument to generate_seq_slices must be a positive integer!")
    seqlen = len(seq)
    # if the seq will fit in one slice, just return that
    if slice_len >= seqlen: 
        yield 0, seq
        return
    for pos in range(0, seqlen-slice_len, step):
        yield pos, seq[pos:pos+slice_len]
    # Special case after the last normal slice, if there's "leftover" sequence:
    #  if the sequence isn't evenly divided into correct-length slices, make one last slice covering the last slice_len 
    #   of the seq (the difference between this slice and the previous one will be <step)
    #  But only if the next slice would be inside seq at all!  
    #   if slice<step, for example 'actg' with len 1 and step 2 should be 'a t', never 'a t g'.
    if (pos+slice_len < seqlen) and (pos+step < seqlen) and extra_last_slice:
        yield seqlen-slice_len, seq[-slice_len:]


####################################### Unit-tests #########################################

class Testing_everything(unittest.TestCase):
    """ Testing all functions/classes.etc. """

    def test__format_base_distance(self):
        assert format_base_distance(0) == "0bp"
        for x in [1,10,11,100,999, 1001, 1234, 13987291876]:
            assert format_base_distance(x) == "%sbp"%x
            assert format_base_distance(x*1000) == "%skb"%x
            assert format_base_distance(x*1000000) == "%sMb"%x
        for x in ['a', '1.23', [], [12]]:
            self.assertRaises((ValueError,TypeError), format_base_distance, x)


    def test__parse_fastq(self):
        # need to actually run through the whole iterator to test it - defining it isn't enough
        def parse_fastq_get_first_last(infile):
            seq_iter = parse_fastq(infile)
            seq1 = seq_iter.next()
            for seq in seq_iter:    
                seqN = seq
            return seq1, seqN
        # checking first and last record of _test_inputs/test.fq
        seq1, seqN = parse_fastq_get_first_last("_test_inputs/test.fq")
        assert seq1 == ("ROCKFORD:4:1:1680:975#0/1", "NCTAATACGCGGCCTGGAGCTGGACGTTGGAACCAA", 
                        "BRRRQWVWVW__b_____^___bbb___b_______")
        assert seqN == ("ROCKFORD:4:1:3367:975#0/1", "NCTAAGGCAGATGGACTCCACTGAGGTTGGAACCAA", 
                        "BQQQNWUWUUbbb_bbbbbbbbb__b_bb_____b_") 
        # non-fastq input files
        self.assertRaises(Exception, parse_fastq_get_first_last, "_test_inputs/test.fa")
        self.assertRaises(Exception, parse_fastq_get_first_last, "_test_inputs/textcmp_file1.txt")


    def test__get_seq_count_from_collapsed_header(self):
        for bad_header in ['aaa','aaa-aa', 'a-3-a', 'a-3a', '3-a','a+3','a:3','a 3','3',3,0,100,None,True,False,[],{}]:
            assert get_seq_count_from_collapsed_header(bad_header, return_1_on_failure=True) == 1
            self.assertRaises(ValueError, get_seq_count_from_collapsed_header, bad_header, return_1_on_failure=False) 
        for header_prefix in ['a','a ','>a','> a','a b c','a-b-c','a-3-100']:
            for count in [0,1,3,5,123214]:
                assert get_seq_count_from_collapsed_header(header_prefix+'-'+str(count)) == count


    def test__chromosome_type(self):
        for chrom in ('chromosome_1 chromosome_12 chromosome_FOO chromosome_1_2_3 chromosome1 chr_3 chr4'.split()):
            assert chromosome_type(chrom) == 'chromosome'
        for chrom in ('scaffold_1 scaffold_12 scaffold_FOO scaffold_1_2_3 scaffold2 scaffold'.split()):
            assert chromosome_type(chrom) == 'scaffold'
        for chrom in ('mitochondrial mitochondrion mito_3 mito3'.split()):
            assert chromosome_type(chrom) == 'mitochondrial'
        for chrom in ('chloroplast chloroplast_4 chloroplastABC chloroplast123 chloro_3 chloro3'.split()):
            assert chromosome_type(chrom) == 'chloroplast'
        for chrom in ('cassette insertion_cassette cassette_pMJ0013 insertion_cassette_A'.split()):
            assert chromosome_type(chrom) == 'cassette'
        for chrom in ('123 something random mitochondrial_chromosome chr_chloro mito_and_chloro mito_cassette'.split()):
            assert chromosome_type(chrom, separate_other=False) == 'other'
        assert chromosome_type('123', separate_other=True) == '123'
        assert chromosome_type('12_3', separate_other=True) == '12'
        assert chromosome_type('something', separate_other=True) == 'something'
        assert chromosome_type('mito_and_chloro', separate_other=True) == 'mito'
        assert chromosome_type('foo_bar_baz', separate_other=True) == 'foo'

    def test__chromosome_sort_key(self):
        chroms_sorted = (
            'chr chromosome_1 chr_2 chr03 chr3 chr_3 chromosome_3 chromosome_12 chromosome_101 chromosome_300 chr301'.split()
            +'31a AAA AAA3 cassette chloroplast insertion_cassettes_31_a some_thing something31a'.split()
            +'scaffold02 scaffold2 scaffold_3 scaffold_3 scaffold_21'.split()
        )
        for _ in range(100):
            chroms = list(chroms_sorted)
            random.shuffle(chroms)
            assert sorted(chroms, key=chromosome_sort_key) == chroms_sorted

    def test__chromosome_color(self):
        assert all([chromosome_color(c) == 'black' for c in 'chromosome_1 chromosome chromosome_A'.split()])
        assert all([chromosome_color(c) == 'blue' for c in 'scaffold_1 scaffold scaffold_A'.split()])
        assert all([chromosome_color(c) == 'green' for c in 'chloroplast chloroplast_A'.split()])
        assert all([chromosome_color(c) == 'red' for c in 'mitochondrial mitochondrial_A'.split()])
        assert all([chromosome_color(c) == '0.6' for c in 'cassette insertion_cassette cassette_pMJ013b'.split()])
        for other_chr in 'some_thing foo_bar mito_and_chloro'.split():
            assert chromosome_color(other_chr, color_for_other='cyan') == 'cyan'
            self.assertRaises(ValueError, chromosome_color, other_chr, color_for_other=None)


    def test__position_test_contains(self):
        # raises an error if start>end for either sequence
        self.assertRaises(ValueError, position_test_contains, 11,10,1,1)
        self.assertRaises(ValueError, position_test_contains, 1,1,11,10)
        # seq1 contains seq2 - True
        assert position_test_contains(1,10,2,3)==True
        assert position_test_contains(1,10,9,10)==True
        assert position_test_contains(1,10,10,10)==True
        assert position_test_contains(10,10,10,10)==True
        # seq2 contains seq1 - False
        assert position_test_contains(2,3,1,10)==False
        assert position_test_contains(9,10,1,10)==False
        assert position_test_contains(10,10,1,10)==False
        # seq2 and seq1 overlap - False
        assert position_test_contains(2,3,3,10)==False
        assert position_test_contains(2,4,3,10)==False
        assert position_test_contains(3,10,2,3)==False
        assert position_test_contains(3,10,2,4)==False
        # seq2 and seq1 don't overlap - False
        assert position_test_contains(2,3,4,10)==False
        assert position_test_contains(2,3,6,10)==False
        assert position_test_contains(4,10,2,3)==False
        assert position_test_contains(6,10,2,3)==False


    def test__position_test_overlap__returns_bool(self):
        # raises an error if start>end for either sequence
        self.assertRaises(ValueError, position_test_overlap, 11,10,1,1, False)
        self.assertRaises(ValueError, position_test_overlap, 1,1,11,10, False)
        # seq1 contains seq2 - True
        assert position_test_overlap(1,10,2,3,False)==True
        assert position_test_overlap(1,10,9,10,False)==True
        assert position_test_overlap(1,10,10,10,False)==True
        assert position_test_overlap(10,10,10,10,False)==True
        # seq2 contains seq1 - True
        assert position_test_overlap(2,3,1,10,False)==True
        assert position_test_overlap(9,10,1,10,False)==True
        assert position_test_overlap(10,10,1,10,False)==True
        # seq2 and seq1 overlap - True
        assert position_test_overlap(2,3,3,10,False)==True
        assert position_test_overlap(2,4,3,10,False)==True
        assert position_test_overlap(3,10,2,3,False)==True
        assert position_test_overlap(3,10,2,4,False)==True
        # seq2 and seq1 don't overlap - False
        assert position_test_overlap(2,3,4,10,False)==False
        assert position_test_overlap(2,3,6,10,False)==False
        assert position_test_overlap(4,10,2,3,False)==False
        assert position_test_overlap(6,10,2,3,False)==False

    def test__position_test_overlap__returns_number(self):
        # raises an error if start>end for either sequence
        self.assertRaises(ValueError, position_test_overlap, 11,10,1,1, True)
        self.assertRaises(ValueError, position_test_overlap, 1,1,11,10, True)
        # seq1 contains seq2 - True
        assert position_test_overlap(1,10,2,3,True)==2
        assert position_test_overlap(1,10,9,10,True)==2
        assert position_test_overlap(1,10,10,10,True)==1
        assert position_test_overlap(10,10,10,10,True)==1
        # seq2 contains seq1 - True
        assert position_test_overlap(2,3,1,10,True)==2
        assert position_test_overlap(9,10,1,10,True)==2
        assert position_test_overlap(10,10,1,10,True)==1
        # seq2 and seq1 overlap - True
        assert position_test_overlap(2,3,3,10,True)==1
        assert position_test_overlap(2,4,3,10,True)==2
        assert position_test_overlap(3,10,2,3,True)==1
        assert position_test_overlap(3,10,2,4,True)==2
        # seq2 and seq1 don't overlap - False
        assert position_test_overlap(2,3,4,10,True)==0
        assert position_test_overlap(2,3,6,10,True)==0
        assert position_test_overlap(4,10,2,3,True)==0
        assert position_test_overlap(6,10,2,3,True)==0


    def test__find_seq_between(self):
        pass
    # TODO implement!


    def test__generate_seq_slices(self):
        # arguments to generate_seq_slices are (seq, slice_len, step);
        #  the output is a list of (pos,slice_seq) tuples, 
        #   which I'm mostly switching to [pos_tuple, seq_tuple] with zip for easier typing.
        # testing step==1
        assert zip(*generate_seq_slices('actg',1,1)) == [tuple(range(4)), tuple('a c t g'.split())]
        assert zip(*generate_seq_slices('actg',2,1)) == [tuple(range(3)), tuple('ac ct tg'.split())]
        assert zip(*generate_seq_slices('actg',3,1)) == [tuple(range(2)), tuple('act ctg'.split())]
        for slice_len in (4,5,10,100,12345):
            assert list(generate_seq_slices('actg',slice_len,1)) == [(0, 'actg')]
        # testing step==2, and the extra_last_slice argument
        assert list(generate_seq_slices('actg',1,2)) == [(0,'a'), (2,'t')]
        assert list(generate_seq_slices('actg',2,2)) == [(0,'ac'), (2,'tg')]
        assert list(generate_seq_slices('actg',3,2,False)) == [(0,'act')]
        assert list(generate_seq_slices('actg',3,2,True)) == [(0,'act'),(1,'ctg')]
        for slice_len in (4,5,10,100,12345):
            assert list(generate_seq_slices('actg',slice_len,1)) == [(0,'actg')]
        ### bad inputs - note that just making the generator doesn't raise an exception, only iterating over it with G.next!
        # step can't be 0, negative, or a float (except that it can be a float if seqlen<=slice_len, since it's not used at all then)
        for slice_len in (1,2,3,4,5,10,100,12345):
            for step in (0, -1, -2, -100, -1235):
                G = generate_seq_slices('actg', slice_len, step)
                self.assertRaises(ValueError, G.next)
        for slice_len in (1,2,3):
            for step in (.5, 1.5, 2.1, 9.0001, 4.999, 1235453534.1):
                G = generate_seq_slices('actg', slice_len, step)
                self.assertRaises(TypeError, G.next)
        # slice_len can't be 0 or negative, or a float (but it can be a float if it's >=seqlen, since then it never gets used)
        for step in (1,2,3,4,5,10,100,12345):
            for slice_len in (0, -1, -2, -100, -1314, 0.5, 2.1, -100.1, -0.01):
                G = generate_seq_slices('actg', slice_len, step)
                self.assertRaises((ValueError,TypeError), G.next)


if __name__=='__main__':
    """ If module is ran directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()

