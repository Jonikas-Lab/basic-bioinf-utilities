#! /usr/bin/env python

"""
Given any number of fasta or fastq files, check if any of the sequences repeat (in one file or between files); depending on options, also check whether any sequence is a subsequence of any other, and whether any reverse-complement sequence is identical to or a subsequence of any of the other sequences.
 -- Weronika Patena, 2012
USAGE: seq_uniqueness_check.py [options] infile1 [infile2 [infile3 ...]]
"""

# standard library
from __future__ import division
import sys, os
from collections import defaultdict
from itertools import combinations
import unittest
# other packages
from Bio import SeqIO
# my modules
from basic_seq_utilities import check_fasta_fastq_format, reverse_complement


def define_option_parser():
    """ Populates and returns an optparse option parser object, with __doc__ as the usage string."""
    from optparse import OptionParser
    parser = OptionParser(__doc__)

    ### test options
    parser.add_option('-t','--test_functionality', action='store_true', default=False, 
                      help="Run the built-in unit test suite (ignores all other options/arguments; default %default).")
    parser.add_option('-T','--test_run', action='store_true', default=False, 
                      help="Run on a test input file, check output against reference files. "
                          + "Ignores all other options/arguments. (default %default).")

    ### functionality options
    parser.add_option('-i','--exact_identity_only', action='store_true', default=False,
                      help="Don't check for seqA being a subsequence of seqB, only for exact identity between sequences " 
                      +"(default %default).")
    parser.add_option('-f','--forward_only', action='store_true', default=False,
                      help="Don't check reverse-complement versions of sequences (default %default).")
    parser.add_option('-e','--ignore_empty_sequences', action='store_true', default=False,
                      help="Ignore empty sequences (i.e. of length 0) (default %default).")
    # MAYBE-TODO could have options governing files - only check for repeats between files, or within one file, etc
    return parser


def check_pair(seqA, seqB, nameA='', nameB='', exact_identity_only=False, forward_only=False, ignore_empty=True):
    """ Check whether seqA and seqB are identical, or one contains the other; same for reverse-complemet; print info. """
    # MAYBE-TODO should empty sequences be allowed?
    if not (seqA and seqB) and ignore_empty:
        return
    if seqA==seqB:              return "%s and %s are identical (%s)"%(nameA, nameB, seqA)
    if not exact_identity_only:
        if seqA in seqB:        return "%s contains %s (%s, %s)"%(nameB, nameA, seqB, seqA)
        if seqB in seqA:        return "%s contains %s (%s, %s)"%(nameA, nameB, seqA, seqB)
    if not forward_only:
        seqBrc = reverse_complement(seqB)
        if seqA==seqBrc:        return "%s and reverse-complement of %s are identical (%s)"%(nameA, nameB, seqA)
        if not exact_identity_only:
            if seqA in seqBrc:  return "reverse-complement of %s contains %s (%s, %s)"%(nameB, nameA, seqBrc, seqA)
            if seqBrc in seqA:  return "%s contains reverse-complement of %s (%s, %s)"%(nameA, nameB, seqA, seqBrc)
    return


def main(infiles, args):
    """ Run the main functionality of the module (see module docstring for more information), excluding testing.
    The options argument should be generated by an optparse parser.
    """

    if not infiles:
        parser.print_help()
        sys.exit("\nError: at least one infile and exactly one outfile are required!")

    all_names_and_seqs = []

    for infile in infiles:
        seq_format = check_fasta_fastq_format(infile)

        with open(infile) as INFILE:
            for sequence in SeqIO.parse(INFILE, seq_format): 
                # using seq.tostring() to convert Biopython Seq objects to plain strings - Seq objects aren't hashable correctly
                all_names_and_seqs.append((sequence.name, sequence.seq.tostring()))

    no_repeats = True
    for (nameA,seqA), (nameB,seqB) in combinations(all_names_and_seqs, 2):
        result = check_pair(seqA, seqB, nameA, nameB, 
                            options.exact_identity_only, options.forward_only, options.ignore_empty_sequences)
        if result:  
            print result
            no_repeats = False
    if no_repeats:
        print "NO REPEATS."


def do_test_run():
    """ Test run: run script on test infile, compare output to reference file."""
    from testing_utilities import run_functional_tests
    test_folder = "test_data"
    sys.exit("NO TESTS DEFINED!")
    # tests in (testname, [test_description,] arg_and_infile_string) format
    test_runs = [ ]
    # argument_converter converts (parser,options,args) to the correct argument order for main
    argument_converter = lambda parser,options,args: (args, options)
    # use my custom function to run all the tests, auto-detect reference files, compare to output.
    return run_functional_tests(test_runs, define_option_parser(), main, test_folder, 
                                argument_converter=argument_converter, append_to_outfilenames='.txt') 
    # LATER-TODO add run-tests!


class Testing(unittest.TestCase):
    """ Unit-tests this module. """

    def test__check_pair(self):
        # Note: I'm not checking the formatting of the return messages, only whether there is one (match) or not.
        # identical sequences:
        assert check_pair('aa','aa', exact_identity_only=True, forward_only=True)
        assert not check_pair('aa','a', exact_identity_only=True, forward_only=False)
        assert not check_pair('aa','aaa', exact_identity_only=True, forward_only=False)
        assert not check_pair('aa','tt', exact_identity_only=False, forward_only=True)
        assert not check_pair('aa','ga', exact_identity_only=False, forward_only=False)
        # seqA contains seqB, or vice versa (but overlaps don't count!)
        assert check_pair('cta', 'ct', exact_identity_only=False, forward_only=True)
        assert check_pair('ct', 'cta', exact_identity_only=False, forward_only=True)
        assert not check_pair('ctc', 'cta', exact_identity_only=False, forward_only=True)
        assert not check_pair('cta', 'tc', exact_identity_only=False, forward_only=True)
        assert not check_pair('cta', 'ag', exact_identity_only=False, forward_only=True)
        # reverse-complement identical 
        assert check_pair('acc','ggt', exact_identity_only=True, forward_only=False)
        assert check_pair('ggt','acc', exact_identity_only=True, forward_only=False)
        assert not check_pair('acc','a', exact_identity_only=True, forward_only=False)
        assert not check_pair('acc','cca', exact_identity_only=True, forward_only=False)
        assert not check_pair('acc','tgg', exact_identity_only=True, forward_only=False)
        # reverse-complement contains
        assert check_pair('acct', 'agg', exact_identity_only=False, forward_only=False)
        assert check_pair('agg', 'acct', exact_identity_only=False, forward_only=False)
        assert not check_pair('acct', 'agg', exact_identity_only=False, forward_only=True)
        assert not check_pair('agg', 'acct', exact_identity_only=False, forward_only=True)
        assert not check_pair('acct', 'agg', exact_identity_only=True, forward_only=False)
        assert not check_pair('agg', 'acct', exact_identity_only=True, forward_only=False)
        assert not check_pair('acct', 'gga', exact_identity_only=False, forward_only=False)
        assert not check_pair('acct', 'tcc', exact_identity_only=False, forward_only=False)
        # ignoring empty sequences are counted or not:
        assert check_pair('','', exact_identity_only=True, forward_only=True, ignore_empty=False)
        assert not check_pair('','', exact_identity_only=True, forward_only=True, ignore_empty=True)
        assert check_pair('aa','', exact_identity_only=False, forward_only=True, ignore_empty=False)
        assert not check_pair('aa','', exact_identity_only=True, forward_only=True, ignore_empty=False)
        assert not check_pair('aa','', exact_identity_only=False, forward_only=True, ignore_empty=True)
    # LATER-TODO more unit-tests?


if __name__=='__main__':
    parser = define_option_parser()
    options,args = parser.parse_args()

    # if run with -t option, do unit tests and quit
    if options.test_functionality:
        print("*** You used the -t option - ignoring all other options/arguments, running the built-in test suite. ***")
        # to run tests for another file, have to use TextTestRunner, not unittest.main -  make a test suite with 
        #   autodetection of all tests (see http://docs.python.org/library/unittest.html#unittest.TestLoader)
        #print("\n * unit-tests for the ______ module")
        #test_suite_1 = unittest.defaultTestLoader.loadTestsFromModule(______
        #unittest.TextTestRunner(verbosity=1).run(test_suite_1)
        # to run tests for current module, just run unittest.main, passing it only the filename 
        #   (by default it takes all of sys.argv and complains about options/arguments it can't recognize)
        print("\n * unit-tests for this module (%s)"%sys.argv[0])
        unittest.main(argv=[sys.argv[0]])   # unittest.main automatically runs sys.exit()

    if options.test_run:
        print("*** You used the -T option - ignoring all other options and running the built-in example test runs. ***")
        test_result = do_test_run()
        sys.exit(test_result)

    # otherwise pass the arguments to the main function
    main(args, options)
