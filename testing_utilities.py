#! /usr/bin/env python

"""
Various testing utilities I wrote (for unit tests, functional tests etc) - see docstring for each function/class for what it does.
 --Weronika Patena, 2012
"""

import itertools
import unittest
import re

debug=0

### output file comparison for functional testing
# The reason for this function is to be able to compare an output file and a reference file but ignore things that are expected to change (like the date in the header) - the reference file can contain lines starting with <REGEX> and then giving a regular expression (.* is fine), in which case the output file lines are compared against the regex instead, which allows wildcards for dates/etc.  See my stackoverflow question http://stackoverflow.com/questions/9726214/testing-full-program-by-comparing-output-file-to-reference-file-whats-it-calle


def _advance_iter_keep_state(iterator, skip_IGNORE_lines=False, skip_empty_lines=False):
    """ Try to advance the iterator; return ('', True) if StopIteration was raised, (next_item, False) otherwise. 
    Skip empty (or whitespace-only) lines and lines starting with <IGNORE> if the appropriate flags are True. 
    """
    stop_iteration = False
    # try to advance the iterator: at least once (that's what first_pass is for), and more times if line should be skipped 
    first_pass = True
    while first_pass or (skip_IGNORE_lines and next_item.startswith('<IGNORE>'))\
                     or (skip_empty_lines and not next_item.strip(' \t\n\r')):
        first_pass = False
        try:
            next_item = iterator.next()
        except StopIteration:
            next_item = ''
            stop_iteration = True
            break
    return next_item, stop_iteration
    # MAYBE-TODO should probably write unit-tests for this too, especially if I end up using it for multiple purposes


def clean_line(curr_line, clean_whitespace=True, make_lowercase=False):
    """ make line lowercase; change between-word whitespace to single spaces, remove start/end whitespace."""
    new_line = curr_line
    if make_lowercase: 
        new_line = new_line.lower()
    # replace all multi-space-tab strings with single space; remove space from start/end
    if clean_whitespace:
        new_line = re.sub('\s+', ' ', new_line)
        new_line = new_line.strip()
    return new_line
    # MAYBE-TODO should probably write unit-tests for this too, especially if I end up using it for multiple purposes


def compare_files_with_regex(iter_1, iter_2, 
                             ignore_empty_lines=False, ignore_whitespace=False, ignore_case=False):
    """ return True if the two string iterators are identical (except for special cases), else description of difference.

    The inputs can be lists of strings, generators, sets, open files, etc.  NO INPUT CHECKING IS DONE.
    The program goes over both iterators in parallel (always looking at element1 and element2), as follows:
     - if neither element starts with either special string (<REGEX> or <IGNORE>), compare them directly:
         if the elements are identical, advance both iterators, otherwise return the two elements
     - if one element starts with <REGEX>, parse the rest of it as a regular expression: 
         (removing whitespace from start/end, and adding ^ and $ to start and end)
         if the other element matches the regex, advance both iterators, otherwise return the two elements
     - if both elements start with <REGEX>, raise an exception - can't compare two regular expressions!
     - if either element starts with <IGNORE>, ignore it and advance its iterator to the next element.
    If the end of both iterators has been reached without finding a mismatch, return True; 
     if the end of one but not the other has been reached, return a string stating that and the first extra line. 

    If ignore_empty_lines is True, all lines containing nothing or just whitespace are ignored (treated like <IGNORE>).
    If ignore_whitespace is True, change all sequences of multiple spaces/tabs to a single space, 
     and remove all space at the beginning and end, before doing comparisons.
    If ignore_case is True, change all elements to lowercase before doing comparisons.
    """
    # MAYBE-TODO add some kind of input-checking? (If I do that, add it to unit-tests too)
    # if the arguments are lists/sets/something, convert them to iterators; this leaves iterators unchanged.
    iter1, iter2 = iter(iter_1), iter(iter_2)

    if debug:   print " *** NEW COMPARISON ***"
    while True:

        # advance both iterators to compare the next line pair (skipping IGNORE lines, and empty lines if requested)
        # if either iterator is stopped, exit immediately
        (line1, iter1_stopped) = _advance_iter_keep_state(iter1,skip_IGNORE_lines=True,skip_empty_lines=ignore_empty_lines)
        (line2, iter2_stopped) = _advance_iter_keep_state(iter2,skip_IGNORE_lines=True,skip_empty_lines=ignore_empty_lines)
        if iter1_stopped or iter2_stopped:  break


        if debug:   print 'raw lines:\n\t1) "%s"\n\t2) "%s"'%(line1, line2)

        # raise exception if both elements are regexes - can't match two regexes!
        if line1.startswith('<REGEX>') and line2.startswith('<REGEX>'):
            raise ValueError("Both elements start with <REGEX>, can't match! %s, %s."%(line1,line2))

        # clean up non-regex lines: clean up whitespace, make lowercase
        for curr_line in line1,line2:
            if not line1.startswith('<REGEX>'):     line1 = clean_line(line1, ignore_whitespace, ignore_case)
            if not line2.startswith('<REGEX>'):     line2 = clean_line(line2, ignore_whitespace, ignore_case)

        if debug:   print 'cleaned lines:\n\t1) "%s"\n\t2) "%s"'%(line1.strip(), line2.strip())

        # if one of the lines is a regex, apply it to the other line; return both lines if they don't match
        flags = re.IGNORECASE if ignore_case else 0
        if line1.startswith('<REGEX>'):
            if not re.match('^%s$'%(line1[7:].strip()), line2, flags=flags):  return (line1.strip(),line2.strip())
        elif line2.startswith('<REGEX>'):
            if not re.match('^%s$'%(line2[7:].strip()), line1, flags=flags):  return (line1.strip(),line2.strip())

        # if neither line is a regex, compare them: return both lines if they don't match
        else:
            if not line1==line2:  return (line1.strip(),line2.strip())

        # if there wasn't a return or exception, just repeat the while loop

    # End condition: if all lines matched and both iterators are empty, return True; 
    #  if one iterator still has non-skipped lines left, return an info line and the next line from the other iterator.
    if iter1_stopped and iter2_stopped: return True
    elif iter1_stopped:                 return ("The first iterator ended. Second iterator next line:\n", line2.strip())
    elif iter2_stopped:                 return ("The second iterator ended. First iterator next line:\n", line1.strip())



############################## unit-tests of the functions in this file ##################################

class Testing__everything(unittest.TestCase):
    """ Testing all functions/classes.etc. """

    def test__compare_files_with_regex(self):
        # MAYBE-TODO in most cases when there shouldn't be a match I'm only checking that the output!=True, not the exact non return value which describes the mismatch - could fix that, but I'm not sure it's worth it.
        # any list should be identical to itself
        for test_list in [[], ['hi','there','thing'], ['1','2','3']]:
            assert compare_files_with_regex(test_list, test_list) == True
        # a list with all IGNORE lines should match the empty list or any other IGNORE-only lists, but not any other list
        for test_list in [['hi','there','thing'], ['1','2','3']]:
            ignore_list_1 = ['<IGNORE> '+l for l in test_list]
            ignore_list_0 = ['<IGNORE> stuff' for l in test_list]
            ignore_list_2 = ['<IGNORE> '+l+l for l in test_list]
            ignore_list_long = ignore_list_0+ignore_list_1+ignore_list_2
            for ignore_list in [ignore_list_1, ignore_list_0, ignore_list_2]:
                assert compare_files_with_regex(ignore_list, ignore_list) == True
                assert compare_files_with_regex(ignore_list, []) == True
                assert compare_files_with_regex([], ignore_list) == True
                # the two below are test cases for when one iterator stopped and the other didn't:
                #  the first line of the return tuple is an info line I don't feel like checking the details of, 
                #   the second is the extra line from whichever iterator was longer (i.e. first line of test_list).
                assert compare_files_with_regex(ignore_list, test_list)[1] == test_list[0]
                assert compare_files_with_regex(test_list, ignore_list)[1] == test_list[0]
            for list1,list2 in itertools.permutations([[],ignore_list_1,ignore_list_0,ignore_list_2,ignore_list_long],2):
                assert compare_files_with_regex(list1, list2) == True
        # a list should match itself no matter how many IGNORE lines are added. 
        for test_list in [['hi','there','thing'], ['1','2','3']]:
            ignore_list = ['<IGNORE> '+l for l in test_list]
            ignore_list_0 = ['<IGNORE> stuff' for l in test_list]
            ignore_list_2 = ['<IGNORE> '+l+l for l in test_list]
            ignore_list_long = ignore_list_0+ignore_list_1+ignore_list_2
            for list1,list2 in itertools.permutations([[], ignore_list, ignore_list_0, ignore_list_2, ignore_list_long],2):
                assert compare_files_with_regex(test_list, list1+test_list) == True
                assert compare_files_with_regex(test_list+list1, test_list+list2) == True
                assert compare_files_with_regex(test_list+list1+test_list, test_list+list2+test_list) == True
        # a list with all REGEX ".*" lines should match any list of the same length, 
        #  but not of a different length (except IGNORE lines); 
        for test_list in [['hi','there','thing'], ['1','2','3']]:
            any_regex_list = ['<REGEX>.*' for l in test_list]
            assert compare_files_with_regex(any_regex_list, test_list) == True
            assert compare_files_with_regex(test_list, any_regex_list) == True
            ignore_list = ['<IGNORE>']*3
            assert compare_files_with_regex(ignore_list+any_regex_list, test_list) == True
            assert compare_files_with_regex(any_regex_list, ignore_list+test_list) == True
            assert compare_files_with_regex(ignore_list+any_regex_list, test_list+ignore_list+ignore_list) == True
        # comparing two REGEX lines should give an exception; 
        #  comparing two lists with REGEX lines in different positions shouldn't
        for test_list in [['hi','there','thing'], ['1','2','3']]:
            any_regex_list = ['<REGEX>.*' for l in test_list]
            self.assertRaises(ValueError, compare_files_with_regex, any_regex_list, any_regex_list)
            first_regex_list = ['<REGEX>.*'] + test_list[1:]
            last_regex_list = test_list[:-1] + ['<REGEX>.*']
            self.assertRaises(ValueError, compare_files_with_regex, first_regex_list, first_regex_list)
            self.assertRaises(ValueError, compare_files_with_regex, last_regex_list, last_regex_list)
            compare_files_with_regex(first_regex_list, last_regex_list) == True
            compare_files_with_regex(last_regex_list, first_regex_list) == True
        # some specific simple ignore and regex testing
        assert compare_files_with_regex(['1','2','3'],['1','2','3']) == True
        assert compare_files_with_regex(['2','1','3'],['1','2','3']) == ('2','1')
        assert compare_files_with_regex(['12','2','3'],['1','2','3']) == ('12','1')
        assert compare_files_with_regex(['<REGEX>[12]','2','3'],['1','2','3']) == True
        assert compare_files_with_regex(['<IGNORE>12','2','3'],['1','2','3']) == ('2','1')
        assert compare_files_with_regex(['<IGNORE>12','2','3'],['2','3']) == True
        assert compare_files_with_regex(['<IGNORE>12','2','3'],['<IGNORE>1','2','3']) == True
        # some specific testing of more complicated regexes
        assert compare_files_with_regex(['<REGEX>.*'],['some text here']) == True
        assert compare_files_with_regex(['<REGEX>some.*'],['some text here']) == True
        assert compare_files_with_regex(['<REGEX>.*some'],['some text here']) == ('<REGEX>.*some','some text here')
        assert compare_files_with_regex(['<REGEX>.*here'],['some text here']) == True
        assert compare_files_with_regex(['<REGEX>some text here'],['some text here']) == True
        assert compare_files_with_regex(['<REGEX>.*some text here'],['some text here']) == True
        assert compare_files_with_regex(['<REGEX>[sometxhr ]*'],['some text here']) == True
        assert compare_files_with_regex(['<REGEX>[eosmtxhr ]'],['some text here'])==('<REGEX>[eosmtxhr ]','some text here')
        assert compare_files_with_regex(['<REGEX>[omes]*\s[xet]*\s[ehr]*'],['some text here']) == True
        assert compare_files_with_regex(['<REGEX>[oes]*\s[xet]*\s[ehr]*'],['some text here']) != True
        assert compare_files_with_regex(['<REGEX>[omes]\s[xet]*\s[ehr]*'],['some text here']) != True
        assert compare_files_with_regex(['<REGEX>[omes]*[xet]*\s[ehr]*'],['some text here']) != True
        # testing ignore_empty_lines
        for other_list in [['', 'some text'], ['some text', '    \t'], ['', '\t', 'some text', '', '']]:
            assert compare_files_with_regex(['some text'], other_list, ignore_empty_lines=False) != True
            assert compare_files_with_regex(['some text'], other_list, ignore_empty_lines=True) == True
        # testing ignore_whitespace
        for other_list in [['some  text'], ['some\ttext'], ['some \t\t text \t']]:
            assert compare_files_with_regex(['some text'], other_list, ignore_whitespace=False) == ('some text', 
                                                                                                    other_list[0].strip())
            assert compare_files_with_regex(['some text'], other_list, ignore_whitespace=True) == True
            assert compare_files_with_regex(['\t\tsome text'], other_list, ignore_whitespace=False) != True
            assert compare_files_with_regex(['\t\tsome text'], other_list, ignore_whitespace=True) == True
            # removing the whitespace completely doesn't work, though
            assert compare_files_with_regex(['sometext'], other_list, ignore_whitespace=False) == ('sometext',
                                                                                                   other_list[0].strip())
            assert compare_files_with_regex(['sometext'], other_list, ignore_whitespace=True) != True
        # testing ignore_case
        for other_list in [['SOME text'], ['Some Text'], ['sOmE tExT']]:
            assert compare_files_with_regex(['some text'], other_list, ignore_case=False) == ('some text',other_list[0])
            assert compare_files_with_regex(['some text'], other_list, ignore_case=True) == True
            assert compare_files_with_regex(['SOME TEXT'], other_list, ignore_case=False) == ('SOME TEXT',other_list[0])
            assert compare_files_with_regex(['SOME TEXT'], other_list, ignore_case=True) == True
        ### testing on files instead of lists - remember to restart the iterators every time!
        if debug:   print " ************* file tests **************** "
        # non-regex files match themselves
        file1,file2,file3 = 'test_inputs/textcmp_file1.txt','test_inputs/textcmp_file2.txt','test_inputs/textcmp_file3.txt'
        with open(file1,'r') as F1:
            with open(file1,'r') as F1_:
                assert compare_files_with_regex(F1, F1_) == True
        with open(file3,'r') as F3:
            with open(file3,'r') as F3_:
                assert compare_files_with_regex(F3, F3_) == True
        # comparing regex file to itself causes an exception
        with open(file2,'r') as F2:
            with open(file2,'r') as F2_:
                self.assertRaises(ValueError, compare_files_with_regex, F2, F2_)
        # comparing a non-regex file to a matching file with regex/ignore matches
        with open(file1,'r') as F1:
            with open(file2,'r') as F2:
                assert compare_files_with_regex(F1, F2) == True
        # comparing to the third file with extra empty lines, whitespace and case differences - all three True/False 
        #  arguments must be true to give a match, no other combination works.
        for first_file in (file1, file2):
            if debug:   print " ****** %s vs %s ******"%(first_file, file3)
            with open(first_file,'r') as FILE1:
                with open(file3,'r') as F3:
                    assert compare_files_with_regex(FILE1, F3, ignore_empty_lines=False, 
                                                    ignore_whitespace=False, ignore_case=False) != True
            with open(first_file,'r') as FILE1:
                with open(file3,'r') as F3:
                    assert compare_files_with_regex(FILE1, F3, ignore_empty_lines=True, 
                                                    ignore_whitespace=True, ignore_case=True) == True
            for v1,v2,v3 in itertools.chain(itertools.permutations([True,False,False],3), 
                                            itertools.permutations([True,True, False],3)):
                with open(first_file,'r') as FILE1:
                    with open(file3,'r') as F3:
                        assert compare_files_with_regex(FILE1, F3, ignore_empty_lines=v1, 
                                                        ignore_whitespace=v2, ignore_case=v3) != True


if __name__=='__main__':
    """ If module is ran directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
