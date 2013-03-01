#! /usr/bin/env python

"""
Various statistical convenience functions I wrote - see docstring for each function for what it does.
Weronika Patena, 2010-2013
"""

# standard library
from __future__ import division 
import sys, os
import unittest
import itertools
from collections import defaultdict
import random
# other packages
import numpy
import scipy.stats
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
R_stats = importr('stats')
# my modules


### HELP FUNCTIONS

def array_1D(x):
    """ Convert to 1D numpy array. """
    return numpy.reshape(numpy.array(x), -1)


### STATISTICAL FUNCTIONS

def chisquare_goodness_of_fit(category_counts, expected_frequencies, dof_subtract=0, return_pvalue_only=True, min_count=50):
    """ Gives p-value for whether a list of category counts is different from the expected frequencies, using the chi-square test.

    Simple wrapper around scipy.stats.chisquare, that calculates the expected counts by normalizing expected_frequencies 
     to have the same total as category_counts.

    Dof_subtract is how many degrees of freedom to SUBTRACT from the default (usually no adjustment needed, 0).

    Raises ValueError if any count is below min_count - the chi-square test shouldn't be used for small numbers.

    If return_pvalue_only is True, returns only pvalue, otherwise (chisquare_statistic, pvalue).
    """
    # convert everything to 1D numpy arrays
    category_counts = array_1D(category_counts)
    expected_frequencies = array_1D(expected_frequencies)
    if numpy.min(category_counts) < min_count:
        raise ValueError("Shouldn't use chisquare_goodness_of_fit with small numbers! Adjust min_count if you have to.")
    # using scipy.sum  and numpy.array in case category_counts/expected_frequencies are multi-dimensional matrices
    norm_factor = sum(category_counts) / sum(expected_frequencies)
    expected_counts_same_total = expected_frequencies * norm_factor
    chisq, pval = scipy.stats.chisquare(category_counts, expected_counts_same_total, ddof=dof_subtract)
    if return_pvalue_only:  return pval
    else:                   return chisq, pval


def chisquare_independence(category_counts_a, category_counts_b, dof_subtract=0, return_pvalue_only=True, min_count=50):
    """ Gives p-value for whether two lists of category counts are different, using the chi-square test.

    Dof_subtract is how many degrees of freedom to SUBTRACT from the default (usually no adjustment needed, 0).

    If return_pvalue_only is True, returns only pvalue, otherwise (chisquare_statistic, pvalue).

    Raises ValueError if any count is below min_count - the chi-square test shouldn't be used for small numbers, 
     use Fisher's exact test instead.

    NOTE ON HOW THIS WORKS: to do a chi-square test of independence, you compare the counts of eiter category to expected counts, 
     where expected counts are CALCULATED FROM BOTH CATEGORIES - you DON'T directly compare one to the other,
      the way you do in a goodness-of-fit test!
        EXAMPLE: if counts_A are 110 and 90, and counts_B are 190 and 10, you don't do a chisquare on [110,90], [190,10], 
         but calculate the expected counts for all four categories from the row/column totals: 
          from A+B (110+190=300 and 90+10=100), so since sum(A) and sum(B) are both 200, the expected is [150,50,150,50]. 
          So compare [110,90,190,10] (A and B together) to [150,50,150,50] using the chi-square goodness-of-fit test, 
          BUT make a degree-of-freedom adjustment: the comparison we're REALLY doing is 2*3 (comparing two 3-length datasets), 
           so the degrees of freedom should be (2-1)*(3-1) = 2.  However, we're transforming it in to a comparison between 
            6-length observed and expected datasets, with (2-1)*(6-1)=5 degrees of freedom by default, 
           so we need to subtract 3 from the degrees-of-freedom for the test (in addition to whatever dof_subtract we already have). 

        EXAMPLE 2: what if the totals in A and B are different?  Say A is [110,190] and B is [90,10]. 
            Then the overall counts are [200,200], so we should be comparing [110,190,90,10] to [150,150,50,50]. 
            Again, with DOF adjustment of 3.

    See https://udel.edu/~mcdonald/statchiind.html for source and description of how it should work, 
     and more examples at http://stattrek.com/chi-square-test/independence.aspx 
     and http://omega.albany.edu:8008/mat108dir/chi2independence/chi2in-m2h.html - it was confusing to me at first.  
    """
    # just convert everything to scipy/numpy arrays
    category_counts_a = array_1D(category_counts_a)
    category_counts_b = array_1D(category_counts_b)

    all_observed = numpy.append(category_counts_a, category_counts_b)

    if min(all_observed) < min_count:
        raise ValueError("Shouldn't use chisquare_goodness_of_fit with small numbers! Adjust min_count if you have to.")

    # calculate the expected frequencies for all categories and both datasets
    both_count_totals = category_counts_a + category_counts_b
    sum_a, sum_b = sum(category_counts_a), sum(category_counts_b)
    full_sum = sum_a + sum_b
    both_count_totals_norm_a = both_count_totals*sum_a/full_sum
    both_count_totals_norm_b = both_count_totals*sum_b/full_sum
    assert sum(both_count_totals_norm_a) == sum_a
    assert sum(both_count_totals_norm_b) == sum_b
    all_expected = numpy.append(both_count_totals_norm_a, both_count_totals_norm_b)
    # Degrees-of-freedom adjustment (see docstring example for details)
    proper_dof = len(category_counts_b) - 1
    transformed_dof = len(all_observed) - 1
    extra_dof_adjustment = proper_dof - transformed_dof
    full_dof_subtract = dof_subtract - extra_dof_adjustment
    # do the chi-square goodness-of-fit test of observed vs expected
    chisq, pval = chisquare_goodness_of_fit(all_observed, all_expected, full_dof_subtract, return_pvalue_only, min_count)
    if return_pvalue_only:  return pval
    else:                   return chisq, pval


######################################## TESTS FOR THIS FILE ########################################

class Testing_everything(unittest.TestCase):
    """ Testing all functions/classes.etc. """

    def test__chisquare_goodness_of_fit(self):
        kwargs = dict(return_pvalue_only=False, min_count=1)
        # test case 1 from https://udel.edu/~mcdonald/statchigof.html
        chisq,pval = chisquare_goodness_of_fit([423,133], [3,1], **kwargs)
        assert round(chisq,2) == 0.35
        assert round(pval,3) == 0.557
        # test case 2 from https://udel.edu/~mcdonald/statchigof.html - note that the numbers here are too low for this test, really
        chisq,pval = chisquare_goodness_of_fit([70, 79, 3, 4], [0.54, 0.4, 0.05, 0.01], **kwargs)
        assert round(chisq,3) == 13.593
        assert round(pval,4) == 0.0035
        # test case 3 from https://udel.edu/~mcdonald/statchigof.html - NOTE that this should have 1 degree of freedom
        chisq,pval = chisquare_goodness_of_fit([14, 21, 25], [0.167, 0.483, 0.350], dof_subtract=1, **kwargs)
        assert round(chisq,2) == 4.54
        assert round(pval,3) == 0.033

    def test__chisquare_independence(self):
        # test case 1 from https://udel.edu/~mcdonald/statchigof.html
        chisq,pval = chisquare_independence([268,199, 42], [807,759,184], dof_subtract=0, return_pvalue_only=False, min_count=1)
        assert round(chisq,2) == 7.26
        assert round(pval,3) == 0.027
        # test case 2 from https://udel.edu/~mcdonald/statchigof.html
        chisq,pval = chisquare_independence([127, 99, 264], [116, 67, 161], dof_subtract=0, return_pvalue_only=False, min_count=1)
        assert round(chisq,2) == 6.26
        assert round(pval,3) == 0.044
        # test case from http://stattrek.com/chi-square-test/independence.aspx
        chisq,pval = chisquare_independence([200, 150, 50], [250, 300, 50], dof_subtract=0, return_pvalue_only=False, min_count=1)
        assert round(chisq,1) == 16.2
        assert round(pval,4) == 0.0003


if __name__=='__main__':
    """ If module is ran directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
