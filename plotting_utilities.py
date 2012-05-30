#! /usr/bin/env python
"""
Various plotting utilities I wrote (usually for matplotlib) - see docstring for each function for what it does.
 -- Weronika Patena, 2010-2012
"""

# standard library
import itertools
import unittest
# other packages
import matplotlib.pyplot as mplt

# For useful tricks see ~/computers_and_programming/matplotlib_notes_and_tricks.txt file.

################################ EASY FUNCTIONS FOR SPECIFIC PLOT TYPES ###################################

def stacked_bar_plot(list_of_sample_category_lists, sample_names=[], bar_width=0.7, colors='bgrcmy'):
    """ Plot list_of_sample_category_lists as a stacked bar plot (all categories per sample on top of each other). 
    Return list of plot_bar objects, to use in legend (like this: "mplt.legend(plot_bar_list, name_list)"). 
    """
    if not len(set([len(category_list) for category_list in list_of_sample_category_lists])) == 1:
        raise ValueError("All lists in list_of_sample_category_lists must be the same length!")
    if sample_names and not len(sample_names)==len(list_of_sample_category_lists):
        raise ValueError("list_of_sample_category_lists and sample_names must be the same length!")
    N_samples = len(list_of_sample_category_lists)
    N_categories = len(list_of_sample_category_lists[0])
    if not sample_names:
        sample_names = ['' for _ in range(N_samples)]
    positions = range(N_samples)
    category_bars_for_legend = []
    bar_bottoms = [0 for _ in sample_names]
    for category_N, color in zip(range(N_categories), itertools.cycle(colors)):
        category_values = [sample_category_list[category_N] for sample_category_list in list_of_sample_category_lists]
        plot_bars = mplt.bar(positions, category_values, bottom=bar_bottoms, color=color, width=bar_width)
        bar_bottoms = [x+y for (x,y) in zip(bar_bottoms,category_values)]
        category_bars_for_legend.append(plot_bars[0])
    mplt.xticks([p + bar_width/2 for p in positions], sample_names)
    mplt.xlim(-(1-bar_width), mplt.xlim()[1])
    return category_bars_for_legend



################################ COSMETIC MODIFICATIONS TO EXISTING PLOTS ###################################

def color_plot_frame(plot_axes, color='grey', color_frame=True, color_ticks=True, color_ticklabels=True): 
    """ Change the color of the frame/ticks/ticklabels of plot_axes (a matplotlib.axes.AxesSubplot object) to color. """
    # source: http://stackoverflow.com/questions/7778954/elegantly-changing-the-color-of-a-plot-frame-in-matplotlib
    if color_frame:         mplt.setp(plot_axes.spines.values(), color=color)
    if color_ticks:         mplt.setp([plot_axes.get_xticklines(), plot_axes.get_yticklines()], color=color)
    if color_ticklabels:    mplt.setp([plot_axes.get_xticklabels(), plot_axes.get_yticklabels()], color=color)


################## OLD FUNCTIONS, moved from general_utilities.py

def plot_function_by_window_size(data,window_size_list,function,figsize=(),title='',xlabel='',xlim=None,ylim=None,yticks=None,ylabel=''):
    """ Return a plot of function on data, with a subplot for each window_size_list value. """
    if figsize: fig = mplt.figure(figsize=figsize)
    else:       fig = mplt.figure()
    for i in range(len(window_size_list)):
        window_size = window_size_list[i]
        mplt.subplot(len(window_size_list),1,i+1)
        if i==0 and title:                          mplt.title(title)
        if i==len(window_size_list)-1 and xlabel:   mplt.xlabel(xlabel)
        function_values, indices = function(data, window_size, return_indices=True)
        mplt.plot(indices, function_values, 'b,',linestyle='none') 
        if ylim:    mplt.ylim(ylim[0],ylim[1])
        if xlim:    mplt.xlim(xlim[0],xlim[1])
        if yticks:  mplt.yticks(yticks[0],yticks[1] if len(yticks)>1 else yticks[0])
        if ylabel:  mplt.ylabel(ylabel%window_size)
    return fig
    # TODO how would one even test this?


### Convert data into "linlog" scale for plotting (linear for some small range around 0, log after that)
# Note: there used to be an old in-progress implementation of this here, called convert_data_to_linlog, but I removed it on 2012-05-30. 
# The right thing to do is use the "symlog" scale in xscale/yscale, with the appropriate linthreshx/linthreshy (also allows negative values):
# - http://stackoverflow.com/questions/3305865/what-is-the-difference-between-log-and-symlog
# - http://matplotlib.sourceforge.net/examples/pylab_examples/symlog_demo.html


if __name__=='__main__':
    """ If module is ran directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    print "NO TESTS FOR THIS FILE"
    #unittest.main()
