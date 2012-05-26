#! /usr/bin/env python
"""
Various plotting utilities I wrote (usually for matplotlib) - see docstring for each function for what it does.
 -- Weronika Patena, 2010-2012
"""

import matplotlib.pyplot as mplt
import unittest


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


### Convert data into linlog scale (slightly weird, for plotting) - OLD CODE, STILL IN PROGRESS, will pick it back up if it's needed for anything.
# for example, if cutoff=10, what we want is: 1->1, 2->2, 5->5, 10->10, 100->20, 1000->30, 10000->40, ...
#  if cutoff=100, what we want is: 1->1, 2->2, 5->5, 10->10, 50->50, 100->100, 1000->200, 10000->300, 100000->400, ...
#  if cutoff=1, what we want is: 1->1, 10->2, 100->3, 1000->4, ...
# TODO this should probably just return a graph instead of the values!  Since there's no conceivable use for this except for graphing.  And then I could make the graph have extra features, like ticks/labels and a lin/log transition mark.
### TODO this example shows the proper way of implementing a custom scale in matplotlib: 
# http://matplotlib.sourceforge.net/examples/api/custom_scale_example.html
# also see StackOverflow: http://stackoverflow.com/questions/6382612/python-equivalent-for-matlabs-normplot/6382851#6382851
def convert_data_to_linlog(dataset,cutoff=10):
    """ Convert the data (list of numbers) for correct plotting on a scale that's linear up to cutoff and log afterward."""
    from math import log10
    cutoff=float(cutoff)
    dataset_linlog = []
    for x in dataset:
        if x <= cutoff:
            dataset_linlog.append(x)
        else:
            dataset_linlog.append((log10(x/cutoff)+1)*cutoff)
    return dataset_linlog
# TODO now how to label it correctly?  Basically first make the graph, get the resulting tick values, and run convert_data_to_linlog on them to get the tick labels.
# the fuctions is mplt.yticks([0,1,2,5,10,20,30,40],['0','1','2','5','10','100','1000','10000'])
# TODO maybe put a transition mark where the scale changes from lin to log?


if __name__=='__main__':
    """ If module is ran directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    print "NO TESTS FOR THIS FILE"
    #unittest.main()
