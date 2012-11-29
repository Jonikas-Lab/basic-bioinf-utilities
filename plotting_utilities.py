#! /usr/bin/env python
"""
Various plotting utilities I wrote (usually for matplotlib) - see docstring for each function for what it does.
 -- Weronika Patena, 2010-2012
"""

# standard library
import itertools
import unittest
import os
# other packages
import numpy
import matplotlib
import matplotlib.pyplot as mplt

# For useful tricks see ~/computers_and_programming/matplotlib_notes_and_tricks.txt file.


def savefig(figname, padding=0.2, dpi=300):
    """ Save current figure as figname, with bbox_inches='tight' and given padding and dpi. """
    mplt.savefig(figname, bbox_inches='tight', pad_inches=padding, dpi=dpi)


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

# NOTE: mplt.gca() is get_current_axis, for when I was to act on it directly with ax.* methods and I don't have an ax arg given;
#       mplt.sca(ax) is set_current_axis, for when I have an ax input and I want to act on it indirectly with mplt.* methods.
# The mplt.* methods are interactive, but the ax.* methods etc need an mplt.draw() to show the results.

# MAYBE-TODO could make the "if ax is None:  ax = mplt.gca(); else:  mplt.sca(ax)" line into a decorator, since it shows up everywhere

def remove_legend(ax=None):
    """ Remove legend for ax or the current axes (detected with gca()). """
    # from Scipy matplotlib cookbook - http://www.scipy.org/Cookbook/Matplotlib/Legend
    if ax is None:  ax = mplt.gca()
    else:           mplt.sca(ax)
    ax.legend_ = None
    # alternative version here - http://stackoverflow.com/questions/5735208/remove-the-legend-on-a-matplotlib-figure
    # ax.legend().set_visible(False)
    # yes, the draw is actually needed in this case!
    mplt.draw()


def remove_xticklabels(ax=None):
    """ Remove x tick labels (leaving the ticks unchanged); acts on ax, or current axes if ax is None. """
    if ax is None:  ax = mplt.gca()
    else:           mplt.sca(ax)
    xlim = mplt.xlim()
    xticks = mplt.xticks()[0]
    mplt.xticks(xticks, [''] * len(xticks))
    mplt.xlim(xlim)

def remove_yticklabels(ax=None):
    """ Remove y tick labels (leaving the ticks unchanged); acts on ax, or current axes if ax is None. """
    if ax is None:  ax = mplt.gca()
    else:           mplt.sca(ax)
    lim = mplt.ylim()
    yticks = mplt.yticks()[0]
    mplt.yticks(yticks, [''] * len(yticks))
    mplt.ylim(lim)


def set_axes_limits(x_min=None, x_max=None, y_min=None, y_max=None, ax=None):
    """ Set whichever limits aren't None, keeping the others the same. """
    if ax is None:  ax = mplt.gca()
    else:           mplt.sca(ax)
    # MAYBE-TODO could add an option to only increase (take max of current and new value), or only decrease, etc...
    if ax is None:  ax = mplt.gca()
    if x_min is None:   x_min = mplt.xlim()[0]
    if x_max is None:   x_max = mplt.xlim()[1]
    mplt.xlim((x_min, x_max))
    if y_min is None:   y_min = mplt.ylim()[0]
    if y_max is None:   y_max = mplt.ylim()[1]
    mplt.ylim((y_min, y_max))


def color_plot_frame(ax=None, color='grey', color_frame=True, color_ticks=True, color_ticklabels=True): 
    """ Change the color of the frame/ticks/ticklabels of ax (a matplotlib.axes.AxesSubplot object) to color. """
    # source: http://stackoverflow.com/questions/7778954/elegantly-changing-the-color-of-a-plot-frame-in-matplotlib
    if ax is None:  ax = mplt.gca()
    else:           mplt.sca(ax)
    if color_frame:         mplt.setp(ax.spines.values(), color=color)
    if color_ticks:         mplt.setp([ax.get_xticklines(), ax.get_yticklines()], color=color)
    if color_ticklabels:    mplt.setp([ax.get_xticklabels(), ax.get_yticklabels()], color=color)

# more on modifying frames here: http://matplotlib.org/examples/pylab_examples/spine_placement_demo.html

def remove_frame(ax=None):
    """ Remove plot frame, including x/y ticks and ticklabels. """
    if ax is None:  ax = mplt.gca()
    else:           mplt.sca(ax)
    ax.set_frame_on(False)
    mplt.xticks([], [])
    mplt.yticks([], [])
    mplt.draw()

def remove_half_frame(ax=None):
    """ Remove the top and right sides of the plot frame, including x/y ticks. """
    # source: http://stackoverflow.com/questions/925024/how-can-i-remove-the-top-and-right-axis-in-matplotlib#925141
    if ax is None:  ax = mplt.gca()
    else:           mplt.sca(ax)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    mplt.draw()


################################ OTHER ADDITIONS ###################################

######### colormaps

def all_colormaps_image(save_as='~/computers_and_programming/colormaps_all_matplotlib.png', close=False):
    """ Make an image with all the colormaps shown; optionally save to file. """
    # code from http://matplotlib.org/examples/pylab_examples/show_colormaps.html
    a = numpy.linspace(0, 1, 256).reshape(1,-1)
    a = numpy.vstack((a,a))
    # Get a list of the colormaps in matplotlib.  Ignore the ones that end with '_r' because these are simply 
    #  reversed versions of ones that don't end with '_r'
    maps = sorted(m for m in mplt.cm.datad if not m.endswith("_r"))
    nmaps = len(maps) + 1
    fig = mplt.figure(figsize=(5,10))
    fig.subplots_adjust(top=0.99, bottom=0.01, left=0.2, right=0.99)
    for i,m in enumerate(maps):
        ax = mplt.subplot(nmaps, 1, i+1)
        mplt.axis("off")
        mplt.imshow(a, aspect='auto', cmap=mplt.get_cmap(m), origin='lower')
        pos = list(ax.get_position().bounds)
        fig.text(pos[0] - 0.01, pos[1], m, fontsize=10, horizontalalignment='right')
    mplt.draw()
    if save_as:     
        savefig(os.path.expanduser(save_as))
    if close:
        mplt.close()

def add_colormap_to_matplotlib(cmap_data, cmap_name, add_reverse_too=True):
    """ Add a colormap to the matplotlib colormap list. """
    # Based on code in https://github.com/dmcdougall/matplotlib/blob/85a26dca55c788c97bbf312b0732db222df55530/lib/matplotlib/cm.py
    # Basically the point is to get matplotlib.cm.get_cmap(cmap_name) to return a cmap instance.
    cm = matplotlib.cm
    cm.datad[cmap_name] = cmap_data
    cm.cmap_d[cmap_name] = cm._generate_cmap(cmap_name, cm.LUTSIZE)
    if add_reverse_too:
        reverse_name = cmap_name + '_r'
        matplotlib.cm.datad[reverse_name] = matplotlib.cm._reverse_cmap_spec(cmap_data)
        cm.cmap_d[reverse_name] = cm._generate_cmap(reverse_name, cm.LUTSIZE)

# CMRmap, from https://github.com/matplotlib/matplotlib/pull/496/files
#  (I think it should be in matplotlib already, but it's not in mine!)
_CMRmap_data = {'red'   : ( (0.000, 0.00, 0.00),
                            (0.125, 0.15, 0.15),
                            (0.250, 0.30, 0.30),
                            (0.375, 0.60, 0.60),
                            (0.500, 1.00, 1.00),
                            (0.625, 0.90, 0.90),
                            (0.750, 0.90, 0.90),
                            (0.875, 0.90, 0.90),
                            (1.000, 1.00, 1.00) ),
                'green' : ( (0.000, 0.00, 0.00),
                            (0.125, 0.15, 0.15),
                            (0.250, 0.15, 0.15),
                            (0.375, 0.20, 0.20),
                            (0.500, 0.25, 0.25),
                            (0.625, 0.50, 0.50),
                            (0.750, 0.75, 0.75),
                            (0.875, 0.90, 0.90),
                            (1.000, 1.00, 1.00) ),
                'blue':   ( (0.000, 0.00, 0.00),
                            (0.125, 0.50, 0.50),
                            (0.250, 0.75, 0.75),
                            (0.375, 0.50, 0.50),
                            (0.500, 0.15, 0.15),
                            (0.625, 0.00, 0.00),
                            (0.750, 0.10, 0.10),
                            (0.875, 0.50, 0.50),
                            (1.000, 1.00, 1.00) )}
add_colormap_to_matplotlib(_CMRmap_data, 'CMRmap')

# MAYBE-TODO add other nice colormaps, from elsewhere or mine:
#  - https://www.mathworks.com/matlabcentral/fileexchange/26026


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
