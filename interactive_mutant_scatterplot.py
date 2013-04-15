#! /usr/bin/env python

""" Hooking up interactive behavior to mutant readcount scatterplots!
Run connect_to_mutant_scatterplot to add interactive behavior an existing scatterplot, given the mutant dataset. 
Current interactive behavior:
 * click events work on the nearest mutant, or multiple mutants if they're all the same distance (usually same position)
    - left-click highlights the mutant and prints its info; left-click outside the axes un-highlights all mutants.
    - right-click adds the mutant to the selected set, or removes it if it was already in the set, and changes its color.
 * rectangle selection (click-and-drag) works on all mutants inside the rectangle, 
    with the same left-click and right-click functionality.
    CAUTION: SLOW AND BUGGY, DOESN'T ALWAYS WORK
 * keyboard shortcuts:
    > or + increases the opacity of all the dots, and < or - decreases it
    U or u removes all mutants from the selected set
    I or i prints information for all mutants in the selected set
 * the selected mutant dataset is available as interactive_mutant_scatterplot.SELECTED_DATASET
More behavior may be added

 -- Weronika Patena, Jonikas Lab, 2013

Example usage (in python interactive shell, on real mutant data):

    >>> import interactive_mutant_scatterplot
    >>> import mutant_analysis_classes
    >>> datafile1 = '/home/weronika/experiments/mutant_pool_screens/1301_Mia-lipid-screen1_basic-analysis/2_mutants/5prime_no-merging/HL1_5prime_mutants_no-merging.pickle'
    >>> datafile2 = '/home/weronika/experiments/mutant_pool_screens/1301_Mia-lipid-screen1_basic-analysis/2_mutants/5prime_no-merging/LL1_5prime_mutants_no-merging.pickle'
    >>> HL_data = mutant_analysis_classes.read_mutant_file(datafile1)
    >>> LL_data = mutant_analysis_classes.read_mutant_file(datafile2)
    >>> both_data = mutant_analysis_classes.Insertional_mutant_pool_dataset(multi_dataset = True)
    >>> both_data.populate_multi_dataset({'HL': HL_data, 'LL': LL_data})
    >>> mutants_sorted = sorted(both_data, key=lambda m: m.position)
    >>> readcounts_HL = [m.by_dataset['HL'].total_read_count for m in mutants_sorted]
    >>> readcounts_LL = [m.by_dataset['LL'].total_read_count for m in mutants_sorted]

    >>> fig = mplt.figure(figsize=(6,6))
    >>> plot = mplt.scatter(readcounts_HL, readcounts_LL, color='k', edgecolor='none')
    >>> axes = mplt.gca()
    >>> mplt.xlabel('high-lipid readcounts')
    >>> mplt.ylabel('low-lipid readcounts')
    >>> mplt.yscale('symlog')
    >>> mplt.xscale('symlog')
    >>> mplt.ylim(-1, 10000)
    >>> mplt.xlim(-1, 10000)
    >>> reload(interactive_mutant_scatterplot)
    >>> interactive_mutant_scatterplot.connect_to_mutant_scatterplot(fig, axes, plot, 'HL', 'LL', both_data)

    >>> interactive_mutant_scatterplot.disconnect_all(fig)
    >>> mplt.close()
"""
# LATER-TODO generalize this to work with other plot types?

# standard library
from __future__ import division
import sys
import time
import math
# other packages
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import widgets
import numpy as np
# my modules
import mutant_analysis_classes


### globals
# MAYBE-TODO make them not global?
CURRENT_CONNECTIONS = []
CURRENT_HIGHLIGHTS = []


######### Help functions to set up for applying this to mutant scatterplots

def find_mutants_from_click_position(click_x, click_y, x_dataset_name, y_dataset_name, joint_mutants):
    """ Given a click position on a readcount scatterplot, return a list of mutants closest to the click."""
    distance_function = lambda mutant: math.sqrt( (mutant.by_dataset[x_dataset_name].total_read_count - click_x)**2 
                                                 +(mutant.by_dataset[y_dataset_name].total_read_count - click_y)**2 )
    lowest_distance = min(distance_function(mutant) for mutant in joint_mutants)
    closest_mutants = [m for m in joint_mutants if distance_function(m)==lowest_distance]
    # MAYBE-TODO fix this to avoid doing float comparison?
    # MAYBE-TODO fix this distance measure to work more predictably on a log-scale plot!
    return closest_mutants
    # this older version just selected the single closest mutant, but sometimes there are multiple equally close ones - pick all!
    #closest_mutant = min(joint_mutants, key = distance_function)

def find_mutants_in_rectangle(x_low, y_low, x_high, y_high, x_dataset_name, y_dataset_name, joint_mutants):
    """ Given the x and y edges of a rectangle on a readcount scatterplot, return a list of mutants inside the rectangle."""
    return [m for m in joint_mutants if x_low <= m.by_dataset[x_dataset_name].total_read_count <= x_high
                                        and y_low <= m.by_dataset[y_dataset_name].total_read_count <= y_high]


def percent_perfect(mutant, dataset):
    """ Return string giving the %perfect reads, or '-' if 0 total reads."""
    total_readcount = mutant.by_dataset[dataset].total_read_count
    if total_readcount:     return '%.0f%%'%(mutant.by_dataset[dataset].perfect_read_count/total_readcount*100)
    else:                   return 'N/A'

def mutant_info_string(mutant):
    m = mutant
    datasets = sorted(m.by_dataset.keys())
    return "%s, readcounts %s (perfect %s), gene %s %s %s"%(m.position, 
                                                        ', '.join('%s %s'%(d, m.by_dataset[d].total_read_count) for d in datasets), 
                                                        ', '.join(percent_perfect(m, d) for d in datasets), 
                                                        m.gene, m.gene_feature, m.orientation)
    # TODO add annotation info, when there is some!
    # MAYBE-TODO write a few different options for the info string, and a function that'll switch between different printing versions based on keypress/something, or even just set the option when doing the interactive connection


######### Interactive behavior functions

def change_selection(method, mutants, scatterplot, joint_mutants, SELECTED_DATASET):
    allowed_methods = 'flip add remove'.split()
    if method not in allowed_methods:
        raise ValueError('change_selection method must be one of %s!'%allowed_methods)
    N_added, N_removed = 0, 0 
    for mutant in mutants:
        if method=='add' or (method=='flip' and mutant.position not in SELECTED_DATASET):
            SELECTED_DATASET.add_mutant(mutant)
            N_added += 1
        elif mutant.position in SELECTED_DATASET:
            SELECTED_DATASET.remove_mutant(mutant)
            N_removed += 1
    # set the colors to match the current SELECTED_DATASET
    # TODO this is sort of slow - make it faster!
    scatterplot.set_facecolors([('dodgerblue' if m.position in SELECTED_DATASET else 'k') 
                                for m in sorted(joint_mutants, key = lambda m: m.position)])
    # MAYBE-TODO pick a better color?  Or make it customizable.
    if len(mutants):
        print "Clicked %s mutants - added %s to selection, removed %s."%(
                        'ALL' if len(mutants)==len(joint_mutants) else len(mutants), N_added, N_removed)


def highlight_new_mutants(axes, x_dataset_name, y_dataset_name, mutants):
        # remove current highlights; if clicked outside the plot, only this is done.
        global CURRENT_HIGHLIGHTS
        for circle in CURRENT_HIGHLIGHTS:  circle.remove()
        # print info - one line for single mutant, otherwise give the total number and the top 5 mutants
        if len(mutants) == 1:
            print mutant_info_string(mutants[0])
        elif len(mutants):
            # LATER-TODO make the number customizable
            N_to_show = 5
            print "Picked %s mutants%s:"%(len(mutants), ' (showing %s)'%N_to_show if len(mutants) > N_to_show else '')
            for mutant in mutants[:N_to_show]:
                print '  - %s'%(mutant_info_string(mutant))
            if len(mutants) > N_to_show:
                print '  ...'
        ### make it highlight the dot to make sure it's the right ones!
        # MAYBE-TODO instead of putting highlights over existing mutants, actually change the dot colors in the plot?
        # draw new highlights; save them so they can be removed later
        if len(mutants):
            circles = axes.plot([m.by_dataset[x_dataset_name].total_read_count for m in mutants], 
                                [m.by_dataset[y_dataset_name].total_read_count for m in mutants],
                                'o', markersize=6, markeredgewidth=2, markeredgecolor='red', markerfacecolor='none')
            CURRENT_HIGHLIGHTS = circles
        else:
            CURRENT_HIGHLIGHTS = []


def on_click(event, fig, axes, scatterplot, x_dataset_name, y_dataset_name, joint_mutants, SELECTED_DATASET):
    x, y = event.xdata, event.ydata
    # this is just so I can look at an example event's attributes/methods interactively
    global some_mouseclick_event
    some_mouseclick_event = event
    # All attributes of a click event:  
    #   button (which mouse button), xdata and ydata (position in data coordinates), x and y (position in ?? coordinates),
    #   canvas, inaxes (which axes the click was in), name ('button_press_event'), 
    #   key (? - empty when I looked), lastevent (? - empty when I looked me), step (? - 0 when I looked), guiEvent (?)
    #   also has a _update_enter_leave method.

    ### find the closest mutants if clicked inside the plot; 
    # if clicked outside, assume no mutants (will only remove previous highlights)
    if event.inaxes is not None:
        closest_mutants = find_mutants_from_click_position(x, y, x_dataset_name, y_dataset_name, joint_mutants)
    else:
        closest_mutants = []
    ### left-click - print mutant info
    if event.button == 1:
        highlight_new_mutants(axes, x_dataset_name, y_dataset_name, closest_mutants)
    ### right-click - add/remove mutants to SELECTED_DATASET, print how many were added/removed
    elif event.button == 3:
        change_selection('flip', closest_mutants, scatterplot, joint_mutants, SELECTED_DATASET)
        # MAYBE-TODO is that the behavior I want - for a right-click to both add and remove?  Otherwise I can use method add/remove.
    fig.canvas.draw()


def on_rectangle_select(start_event, end_event, fig, axes, scatterplot, x_dataset_name, y_dataset_name, 
                        joint_mutants, SELECTED_DATASET):
    """ Rectangle selection behavior - see module docstring for detail. """
    x_low, x_high = sorted((start_event.xdata, end_event.xdata))
    y_low, y_high = sorted((start_event.ydata, end_event.ydata))
    # TODO this is extremely slow and weird and often doesn't work!!  Not sure why...  Interference between click and RectangleSelect?  Test it, read the docs!  I could make a keybinding switch between the two...  

    ### find the mutants inside the rectangle (I think the rectangle won't work outside the plot):
    mutants_in_rectangle = find_mutants_in_rectangle(x_low, y_low, x_high, y_high, x_dataset_name, y_dataset_name, joint_mutants)

    ### left-click - print mutant info
    if start_event.button == 1:
        highlight_new_mutants(axes, x_dataset_name, y_dataset_name, mutants_in_rectangle)
    ### right-click - add/remove mutants to SELECTED_DATASET, print how many were added/removed
    elif start_event.button == 3:
        change_selection('flip', mutants_in_rectangle, scatterplot, joint_mutants, SELECTED_DATASET)
    fig.canvas.draw()
    # MAYBE-TODO can't extend the rectangle past the axes border - that's sort of annoying... Is there an easy fix?

# MAYBE-TODO implement Lasso selection?
    # Actually Lasso isn't what I want, because it has a defined start point...  I want LassoSelector, which exists in newer matplotlib (AND has better documentation) - http://matplotlib.org/api/widgets_api.html#matplotlib.widgets.LassoSelector
    # To switch between Rectangle and Lasso, do something like this:
    #    if not on_key_press.RectangleSelect.active:
    #        print ' RectangleSelector activated.'
    #        on_key_press.RectangleSelect.set_active(True)


def on_key_press(event, fig, scatterplot, joint_mutants, SELECTED_DATASET):
    """ Keyboard key-press triggers - see module docstring for detail. """
    # Note: there are also some preexisting matplotlib keybindings - I don't know how to turn them off, not sure if I want to.  Sometimes these will interact with them.
    ### this is just to ignore it when the user presses something non-alphanumeric (otherwise it prints warnings)
    if event.key is None:
        return
    ### Increasing/decreasing transparency
    # MAYBE-TODO customize the step?
    step = 0.1
    curr_alpha = scatterplot.get_alpha()
    if curr_alpha is None:  curr_alpha = 1
    if event.key in '+>':
        scatterplot.set_alpha(min(curr_alpha + step, 1))
    if event.key in '-<':
        scatterplot.set_alpha(max(curr_alpha - step, 0))
    ### Printing or unselecting selected circles
    if event.key in 'Uu':
        change_selection('remove', joint_mutants, scatterplot, joint_mutants, SELECTED_DATASET)
    if event.key in 'Ii':
        print "%s selected mutants:"%len([m for m in SELECTED_DATASET])
        # TODO have to use "len([m for m in SELECTED_DATASET])" instead of just "len(SELECTED_DATASET)" because the latter gives an error on multi-dataset mutants - fix that?
        for mutant in SELECTED_DATASET:
            print '  - %s'%(mutant_info_string(mutant))
    ### other possible functionality
    # manipulating zorder (randomize, or put all selected dots on top) - can't do that for scatterplot data, I think!
    ### this is just so I can look at an example event's attributes/methods interactively
    global some_keypress_event
    some_keypress_event = event
    ### draw the new display
    fig.canvas.draw()


def connect_to_mutant_scatterplot(fig, axes, scatterplot, x_dataset_name, y_dataset_name, joint_mutants):
    global SELECTED_DATASET
    # MAYBE-TODO make SELECTED_DATASET not global?  Right now it's global because I don't know how else to get it from the outside.

    # set up rectangle selection - nonzero minspanx keeps a simple click from being interpreted as a tiny rectangle
    # LATER-TODO why is RectangleSelect an attribute of on_key_press, anyway?
    # TODO minspanx should probably be customizable...
    on_key_press.RectangleSelect = widgets.RectangleSelector(axes, lambda start,end: on_rectangle_select(start, end, fig, axes, scatterplot, x_dataset_name, y_dataset_name, joint_mutants, SELECTED_DATASET), minspanx=0.01)
    # connect event handlers to figure
    global CURRENT_CONNECTIONS
    SELECTED_DATASET = mutant_analysis_classes.Insertional_mutant_pool_dataset()
    cID = fig.canvas.mpl_connect('key_press_event', lambda event: on_key_press(event, fig, scatterplot, 
                                                                               joint_mutants, SELECTED_DATASET))
    CURRENT_CONNECTIONS.append(cID) 
    cID = fig.canvas.mpl_connect('button_press_event', lambda event: on_click(event, fig, axes, scatterplot, 
                                                              x_dataset_name, y_dataset_name, joint_mutants, SELECTED_DATASET))
    CURRENT_CONNECTIONS.append(cID) 
    # MAYBE-TODO another option is to do something with pick_events...
    #fig.canvas.mpl_connect('pick_event', lambda event: on_pick_normal(event, fig, circle_objects, circle_data, IDs_picked))


def disconnect_all(fig):
    """ Disconnect all connections; remove current mutant highlights. """
    # remove all highlights - can pass None for axes and dataset names, since they're not needed, and [] for new highlighted mutants
    highlight_new_mutants(None, None, None, mutants=[])
    # disconnect all current connection IDs
    global CURRENT_CONNECTIONS
    for connection_ID in CURRENT_CONNECTIONS:
        fig.canvas.mpl_disconnect(connection_ID)
    CURRENT_CONNECTIONS = []
    # That doesn't seem to WORK...  LATER-TODO fix?
    # Alternative version, temporary, to just disconnect all connection IDs (seems to work):
    for x in range(10000): fig.canvas.mpl_disconnect(x)


# TODO write a "if __name__=='__main__' section for this to work as a shell-script that takes the two mutant files and shows the plot with all the interactive features, printing data to the terminal and saving the selected mutants to a new file if desired.  It'd be good to be able to email people a script they can run and do this interactive thing themselves!

# Does this interactive stuff even work on a Mac?  TODO check!  Maybe try it in ipython notebook (have to use the other pylab setting, not inline) - that's what the guy at PyData used, I think.  That would be harder to put in a standalone script, though...

# It would still require people to have python installed, but they mostly do.  Unless I did py2exe or something...  MAYBE-TODO look into that...
