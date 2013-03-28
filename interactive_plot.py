#! /usr/bin/env python

"""
More playing around with making interactive plots - getting close to something useful!

Normally, clicking on an object adds it to the selected set, or removes it if it was already selected (selected objects are red, others are blue); you can click-and-drag to make a rectangle and toggle the selected status of all points in the rectangle.

Right-clicking prints the info on the given object without selecting/deselecting it (and makes its edge thicker and green). 
Pressing the C key resets all the clicked circles to normal appearance (again, this has nothing to do with selecting/deselecting).

Pressing +/- or </> on the keyboard changes the transparency of all objects. 
Pressing * randomizes the z-order (so that different ones are on top, if they were overlapping); pressing ^ makes all the selected objects be on top of the non-selected ones.

 -- Weronika Patena, 2013

USAGE: data_selection_methods.py [options]
"""

# standard library
from __future__ import division
import sys
import time
# other packages
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import widgets
import numpy as np
# my modules

def define_option_parser():
    """ Populates and returns an optparse option parser object, with __doc__ as the usage string."""
    from optparse import OptionParser
    parser = OptionParser(__doc__)
    parser.add_option('-n', '--N_circles', type='int', default='50', metavar='N', 
                      help="(default %default).")
    parser.add_option('-t', '--starting_transparency', type='float', default='0.5', metavar='X', 
                      help="(default %default).")
    return parser

### making circles and printing them

def make_random_circles(total=50, starting_opacity=1):
    fig, ax = plt.subplots()
    circle_objects = {}
    circle_data = {}
    for N in range(total):
        x, y, r = np.random.random(3)
        circ = plt.Circle((x,y), r/10, facecolor='blue', alpha=starting_opacity, picker=True)
        ax.add_patch(circ)
        circle_objects[id(circ)] = circ
        circle_data[id(circ)] = N
    fig.canvas.draw()
    return circle_objects, circle_data, fig, ax

def circle_info_string(circle_artist):
    return "radius %.3f, x/y position %.3f/%.3f"%( circle_artist.get_radius(), circle_artist.center[0], circle_artist.center[1])

def print_circle_set_info(circles_by_ID, circle_data, ID_set, set_name='in set'):
    picked_circle_info = ["circle number %s: %s"%(circle_data[ID], circle_info_string(circles_by_ID[ID])) for ID in IDs_picked]
    if picked_circle_info:  picked_circle_string = '\n - ' + '\n - '.join(picked_circle_info)
    else:                   picked_circle_string = '-'
    return "%s circles %s: %s"%(len(IDs_picked), set_name, picked_circle_string)

def print_picked_circle_data(circles_by_ID, circle_data, IDs_picked):
    print print_circle_set_info(circles_by_ID, circle_data, IDs_picked, 'picked')
    print "%s circles not picked."%(len(set(circles_by_ID) - IDs_picked))

### interactive behavior options

def select_artist(artist, IDs_picked):
    artist.set_facecolor('red')
    IDs_picked.add(id(artist))

def unselect_artist(artist, IDs_picked):
    artist.set_facecolor('blue')
    IDs_picked.discard(id(artist))

def toggle_select_artist(artist, IDs_picked):
    if id(artist) in IDs_picked:    unselect_artist(artist, IDs_picked)
    else:                           select_artist(artist, IDs_picked)

def highlight_artist(artist):
    # MAYBE-TODO instead of changing properties, draw a separate circle/something around the object to highlight it?
    artist.set_edgecolor('green')
    artist.set_linewidth(10)

def unhighlight_artist(artist):
    artist.set_edgecolor('black')
    artist.set_linewidth(1)

def on_key_press(event, fig, objects_by_ID, object_data, IDs_picked):
    """ Keyboard key-press triggers - see module docstring for detail. """
    # Note: there are also some preexisting matplotlib keybindings - I don't know how to turn them off, not sure if I want to.  Sometimes these will interact with them.
    ### this is just to ignore it when the user presses something non-alphanumeric (otherwise it prints warnings)
    if event.key is None:
        return
    ### manipulating transparency
    curr_alpha = objects_by_ID.values()[0].get_alpha()
    if event.key in '+>':
        for artist in objects_by_ID.values():  artist.set_alpha(min(curr_alpha + 0.2, 1))
    if event.key in '-<':
        for artist in objects_by_ID.values():  artist.set_alpha(max(curr_alpha - 0.2, 0))
    ### manipulating zorder
    if event.key in '*':
        N_objects = len(objects_by_ID)
        for artist in objects_by_ID.values():   artist.set_zorder(np.random.randint(100)+1)
    if event.key in '^':
        N_objects = len(objects_by_ID)
        for artist in objects_by_ID.values():   artist.set_zorder(1)
        for ID in IDs_picked:                   objects_by_ID[ID].set_zorder(2)
    ### manipulating clicked/unclicked or selected/unselected circles
    if event.key in 'Cc':
        print "Unclicked all circles"
        for artist in objects_by_ID.values():   unhighlight_artist(artist)
    if event.key in 'Uu':
        print "Unselected all circles"
        for artist in objects_by_ID.values():   unselect_artist(artist, IDs_picked)
    if event.key in 'Pp':
        print_picked_circle_data(objects_by_ID, object_data, IDs_picked)
    ### this is just so I can look at an example event's attributes/methods interactively
    global some_keypress_event
    some_keypress_event = event
    ### draw the new display
    fig.canvas.draw()

def on_pick_normal(event, fig, objects_by_ID, object_data, IDs_picked):
    """ Mouse-button-press behavior - see module docstring for detail. """
    artist = event.artist
    ID = id(artist)
    try:                
        if event.mouseevent.button==1:
            toggle_select_artist(artist, IDs_picked)
        else:                               # button 2 is middle-click
            # if this was the previously clicked object, just unclick it
            if artist.get_linewidth() == 10:
                print "Unclicked circle number %s"%(object_data[id(artist)])
                unhighlight_artist(artist)
            # otherwise mark it with a wide green border and print info
            else:
                print "Clicked on circle number %s: %s"%(object_data[id(artist)], circle_info_string(artist))
                highlight_artist(artist)
                # MAYBE-TODO also un-click all the other dots?  Actually keeping them clicked isn't bad...
    except KeyError:
        print "Clicked on unidentified object (id %s) - not sure what to do!"%(id(artist))
    fig.canvas.draw()
    # MAYBE-TODO right-click on a non-object location (background) should unclick all the green-border dots, same as pressing C does!  Not sure how to do that.
    # this is just so I can look at an example event's attributes/methods interactively
    global some_mouseclick_event
    some_mouseclick_event = event

def on_rectangle_select(start_event, end_event, objects_by_ID, IDs_picked):
    """ Rectangle selection behavior - see module docstring for detail. """
    x_low, x_high = sorted((start_event.xdata, end_event.xdata))
    y_low, y_high = sorted((start_event.ydata, end_event.ydata))
    N_changed = 0
    for N,artist in objects_by_ID.items():
        x,y = artist.center
        if x_low <= x <= x_high and y_low <= y <= y_high:
            toggle_select_artist(artist, IDs_picked)
            N_changed += 1
    print "Toggled selection status for %s circles."%N_changed
    fig.canvas.draw()
    # TODO we'd like some rectangle de-select method!  Either make the rectangle toggle selection, or make left-click rectangle deselect, or make two separate behaviors for rectangle-select and rectangle-deselect with separate keybindings...
    # MAYBE-TODO can't extend the rectangle past the axes border - that's sort of annoying... Is there an easy fix?
    # LATER-TODO do I have to use the raw positions here, or can I connect the RectangleSelector to a picker?  Although I guess for scatterplot data I'll have to work with positions ANYWAY, so this is probably fine. 

# LATER-TODO implement the Lasso selection!
# Actually Lasso isn't what I want, because it has a defined start point...  I want LassoSelector, which exists in newer matplotlib (AND has better documentation) - http://matplotlib.org/api/widgets_api.html#matplotlib.widgets.LassoSelector
# To switch between Rectangle and Lasso, do something like this:
#    if not on_key_press.RectangleSelect.active:
#        print ' RectangleSelector activated.'
#        on_key_press.RectangleSelect.set_active(True)

if __name__=='__main__':
    parser = define_option_parser()
    options,args = parser.parse_args()
    if args:    sys.exit("\nError: this function takes no positional arguments, only options!")

    # draw the circles; starting with an empty set of picked IDs
    circle_objects, circle_data, fig, ax = make_random_circles(options.N_circles, options.starting_transparency)
    IDs_picked = set()

    ### set desired interactive behavior
    # set up rectangle selection - nonzero minspanx keeps a simple click from being interpreted as a tiny rectangle
    on_key_press.RectangleSelect = widgets.RectangleSelector(ax, lambda c,r: on_rectangle_select(c,r, circle_objects, IDs_picked), 
                                                             button=[1], minspanx=0.01)
    # connect event handlers to figure
    fig.canvas.mpl_connect('pick_event', lambda event: on_pick_normal(event, fig, circle_objects, circle_data, IDs_picked))
    fig.canvas.mpl_connect('key_press_event', lambda event: on_key_press(event, fig, circle_objects, circle_data, IDs_picked))

    # Wait for user to finish playing around, then print the picked dot info
    raw_input("Press Enter when done with the plot (workaround for fvwm closing plot windows prematurely)...\n")
    print_picked_circle_data(circle_objects, circle_data, IDs_picked)
