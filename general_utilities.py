#! /usr/bin/env python

"""
Various general programming/math utilities I wrote - see docstring for each function for what it does.
Weronika Patena, 2010-2011
"""

from __future__ import division 
import sys, os
import subprocess
from collections import defaultdict
import unittest


######################################## STRUCTURES / CLASSES / ETC ########################################

### list/dictionary/etc utilities

def compare_lists_unordered(list1,list2):
    """ Given two lists, return True if they contain the same elements, False otherwise.
    The same could be done with bool(set(list1)==set(list2)), except when list elements are unhashable. """
    if len(list1)!=len(list2):  return False
    for x in list1:
        if x not in list2:      return False
    return True

def reduce_dicts_to_overlaps(dict_list, exception_on_different_values=False):
    """ Given a list of dictionaries, return a similar new list of dicts with only the keys shared by all the dictionaries.
    Caution: returns a LIST OF DICTS, NOT a new dictionary that contains only the overlap! Easy to get confused."""
    # get the set of keys present in all dictionaries (set intersection)
    if len(dict_list)==0:   return []
    # generate list of keys that are in ALL dictionaries (starting from the list of the first one and reducing)
    overlap_keys = set(dict_list[0].keys())
    for d in dict_list:     overlap_keys &= set(d.keys())
    # make a new dictionary list
    new_dict_list = []
    for d in dict_list:     
        new_reduced_dict = dict([(k,v) for (k,v) in d.iteritems() if k in overlap_keys])
        new_dict_list.append(new_reduced_dict)
    return new_dict_list

def count_list_values(input_list):
    """ Given a list, return a value:number_of_times_value_occurred_in_list dictionary. """
    value_counts = defaultdict(lambda: 0)
    for value in input_list:
        value_counts[value] += 1
    return dict(value_counts)
    # Note: this is basically a very simplified version of collections.Counter - which was only added in python 2.7, and I'm still on 2.6, so I can't use it for now.  Should switch to that once I'm on 2.7 though.. 

def add_dicts_of_ints(dict1, dict2, recursive=False, 
                      allow_nonconflicting_values=False, remove_nonallowed_values=False):
    """ Add int-valued dicts together by adding the values: given {a:1, b:1} and {a:2, c:2}, return {a:3, b:1, c:2}. 

    If recursive is True, if both dicts have another dict as a value for a given key, recursively apply 
     add_dicts_of_ints to those values.
    Normally only integer values are allowed (and possibly dictionary values, if recursive is True), and other values
     cause a ValueError. If remove_nonallowed_values is True, other values are silently ignored instead; 
     if allow_nonconflictiong_values is True, other values are kept if the key is only present in one dict
      or if the value is the same in both dicts.
    This uses explicit type checks for ints, not duck typing, sorry - don't want to add strings/lists/etc! 
    """
    # MAYBE-TODO add an argument that'd be a list of other types that should be treated like ints (i.e. added)?
    keys_1_only = set(dict1.keys()) - set(dict2.keys())
    keys_2_only = set(dict2.keys()) - set(dict1.keys())
    keys_both = set(dict2.keys()) & set(dict1.keys())
    new_dict = {}
    ### for keys in both dicts:
    for key in keys_both:
        # if both values are ints, use their sum as new value
        if type(dict1[key]) == type(dict2[key]) == int:
            new_dict[key] = dict1[key] + dict2[key]
        # if recursive is True and both values are dicts, run add_dicts_of_ints on them with same options for new value
        elif recursive and type(dict1[key]) == type(dict2[key]) == dict:
            new_dict[key] = add_dicts_of_ints(dict1[key], dict2[key], recursive=recursive, 
                                              allow_nonconflicting_values=allow_nonconflicting_values, 
                                              remove_nonallowed_values=remove_nonallowed_values)
        # if allow_nonconflicting_values is True and both values are the same, use that for new value
        elif allow_nonconflicting_values and dict1[key] == dict2[key]:
            new_dict[key] = dict1[key]
        # if the values aren't any of the above cases and remove_nonallowed_values is True:
        elif remove_nonallowed_values:
            # if one is an allowed type and the other isn't, keep the allowed one and ignore the other
            #  (really the allowed one should just be moved to the keys_*_only set, in case further checking is needed)
            if recursive:       allowed_types = [int,dict]
            else:               allowed_types = [int]
            if type(dict1[key]) in allowed_types and type(dict2[key]) not in allowed_types:
                keys_1_only.add(key)
            elif type(dict2[key]) in allowed_types and type(dict1[key]) not in allowed_types:
                keys_2_only.add(key)
            # otherwise (if neither are an allowed type, or both are but they're different, like int/dict) just ignore both
            else:
                continue
        # otherwise raise exception
        else:
            raise ValueError("Encountered non-allowed value in one of the dictionaries! See docstring for options.")
    ### for keys that are only in one of the dicts:
    for key_set, curr_dict in [(keys_1_only, dict1), (keys_2_only, dict2)]:
        for key in key_set:     
            # copy ints to new_dict; if allow_nonconflicting_values is True, copy any values without checking
            if allow_nonconflicting_values or type(curr_dict[key])==int:
                new_dict[key] = curr_dict[key]
            # if recursive is True and the value is a dict, apply add_dicts_of_ints to it and {} to get new value
            elif recursive and type(curr_dict[key])==dict:
                new_dict[key] = add_dicts_of_ints(curr_dict[key], {}, recursive=recursive, 
                                                  allow_nonconflicting_values=allow_nonconflicting_values, 
                                                  remove_nonallowed_values=remove_nonallowed_values)
            # if the value isn't any of the above and remove_nonallowed_values is True, just ignore it
            elif remove_nonallowed_values:
                continue
            # otherwise raise exception
            else:
                raise ValueError("Encountered non-allowed value in one of the dictionaries! See docstring for options.")
    return new_dict


def invert_list_to_dict(input_list):
    """ Given a list with no duplicates, return a dict mapping the values to list positions: [a,b,c] -> {a:1,b:2,c:3}."""
    if not len(set(input_list)) == len(input_list):
        raise ValueError("Can't reliably invert a list with duplicate elements!")
    return dict([(value,index) for (index,value) in enumerate(input_list)])

def invert_dict_nodups(input_dict):
    """ Given a dict with no duplicate values, return value:key dict: {a:1,b:2]} -> {1:a,2:b}."""
    if not len(set(input_dict.values())) == len(input_dict.values()):
        raise ValueError("This method can't invert a dictionary with duplicate values! Use invert_dict_tolists.")
    return dict([(value,key) for (key,value) in input_dict.iteritems()])

def invert_dict_tolists(input_dict):
    """ Given a dict (duplicate values allowed), return value:key_list dict: {a:1,b:2,c:1]} -> {1:[a,c],2:b}."""
    inverted_dict = defaultdict(lambda: set())
    for key,value in input_dict.iteritems():
        inverted_dict[value].add(key)
    return dict(inverted_dict)      # changing defaultdict to plain dict to avoid surprises

def invert_listdict_nodups(input_dict):
    """ Given a dict with non-overlapping list/set values, return single_value:key dict: {a:[1,2],b:[3]} -> {1:a,2:a,3:b}.
    """
    inverted_dict = {}
    for key,value_list in input_dict.iteritems():
        try:
            for value in value_list:
                if value in inverted_dict:
                    raise ValueError("This method can't invert a dictionary with duplicate values! "
                                     +"Use invert_listdict_tolists.")
                inverted_dict[value] = key
        except TypeError:
            raise ValueError("invert_listdict_nodups expects all input_dict values to be lists/sets/etc!")
    return inverted_dict

def invert_listdict_tolists(input_dict):
    """ Given a dict with list/set values, return single_value:key_list dict: {a:[1,2],b:[2]]} -> {1:[a],2:[a,b]}."""
    inverted_dict = defaultdict(lambda: set())
    for key,value_list in input_dict.iteritems():
        try:
            for value in value_list:
                inverted_dict[value].add(key)
        except TypeError:
            raise ValueError("invert_listdict_tolists expects all input_dict values to be lists/sets/etc!")
    return dict(inverted_dict)      # changing defaultdict to plain dict to avoid surprises

# MAYBe-TODO refactor to avoid code duplication between the *_tolists and *_nodups pairs above?

class keybased_defaultdict(defaultdict):
    """ A defaultdict equivalent that passes the key as an argument to the default-value factory, instead of no argumnets.
    Normal defaultdict(int)[9] is 0, because no argument is passed to the int function, and int() is 0.  
    On the other hand keybased_defaultdict(int)[9] would be 9, keybased_defaultdict(bool)[9] would be True, etc.  """
    # from http://stackoverflow.com/questions/2912231/is-there-a-clever-way-to-pass-the-key-to-defaultdicts-default-factory
    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        else:
            value = self[key] = self.default_factory(key)
            return value

### Useful class mix-ins


class FrozenClass(object):
    """ Class that allows prevention of adding new attributes at some point after creation - NOT full immutability.

    Use by inheriting from this, and then running self._freeze() at the end of __init__, like this:
        class Test_freezing(FrozenClass):
            def __init__(self, x, y):
                self.x = x
                self.y = y
                self._freeze() # no new attributes after this point.
        a = Test_freezing(2,3)
        a.x = 10    # can modify existing attributes after creation
        a.z = 10    # fails - cannot add NEW atributes after creation
    """
    # by Jochen Ritzel on Stackoverflow 
    #   (http://stackoverflow.com/questions/3603502/prevent-creating-new-attributes-outside-init)
    # This has to be a new-style class (i.e. inheriting from object), old-style classes don't seem to have __setattr__?
    # MAYBE-TODO could/should this be done as a decorators instead?

    __isfrozen = False

    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError( "%r is a frozen class" % self )
        object.__setattr__(self, key, value)
    
    def _freeze(self):
        self.__isfrozen = True


def reversibly_immutable(wrapped_class):
    """ NOT FINISHED!!!
    
    DECORATOR: gives a class a make_immutable() and make_mutable() method; also makes hashable if _hash() is defined.

    _____
    """
    # MAYBE-TODO not finished!!!  see ~/experiments/mutant_pool_screens/mutant_deepseq_analysis/code/mutant_analysis_classes.py Insertion_position class for one implementation of this - I think it's too class-specific for a general decorator to be easy.

    if not isinstance(object,wrapped_class):
        raise ValueError("Can only apply the reversibly_immutable decorator to new-style classes based on object!")

    def changes_not_allowed(self, *args, **kwargs):
        raise AttributeError("'%s' object is currently immutable, can't make changes!"%type(self))

    # TODO what methods besides __setitem__/__delitem__ do we need to override?  Probably depends on class... COMPLICATED!!
    # MAYBE-TODO this overrides the general methods like __setitem__ etc, but what if the class has other custom methods that modify things? 
    def make_immutable(self):
        # TODO implement removing mutability methods!
        raise NotImplementedError()
        self.__setitem__ = changes_not_allowed
        self.__delitem__ = changes_not_allowed
        # make hashable by setting __hash__ to the private _hash method
        if hasattr(self,'_hash'):
            self.__hash__ = self._hash

    def make_mutable(self):
        # TODO implement re-adding mutability methods!
        raise NotImplementedError()
        # since class is mutable now, it can't be hashable - remove __hash__, 
        #  but keep private _hash method in case we want to make it immutable/hashable again
        if hasattr(self,'__hash__'):
            del self.__hash__

    wrapped_class.make_immutable = make_immutable
    wrapped_class.make_mutable = make_mutable
    return wrapped_class

    
######################################## FILE READING/WRITING, COMMAND-LINE STUFF  ########################################

### Read various file types

def read_two_column_file(filename,numerical_values=True):
    """ Read in a two-column (name,value) tab-separated file (ignore #-comment lines), return data:float(value) dict. """
    data_dict = {}
    for line in open(filename):
        if line[0]=='#':    continue
        name,value = line.strip().split('\t')
        if numerical_values:    data_dict[name] = float(value)
        else:                   data_dict[name] = value
    return data_dict

def read_tab_separated_file(filename, ignore_comments=True, separator='\t'):
    """ Read in a tab-separated file (ignore #-comment lines), return a list for each column. """
    for line in open(filename):
        if ignore_comments and line[0]=='#':    continue
        fields = line.strip('\n').split(separator)
        try:
            for i in range(len(fields)): data_list[i].append(fields[i])
        except NameError:
            data_list = [[x] for x in fields]
        if not len(fields)==len(data_list):
            sys.exit("Error: Found line with %s columns, but first line had %s columns! (line %s)"%(len(fields), 
                                                                                                    len(data_list), line))
    return data_list
    # TODO add to unit-tests? Or some kind of test.

def read_tab_separated_file_with_headers(filename, ID_column=0, ignore_comments=True, all_values_are_numbers=False):
    """ Read a tab-separated file with column headers in the first line; ignore comment lines (starting with #) unless specified otherwise. Return a list of IDs from ID_column (in order), and a column_header:column_data dictionary, with one entry per non-ID column, where column_data is an ID:value dictionary with values taken from the column with that header."""
    list_of_column_data = read_tab_separated_file(filename,ignore_comments)
    ID_list = list_of_column_data.pop(ID_column)    # grab the ID column and remove it from data column list
    ID_list.pop(0)                                  # discard header from ID column
    data_dict_by_header = {}
    for column_data in list_of_column_data:
        # grab header, put the rest in a dictionary by ID, in order
        colheader = column_data.pop(0)
        if all_values_are_numbers:
            data_dict_by_header[colheader] = dict([(ID_list[i],float(column_data[i])) for i in range(len(column_data))])
        else:
            data_dict_by_header[colheader] = dict([(ID_list[i],column_data[i]) for i in range(len(column_data))])
    return ID_list, data_dict_by_header
    # TODO add to unit-tests? Or some kind of test.


### Writing to files

# DECORATOR
def replaces_infile_with_outfile(function_to_be_decorated):
    """ DECORATOR, takes a function that takes an infile and outfile and makes it replace infile with outfile instead. 
    The decorated function must still have an argument named outfile, but its value will never be used."""
    def wrapped_function(infile, *args, **kwargs):
        outfile = '.__tmp__'+infile
        kwargs['outfile'] = outfile
        return_val = function_to_be_decorated(infile, *args, **kwargs)
        if os.path.exists(infile):
            os.remove(infile)
        os.rename(outfile,infile)
        return return_val
    return wrapped_function
    # MAYBE-TODO is the way I'm doing this decorator really the best?  Not bad, but requires adding the fake "outfile" argument to all the functions... If I was somehow the function_to_be_decorated's variable list directly instead (with f.__dict__ or f.__setattr__), that wouldn't be an issue... Ask on stackoverflow?
    # TODO add to unit-tests? Or some kind of test.

def save_line_list_as_file(line_list, filename, header="", add_newlines=True):
    """ Given a list of lines, a filename, and an optional header, open file, write header and all lines, close. """
    line_end = "\n" if add_newlines else ""
    # using the with-as syntax automatically closes the file afterward
    with open(filename,'w') as OUTFILE:
        if header:               OUTFILE.write(header+line_end)
        for line in line_list:   OUTFILE.write(line+line_end)

def write_header_data(OUTFILE,options=None):
    """ Print general script run data (command, path, date/time, full optparse options, etc) to given open file object."""
    import pwd,time,socket
    OUTFILE.write("# Command line this file was generated with: %s\n"%(' '.join(sys.argv)))
    OUTFILE.write("# Path: %s\n"%(os.getcwd()))
    OUTFILE.write("# Date: %s,\t\tUser: %s,\t\tSystem: %s\n"%(time.ctime(), pwd.getpwuid(os.getuid())[0], 
                                                              socket.gethostname()))
    if options:     OUTFILE.write("# Full options: %s\n"%options)

def print_text_from_file(infile, OUTFILE=None, printing=True, add_newlines=0):
    """ Write all text from infile to OUTFILE (if not None), also print to stdout if printing is set. 
    Return line counts.  Infile should be a filename; OUTFILE should be an open file object. """
    line_count = 0
    for line in open(infile):
        if OUTFILE is not None:     OUTFILE.write(line)
        if printing:                print line,
        line_count += 1
    if add_newlines:
        if OUTFILE is not None:     OUTFILE.write('\n'*add_newlines)
        if printing:                print '\n'*add_newlines,
    return line_count


### Command/line utilities (running processes, getting output, etc)

def run_command_and_print_info(command, LOGFILE=None, printing=True, shell=True, program_name=None):
    """ Run command using subprocess.call; first print a line describing that to LOGFILE and/or stdout.
    The shell arg to subprocess.call is given by shell; LOGFILE should be an open file object; 
    program_name is only used for printing, and the first word of the command will be used by default. """
    if program_name is None:
        program_name = command.split(' ')[0]
    output = "### Running %s: %s"%(program_name, command)
    if LOGFILE is not None:     LOGFILE.write(output+'\n')
    if printing:                print output
    subprocess.call([command], shell=shell)


######################################## NUMERIC DATA MANIPULATION ########################################

def int_or_float(x):
    """ Return int(x) if it's equal to x, else return x (useful to convert whole-number float to int). """
    if int(x) == x:     return int(x)
    else:               return x

### Get rid of nan/inf numbers singly or in lists/dicts, replace by input

def clean_number(val,replace_NaN,replace_Inf,replace_NegInf,make_positive=False):
    """ Replace NaN/Inf/NegInf value with arguments provided, otherwise return input. Also make val positive if asked. """
    from numpy import isnan,isinf,isneginf
    if isnan(val):            return replace_NaN
    elif isinf(val):          return replace_Inf
    elif isneginf(val):       return replace_NegInf
    elif make_positive and val<=0:  return replace_NegInf
    else:                           return val
    # TODO add to unit-tests

def clean_data(input_data,replace_NaN,replace_Inf=None,replace_NegInf=None,make_positive=False):
    """ Take a list/tuple/set/dict, return same with NaN/Inf/negative values replaced with the respective arguments. """
    from numpy import isnan,isinf,isneginf
    # if list/tuple/set, just treat as a list and transform back at the end
    if type(input_data) in (list,tuple,set):
        input_list = list(input_data)
    # if dictionary, take the key and value lists, clean the value list, then transform back to dict at the end
    elif type(input_data)==dict:
        key_list = input_data.keys()
        input_list = [input_data[x] for x in key_list]
    else:
        raise ValueError("input_data argument must be a list, tuple, set or dictionary")
    # if replace_NegInf/replace_Inf weren't set, set them to lower/higher than current min/max value
    if replace_NegInf==None:    
        if make_positive:   replace_NegInf = min([x for x in input_list if x>0 and not isnan(x)])
        else:               replace_NegInf = min([x for x in input_list if not isneginf(x) and not isnan(x)])
        if replace_NegInf<0:    replace_NegInf *=1.5
        else:           replace_NegInf *= 0.5
    if replace_Inf==None:    replace_Inf = 1.5 * max([x for x in input_list if not isinf(x) and not isnan(x)])
    # use clean_number to make a cleaned value list
    new_list = [clean_number(x,replace_NaN,replace_Inf,replace_NegInf,make_positive) for x in input_list]
    # output - switch back to original input type
    if type(input_data)==list:      return new_list
    elif type(input_data)==tuple:   return tuple(new_list)
    elif type(input_data)==set:     return set(new_list)
    elif type(input_data)==dict:    return dict(zip(key_list,new_list))
    # TODO add to unit-tests

def clean_data_remove(input_data,remove_NaN=True,remove_Inf=True,remove_NegInf=True,remove_negative=False,remove_zero=False):
    """ Take a list/tuple/set/dict, return same with NaN/Inf/NegInfo/negative values removed accordint to arguments. """
    from numpy import isnan,isinf,isneginf
    def bad_value(x):
        if remove_NaN and isnan(x):         return True
        elif remove_Inf and isinf(x):       return True
        elif remove_NegInf and isneginf(x): return True
        elif remove_negative and x<0:       return True
        elif remove_zero and x==0:          return True
        return False
    # if list/tuple/set, just treat as a list and transform back at the end
    if type(input_data) in (list,tuple,set):
        new_data = [x for x in input_data if not bad_value(x)]
        if type(input_data)==list:      return new_data
        elif type(input_data)==tuple:   return tuple(new_data)
        elif type(input_data)==set:     return set(new_data)
    # if dictionary, have to use a different removal method
    elif type(input_data)==dict:
        return dict([(k,x) for (k,x) in input_data.iteritems() if not bad_value(x)])
    else:
        raise ValueError("input_data argument must be a list, tuple, set or dictionary")
    # TODO add to unit-tests


### cutoffs, ranges, sets, etc 

def parse_cutoffs_into_ranges(cutoff_list, if_overlapping):
    """ Given a cutoff_list [0,2,5], return [(0,inf),(2,inf),(5,inf)] if if_overlapping, else [(0,2),(2,5),(5,inf)]. """ 
    if not if_overlapping:
        return [(x,float('inf')) for x in cutoff_list]
    ranges = []
    for i in range(len(cutoff_list)-1):
        ranges.append((cutoff_list[i],cutoff_list[i+1]))
    ranges.append((cutoff_list[-1],float('inf')))
    return ranges
    # TODO add to unit-tests

def get_sets_from_cutoffs(value_dict, value_ranges):
    """ Given a (name:val) dictionary and a list of (min,max) value ranges, return a (range:names_in_range) dict. """
    ranges_to_sets = {}
    for (vmin,vmax) in value_ranges:
        ranges_to_sets[(vmin,vmax)] = [name for (name,val) in value_dict.items() if val>=vmin and val<vmax]
    return ranges_to_sets
    # TODO add to unit-tests

### local min/max finding

def find_local_maxima_by_width(data, N_surrounding_points=1, include_start_end=True):
    """ Return list of local maxima (value,index) tuples; deals with noise by requiring N lower values on each side.

    Go over the data list searching for a point that's higher than N_surrounding_points on both sides;
     return the list of (value,index) tuples for all such points.
    The point of the N_surrounding_points is to filter out noise by excluding local maxima that are just due to small 
     random changes in what without the noise would be a monotonic line.  There are other approaches to noise problems.

    Dealing with runs of identical local maximum values: the program requires N_surrounding_points on both sides 
      of the entire run of identical values, and only returns the first value/index of the identical run.

    The include_start_end argument governs how to deal with the first and last N_surrounding_points of data:
      - if False, these points will be ignore, since they cannot have N_surrounding_points on both sides
      - if True, these points will be used - the N_surrounding_points requirement will be relaxed to allow fewer
                                                surrounding points if they go up against the data start/end.
    """
    # I also have a full command-line program with options that's a wrapper around this function:
    #    find_local_maxima.py in ~/experiments/other_projects/local_maxima_detector_for_Ute
    if include_start_end:
        padding = [min(data)] * N_surrounding_points
        data = padding + data + padding
    else:   padding = []
    local_maxima = []
    for i in range(N_surrounding_points, len(data)-N_surrounding_points):
        curr_value = data[i]
        # if there are multiple identical values, ignore any but the first
        if data[i-1] == data[i]:
            continue
        # grab N_surrounding_points values before the current value
        values_before_curr = data[i-N_surrounding_points : i]
        # grab N_surrounding_points values after the current value, IGNORING any identical ones
        N_plateau_values = 0
        for x in data[i+1:]:
            if x != curr_value:
                break
            N_plateau_values += 1
        values_after_curr = data[i+1+N_plateau_values : i+1+N_plateau_values+N_surrounding_points]
        # if the current value is higher than N_surrounding_points on both sides, add it to local_maxima list
        if curr_value > max(values_before_curr) and curr_value > max(values_after_curr):
            # (if the data was padded, subtract the padding size from the index to get index in original data)
            curr_index = i - len(padding)
            local_maxima.append((curr_value, curr_index))
    return local_maxima

# MAYBE-TODO an alternative implementation would be find_local_maxima_by_height, where we'd take any point that was at least Kx higher than the surrounding M or fewer points, even if M was just 2.  That would also filter out noise reasonably well, assuming the noise is small in amplitude.
# MAYBE-TODO for more flexible dealing with noise, there may be something in scipy... See http://stackoverflow.com/questions/1713335/peak-finding-algorithm-for-python-scipy

### moving average and moving median, plus some help functions (anything starting with _ won't be imported)

def _check_input(data,window_size):
    if window_size<1 or window_size>len(data):   raise ValueError, "window size must be between 1 and data length."

def _make_window_indices(data,window_size):
    half_window = int(window_size/2)
    index_values = [i+half_window for i in range(len(data)-window_size+1)]
    return index_values


def moving_average(data,window_size=10,return_indices=False):
    """ Return the moving average of the input data with the given window size. 

    Optionally also returns the window center indices for plotting. 
    Reasonably efficient implementation, doesn't calculate the whole average for each window.
    """
    _check_input(data,window_size)
    first_average = sum(data[:window_size])/window_size
    average_values = [first_average]
    for i in range(len(data)-window_size):
        # for each new average, just take the previous one, subtract the oldest data value and add the new data value
        new_average = average_values[-1] + (data[i+window_size]-data[i])/window_size
        average_values.append(new_average)
    if return_indices:  return average_values, _make_window_indices(data,window_size)
    else:               return average_values 
    # TODO add to unit-tests


def moving_median(data,window_size=10,return_indices=False):
    """ Return the moving median of the input data with the given window size. 

    Optionally also returns the window center indices for plotting. 
    """
    _check_input(data,window_size)
    from numpy import median
    median_values = []
    for i in range(len(data)-window_size+1):
        median_values.append(median(data[i:i+window_size]))
    if return_indices:  return median_values, _make_window_indices(data,window_size)
    else:               return median_values 
    # TODO add to unit-tests


### Other specialized functions 

def split_into_N_sets_by_counts(ID_counts, N):
    """ Given an ID:count dictionary, return a list of sets of IDs with total counts balanced between the sets. """
    # make a sorted (high to low) list of (count,ID) tuples
    counts_IDs = sorted([(count,ID) for (ID,count) in ID_counts.iteritems()], reverse=True)
    output_counts_sets = [[0,set()] for i in range(N)]
    # now go over all IDs, adding an ID (and the corresponding count) to the smallest set each time
    for (count,ID) in counts_IDs:
        output_counts_sets[0][1].add(ID)
        output_counts_sets[0][0] += count
        output_counts_sets.sort()
    return [ID_set for [count,ID_set] in output_counts_sets]



######################################## TESTS FOR THIS FILE ########################################

class Testing_everything(unittest.TestCase):
    """ Testing all functions/classes.etc. """

    def test__compare_lists_unordered(self):
        from itertools import permutations
        for input_list in [[1,2,3], [True,False,True], ['a','bb',''], [1,True,'str',321314,None]]:
            assert compare_lists_unordered(input_list, input_list) == True
            assert compare_lists_unordered(input_list, input_list*2) == False
            assert compare_lists_unordered(input_list, []) == False
            for permuted_list in permutations(input_list, len(input_list)):
                assert compare_lists_unordered(input_list, permuted_list) == True
                assert compare_lists_unordered(input_list, permuted_list[:-1]) == False
                assert compare_lists_unordered(input_list[:-1], permuted_list) == False
            for permuted_list in permutations(input_list, len(input_list)-1):
                assert compare_lists_unordered(input_list, permuted_list) == False

    def test__reduce_dicts_to_overlaps(self):
        d1 = {1:1}
        d2 = {1:2, 2:2}
        d3 = {2:3, 3:3}
        assert reduce_dicts_to_overlaps([d1,d2]) == [{1:1},{1:2}] 
        assert reduce_dicts_to_overlaps([d2,d3]) == [{2:2},{2:3}] 
        assert reduce_dicts_to_overlaps([d1,d3]) == [{},{}] 
        assert reduce_dicts_to_overlaps([d1,d2,d3]) == [{},{},{}] 
        assert reduce_dicts_to_overlaps([]) == [] 

    def test__count_list_values(self):
        assert count_list_values([]) == {}
        assert count_list_values([10,12,11]) == {10:1, 12:1, 11:1}
        assert count_list_values([10,12,10]) == {10:2, 12:1}
        assert count_list_values([1]*100) == {1:100}
        assert count_list_values(['a',None,11]) == {'a':1, None:1, 11:1}
        assert count_list_values([None,None]) == {None:2}

    def test__add_dicts_of_ints(self):
        ### Basic int-value-only cases should be the same regardless of options
        for recursive in True,False:
            for allow_nonconflicting in True,False:
                for remove_nonallowed in True,False:
                    kwargs = {'recursive':recursive, 'allow_nonconflicting_values':allow_nonconflicting, 
                              'remove_nonallowed_values':remove_nonallowed}
                    for dict2 in [{}, {1:1}, {'a':1, 'b':100, 'c':-200}]:
                        assert add_dicts_of_ints({}, dict2, **kwargs) == dict2
                    assert add_dicts_of_ints({1:1}, {1:2}, **kwargs) == {1:3}
                    assert add_dicts_of_ints({1.00:1}, {1.00:2}, **kwargs) == {1.00:3}
                    assert add_dicts_of_ints({1:1}, {2:2}, **kwargs) == {1:1, 2:2}
                    assert add_dicts_of_ints({'one':1}, {'two':2}, **kwargs) == {'one':1, 'two':2}
                    assert add_dicts_of_ints({1:1, 2:1}, {2:2, 3:2}, **kwargs) == {1:1, 2:3, 3:2}
                    assert add_dicts_of_ints({(1,):1, (2,):1}, {(2,):2, (3,):2}, **kwargs) == {(1,):1, (2,):3, (3,):2}
        ### Simplest version - no recursion, no non-int values allowed: all "weird" cases raise exceptions 
        kwargs = {'recursive':False, 'allow_nonconflicting_values':False, 'remove_nonallowed_values':False}
        for bad_value in ['a', [], {}, {1:1}, {'a':1}, {1:'a'}, None, True, False, [1,2,3], set([1,2,3]), 'abadf', int]:
            self.assertRaises(ValueError, add_dicts_of_ints, {1:bad_value}, {}, **kwargs)
            self.assertRaises(ValueError, add_dicts_of_ints, {1:bad_value}, {1:bad_value}, **kwargs)
        # if nonallowed values are removed, exceptions aren't raised, the bad values are just ignored
        kwargs = {'recursive':False, 'allow_nonconflicting_values':False, 'remove_nonallowed_values':True}
        for bad_value in ['a', [], {}, {1:1}, {'a':1}, {1:'a'}, None, True, False, [1,2,3], set([1,2,3]), 'abadf', int]:
            # if there's only a bad value, or two of them, the result is empty
            assert add_dicts_of_ints({1:bad_value}, {}, **kwargs) == {}
            assert add_dicts_of_ints({1:bad_value}, {1:bad_value}, **kwargs) == {}
            # same if there's also a good value for a different key - only that is kept
            assert add_dicts_of_ints({1:bad_value}, {2:1}, **kwargs) == {2:1}
            # if there's a bad and good value for the same key, only the good one is kept
            assert add_dicts_of_ints({1:bad_value}, {1:1}, **kwargs) == {1:1}
        ### Allowing recursion: 
        kwargs = {'recursive':True, 'allow_nonconflicting_values':False, 'remove_nonallowed_values':False}
        # still fails for non-dict "bad" values
        for bad_value in ['a', [], {1:'a'}, None, True, False, [1,2,3], set([1,2,3]), 'abadf', int]:
            self.assertRaises(ValueError, add_dicts_of_ints, {1:bad_value}, {}, **kwargs)
            self.assertRaises(ValueError, add_dicts_of_ints, {1:bad_value}, {1:bad_value}, **kwargs)
            self.assertRaises(ValueError, add_dicts_of_ints, {1:bad_value}, {1:1}, **kwargs)
        # if there's a dictionary-type value in one dictionary and nothing in the other, it's kept
        for dict_value in [{}, {1:1}, {'a':1}]:
            assert add_dicts_of_ints({1:dict_value}, {}, **kwargs) == {1:dict_value}
        # if there's a dictionary-type value in both dictionaries, recursively add them up
        assert add_dicts_of_ints({1:{}}, {1:{}}, **kwargs) == {1:{}}
        assert add_dicts_of_ints({1:{1:1}}, {1:{1:1}}, **kwargs) == {1:{1:2}}
        assert add_dicts_of_ints({1:{'a':1}}, {1:{'a':1}}, **kwargs) == {1:{'a':2}}
        assert add_dicts_of_ints({1:1, 2:{2:2, 3:{3:3}}}, {}, **kwargs) == {1:1, 2:{2:2, 3:{3:3}}}
        assert add_dicts_of_ints({1:1, 2:{2:2, 3:{3:3}}}, {1:1, 2:{}}, **kwargs) == {1:2, 2:{2:2, 3:{3:3}}}
        assert add_dicts_of_ints({1:1, 2:{2:2, 3:{3:3}}}, {1:1, 2:{4:4}}, **kwargs) == {1:2, 2:{2:2, 3:{3:3}, 4:4}}
        # if there's a dictionary-type value in one dictionary and an int in the other, raise an exception
        for dict_value in [{}, {1:1}, {'a':1}]:
            self.assertRaises(ValueError, add_dicts_of_ints, {1:dict_value}, {1:1}, **kwargs)
        self.assertRaises(ValueError, add_dicts_of_ints, {1:1, 2:{2:2, 3:{3:3}}}, {1:1, 2:1}, **kwargs)
        ## if nonallowed values are removed, exceptions aren't raised: the bad values are just ignored
        kwargs = {'recursive':True, 'allow_nonconflicting_values':False, 'remove_nonallowed_values':True}
        # for bad values that don't involve dictionaries, the outcomes are the same as in the non-recursive version
        for bad_value in ['a', [], None, True, False, [1,2,3], set([1,2,3]), 'abadf', int]:
            assert add_dicts_of_ints({1:bad_value}, {}, **kwargs) == {}
            assert add_dicts_of_ints({1:bad_value}, {1:bad_value}, **kwargs) == {}
            assert add_dicts_of_ints({1:bad_value}, {2:1}, **kwargs) == {2:1}
            assert add_dicts_of_ints({1:bad_value}, {1:1}, **kwargs) == {1:1}
        # if the two values are mismatched but they're both "good" (like an int and a dict), both are ignored
        for dict_value in [{}, {1:1}, {'a':1}]:
            assert add_dicts_of_ints({1:dict_value}, {1:1}, **kwargs) == {}
        assert add_dicts_of_ints({1:1, 2:{2:2, 3:{3:3}}}, {1:1, 2:1}, **kwargs) == {1:2}
        ## another special case: if a value is a dictionary but it's weird, it gets recursively passed down:
        #   if remove_nonallowed_values is False, that leads to a ValueError a level later (which was tested above), 
        #   but if it's True, it leads to the dictionary being kept in an empty form (since it contained a bad value)
        weird_value = {1:'a'}
        assert add_dicts_of_ints({1:weird_value}, {}, **kwargs) == {1:{}}
        assert add_dicts_of_ints({1:weird_value}, {1:weird_value}, **kwargs) == {1:{}}
        assert add_dicts_of_ints({1:weird_value}, {2:1}, **kwargs) == {1:{}, 2:1}
        assert add_dicts_of_ints({1:weird_value}, {1:1}, **kwargs) == {}
        ### With allow_nonconflicting_values turned on:
        kwargs = {'recursive':False, 'allow_nonconflicting_values':True, 'remove_nonallowed_values':False}
        # this time if the "bad" value is only in one of the dictionaries and nothing in the other, 
        #  or is identical in both dictionaries, it's kept. 
        for bad_value in ['a', [], {}, {1:1}, {'a':1}, {1:'a'}, None, True, False, [1,2,3], set([1,2,3]), 'abadf', int]:
            assert add_dicts_of_ints({1:bad_value}, {}, **kwargs) == {1:bad_value}
            assert add_dicts_of_ints({1:bad_value}, {1:bad_value}, **kwargs) == {1:bad_value}
            assert add_dicts_of_ints({1:bad_value}, {2:1}, **kwargs) == {1:bad_value, 2:1}
        # but if there's a "good" value in one dictionary and a "bad" one in the other, an exception is still raised.
        for bad_value in ['a', [], {}, {1:1}, {'a':1}, {1:'a'}, None, False, [1,2,3], set([1,2,3]), 'abadf', int]:
            self.assertRaises(ValueError, add_dicts_of_ints, {1:bad_value}, {1:1}, **kwargs)
            # except unfortunately this doesn't work right for True, because 1==True... TODO try to fix that somehow?
        ## if nonallowed values are removed, the first cases are the same, in the last case the good value is kept
        kwargs = {'recursive':False, 'allow_nonconflicting_values':True, 'remove_nonallowed_values':True}
        for bad_value in ['a', [], {}, {1:1}, {'a':1}, {1:'a'}, None, True, False, [1,2,3], set([1,2,3]), 'abadf', int]:
            assert add_dicts_of_ints({1:bad_value}, {}, **kwargs) == {1:bad_value}
            assert add_dicts_of_ints({1:bad_value}, {1:bad_value}, **kwargs) == {1:bad_value}
            assert add_dicts_of_ints({1:bad_value}, {2:1}, **kwargs) == {1:bad_value, 2:1}
            assert add_dicts_of_ints({1:bad_value}, {1:1}, **kwargs) == {1:1}
        ### Recursive takes priority over allow_nonconflicting_values if both apply:
        # if recursive is False, if both values are identical dictionaries, they always get passed on unchanged.
        kwargs = {'recursive':False, 'allow_nonconflicting_values':True, 'remove_nonallowed_values':False}
        assert add_dicts_of_ints({1:{1:2}}, {1:{1:2}}, **kwargs) == {1:{1:2}}
        assert add_dicts_of_ints({1:{1:'a'}}, {1:{1:'a'}}, **kwargs) == {1:{1:'a'}}
        # if recursive is True, the same still happens with bad-valued dictionaries, 
        #  but for int-valued dictionaries the new value is the sum of the originals instead of a copy
        kwargs = {'recursive':True, 'allow_nonconflicting_values':True, 'remove_nonallowed_values':False}
        assert add_dicts_of_ints({1:{1:2}}, {1:{1:2}}, **kwargs) == {1:{1:4}}
        assert add_dicts_of_ints({1:{1:'a'}}, {1:{1:'a'}}, **kwargs) == {1:{1:'a'}}

    def test__invert_list_to_dict(self):
        assert invert_list_to_dict([]) == {}
        assert invert_list_to_dict([10,12,11]) == {10:0, 12:1, 11:2}
        assert invert_list_to_dict(['a',None,11]) == {'a':0, None:1, 11:2}
        self.assertRaises(ValueError, invert_list_to_dict, [10,11,10])
        self.assertRaises(ValueError, invert_list_to_dict, [None,None])

    def test__invert_dict_nodups(self):
        assert invert_dict_nodups({}) == {}
        assert invert_dict_nodups({1:2,3:4}) == {2:1,4:3}
        self.assertRaises(ValueError, invert_dict_nodups, {1:2,3:2})

    def test__invert_dict_tolists(self):
        assert invert_dict_tolists({}) == {}
        assert invert_dict_tolists({1:2,3:4}) == {2:set([1]),4:set([3])}
        assert invert_dict_tolists({1:2,3:2}) == {2:set([1,3])}

    def test__invert_listdict_nodups(self):
        assert invert_listdict_nodups({}) == {}
        assert invert_listdict_nodups({1:[2],3:[4]}) == {2:1,4:3}
        assert invert_listdict_nodups({1:[2],3:[4,6]}) == {2:1,4:3,6:3}
        assert invert_listdict_nodups({1:[2,8],3:[4,6]}) == {2:1,8:1,4:3,6:3}
        # no duplicates
        self.assertRaises(ValueError, invert_listdict_nodups, {1:[2,2],3:[4]})
        self.assertRaises(ValueError, invert_listdict_nodups, {1:[2,4],3:[2]})
        # values must be lists/sets
        self.assertRaises(ValueError, invert_listdict_nodups, {1:2,3:2})

    def test__invert_listdict_tolists(self):
        assert invert_listdict_tolists({}) == {}
        assert invert_listdict_tolists({1:[2],3:[4]}) == {2:set([1]),4:set([3])}
        assert invert_listdict_tolists({1:[2],3:[2]}) == {2:set([1,3])}
        assert invert_listdict_tolists({1:[2,4],3:[2]}) == {2:set([1,3]),4:set([1])}
        self.assertRaises(ValueError, invert_listdict_tolists, {1:2,3:2})

    def test__keybased_defaultdict(self):
        D_nodefault = keybased_defaultdict(None)
        self.assertRaises(KeyError, lambda: D_nodefault[1])
        D_constantdefault = keybased_defaultdict(lambda x: 0)
        assert D_constantdefault[1] == 0
        assert D_constantdefault[2] == 0
        D_variabledefault = keybased_defaultdict(lambda x: 2*x)
        assert D_variabledefault[1] == 2
        assert D_variabledefault[2] == 4

    def test__FrozenClass(self):
        class Test_freezing(FrozenClass):
            def __init__(self, x, y):
                self.x = x
                self.y = y
                self._freeze() # no new attributes after this point.
        a = Test_freezing(2,3)
        assert a.x == 2
        a.x = 10    # can modify existing attributes after creation - shouldn't raise an exception
        assert a.x == 10
        # doing a.z = 10 should fail - cannot add NEW atributes after creation
        #  testing this two ways: with __setattr__ as a function, and with writing a test function to explicitly 
        #  test the "a.z = 1" statement (which of course should use __setattr__ anyway, but may as well check)
        self.assertRaises(TypeError, a.__setattr__, 'z', 10)
        def test_function(obj,val):  
            obj.z = val
        self.assertRaises(TypeError, test_function, a, 10)

    def test__int_or_float(self):
        for bad_input in ['a', [], [1,2,3]]:
            self.assertRaises((ValueError,TypeError), int_or_float, bad_input)
        assert int_or_float(3) == 3
        assert int_or_float(3.0) == 3
        assert int_or_float(3.5) == 3.5

    def test__find_local_maxima_by_width(self):
        # basic functionality - find the local maximum
        for N_surrounding_points in [1,2,3]:
            assert find_local_maxima_by_width([1,2,3,4,3,2,1], N_surrounding_points) == [(4,3)]
            assert find_local_maxima_by_width([1,2,3,4,3,2,1,2,3,4,3,2,1], N_surrounding_points) == [(4,3),(4,9)]
        # require at least N_surrounding_points lower points on both sides of the maximum
        assert find_local_maxima_by_width([1,2,1,2,1,2,1], N_surrounding_points=1) == [(2,1),(2,3),(2,5)]
        for N_surrounding_points in [2,3]:
            assert find_local_maxima_by_width([1,2,1,2,1,2,1], N_surrounding_points) == []
        # require N_surrounding_points on either end 
        assert find_local_maxima_by_width([1,1,2,1,2], N_surrounding_points=2) == []
        assert find_local_maxima_by_width([1,1,2,1,2], N_surrounding_points=1) == [(2,2),(2,4)]
        assert find_local_maxima_by_width([2,1,2,1,1], N_surrounding_points=2) == []
        assert find_local_maxima_by_width([2,1,2,1,1], N_surrounding_points=1) == [(2,0),(2,2)]
        # how to deal with points in the first/last N_surrounding_points depends on include_start_end:
        assert find_local_maxima_by_width([2,1,1,2], N_surrounding_points=1, include_start_end=False) == []
        assert find_local_maxima_by_width([2,1,1,2], N_surrounding_points=1, include_start_end=True) == [(2,0),(2,3)]
        # when there are adjacent points with same value, return either one
        assert find_local_maxima_by_width([1,2,2,1],1) in ( [(2,1)], [(2,2)] )

    def test__split_into_N_sets_by_counts(self):
        input1 = {'a':1000}
        for N in range(1,10):
            assert compare_lists_unordered(split_into_N_sets_by_counts(input1,N), 
                                           [set(['a'])] + [set() for i in range(N-1)])

        input2 = {'a':1002, 'b':1001, 'c':1000}
        assert compare_lists_unordered(split_into_N_sets_by_counts(input2,1), [set(['a','b','c'])])
        assert compare_lists_unordered(split_into_N_sets_by_counts(input2,2), [set(['a']),set(['b','c'])])
        for N in range(3,10):
            assert compare_lists_unordered(split_into_N_sets_by_counts(input2,N), 
                                           [set(['a']), set(['b']), set(['c'])] + [set() for i in range(N-3)])

        input3 = {'a':5, 'b':4, 'c':3, 'd':2, 'e':1}
        assert compare_lists_unordered(split_into_N_sets_by_counts(input3,1), [set(['a','b','c','d','e'])])
        assert compare_lists_unordered(split_into_N_sets_by_counts(input3,2), [set(['a','d','e']), set(['b','c'])])
        assert compare_lists_unordered(split_into_N_sets_by_counts(input3,3), [set(['a']), set(['b','e']), set(['c','d'])])
        assert compare_lists_unordered(split_into_N_sets_by_counts(input3,4), 
                                       [set(['a']), set(['b']), set(['c']), set(['d','e'])])
        assert compare_lists_unordered(split_into_N_sets_by_counts(input3,5), 
                                       [set(['a']), set(['b']), set(['c']), set(['d']), set(['e'])])
        assert compare_lists_unordered(split_into_N_sets_by_counts(input3,6), 
                                       [set(['a']), set(['b']), set(['c']), set(['d']), set(['e']), set()])
        assert compare_lists_unordered(split_into_N_sets_by_counts(input3,7), 
                                       [set(['a']), set(['b']), set(['c']), set(['d']), set(['e']), set(), set()])
        assert compare_lists_unordered(split_into_N_sets_by_counts(input3,8), 
                                       [set(['a']), set(['b']), set(['c']), set(['d']), set(['e']), set(), set(), set()])

    # TODO add tests for everything else


if __name__=='__main__':
    """ If module is ran directly, run tests. """
    print "This is a module for import by other programs - it doesn't do anything on its own.  Running tests..."
    unittest.main()
