#! /usr/bin/env python2
"""
Testing utilities. 
 -- Weronika Patena, Jonikas lab, Carnegie Institution, July 2011
"""

class TestingError(Exception):
    pass

def call_should_fail(function, args, Error, message="test failed!", NewError=TestingError):
    """ Make sure that funtion(*args) raises Error; work silently if it does, raise NewError(message) if it doesn't."""
    try:            function(*args)
    except Error:   pass
    else:           raise NewError(message)


if __name__ == '__main__':
    print "Testing call_should_fail function..."
    # list(3) SHOULD return TypeError, so call_should_fail on that should work and not raise an exception!
    try:                    call_should_fail(list, [3], TypeError)
    except TestingError:    raise Exception("call_should_fail(list(3)) (invalid arg) should have worked, but didn't!")
    # on the other hand call_should_fail on a valid argument to list should raise an exception.
    try:                    call_should_fail(list, [[1,2,3]], TypeError)
    except TestingError:    pass
    else:                   raise Exception("call_should_fail(list([1,2,3])) (valid arg) should have failed, but didn't!")
    # similarly with int('t') vs int(3)  (and this time passing ValueError as NewError
    try:                call_should_fail(int, ['t'], ValueError, NewError=ValueError)
    except ValueError:   raise Exception("call_should_fail(int('t')) (invalid arg) should have worked, but didn't!")
    try:                call_should_fail(int, ['3'], ValueError, NewError=ValueError)
    except ValueError:   pass
    else:               raise Exception("call_should_fail(int('3')) (valid arg) should have failed, but didn't!")
    print "...DONE"

