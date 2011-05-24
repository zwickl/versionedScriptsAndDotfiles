#!/usr/bin/env python

def safe_open(filename, flags = "rb"):
    '''
    
    >>> safe_open("ajsdlfj.poo")
    SystemExit: problem opening file ajsdlfj.poo
    
    '''
    if filename is None:
        exit("call safe_open with a filename")
    try:
        f = open(filename, flags)
    except IOError:
        raise IOError("problem opening file %s" % filename)
        #exit("problem opening file %s" % filename)
    except:
        exit("unknown problem opening file %s with flags %s" % ( filename, flags))
    return f

if __name__ == "__main__":
    import doctest
    doctest.testmod()

def unsafe_open(filename, flags = "rb"):
    if filename is None:
        exit("call safe_open with a filename")
    f = open(filename, flags)

