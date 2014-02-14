import math
cimport cython
from itertools import tee
from libc.stdlib cimport malloc, free
import random

cdef extern from "WellingtonC.h":
    struct tuple2:
         float * fpscores
         unsigned int * mles
    double bdtrc(int, int, float)
    double bdtr(int, int, float)
    tuple2 * wellington(unsigned int *, unsigned int *, unsigned int,unsigned int *, unsigned int,unsigned int *, unsigned int)

def logsf(a,b,c):
    return bdtrc(a, b, c)
def logcdf(a,b,c):
    return bdtr(a,b,c)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def window(iterable, int size):
    """
    Takes reads list (iterable) and returns reads list, each of length size, of rolling windows.
    >>> [i for i in window(range(0,12,2), 3)]
    [(0, 2, 4), (2, 4, 6), (4, 6, 8), (6, 8, 10)]
    """
    iters = tee(iterable, size)
    cdef unsigned int i
    for i in range(1, size):
        for each in iters[i:]:
            next(each, None)
    return zip(*iters)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef percentile(list N, float percent):
    """
    Find the percentile of a list of values.

    @parameter N - is a list of values.
    @parameter percent - a float value from 0.0 to 1.0.

    @return - the percentile of the values as a float
    """
    cdef float k,d0,d1
    cdef unsigned int f, c
    if not N:
        return None
    N = sorted(N)
    k = (len(N)-1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return float(N[int(k)])
    d0 = N[int(f)] * (c-k)
    d1 = N[int(c)] * (k-f)
    return float(d0+d1)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def calculate(FDR,list forwardArray,list backwardArray, footprint_sizes, shoulder_sizes, bonferroni_factor):

    if FDR:
        random.shuffle(forwardArray)
        random.shuffle(backwardArray)

    cdef unsigned asize = len(forwardArray)


    #Copy the forward and reverse reads and the shoulder and fp sizes into the C arrays
    cdef unsigned int *farr  = <unsigned int *> malloc(asize * sizeof(unsigned int))
    cdef unsigned int *rarr  = <unsigned int *> malloc(asize * sizeof(unsigned int))
    cdef unsigned int *shsize  = <unsigned int *> malloc(len(shoulder_sizes) * sizeof(unsigned int))
    cdef unsigned int *fpsize  = <unsigned int *> malloc(len(footprint_sizes) * sizeof(unsigned int))
    for k in range(asize):
        farr[k] = forwardArray[k]
        rarr[k] = backwardArray[k]

    for index, size in enumerate(shoulder_sizes):
        shsize[index] = size
    for index, size in enumerate(footprint_sizes):
        fpsize[index] = size

    #Farm off the Computation to C
    cdef tuple2 * test = wellington(farr,rarr,asize,shsize,len(shoulder_sizes),fpsize,len(footprint_sizes))
    #Free the arrays
    free(farr)
    free(rarr)

    #Push the values from C arrays back into Python lists
    cdef list m,f
    m, f  = [], []
    cdef unsigned int j
    for j in range(asize):
        if bonferroni_factor:
            f.append(min(0,test.fpscores[j] - bonferroni_factor))
        else:
            f.append(test.fpscores[j])
        m.append(test.mles[j])

    #Clean everything up
    free(test.mles)
    free(test.fpscores)
    free(test)

    #Return the Scores and the FP parameters
    return(f,m)

