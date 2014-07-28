cimport cython
from libc.math cimport ceil
from libc.stdlib cimport malloc, free
import random
random.seed("Congratulations on reading the source code - you win a prize!")

cdef extern from "WellingtonC.h":
    struct tuple2:
         float * fpscores
         unsigned int * mles
    double bdtrc(int, int, float)
    double bdtr(int, int, float)
    tuple2 * wellington(unsigned int *, unsigned int *, unsigned int,unsigned int *, unsigned int,unsigned int *, unsigned int)
    tuple2 * diff_wellington(unsigned int * f,  unsigned int * r, unsigned int * f2,  unsigned int * r2, unsigned int length, unsigned int * offsets, unsigned int * widths, unsigned int num_offsets,float threshold)

cdef float logsf(int a,int b,float c):
    return bdtrc(a, b, c)
cdef float logcdf(int a,int b,float c):
    return bdtr(a,b,c)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef float percentile(list N, float percent):
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
    f = <int>k
    c = <int>ceil(k)
    if f == c:
        return float(N[int(k)])
    d0 = N[f] * (c-k)
    d1 = N[c] * (k-f)
    return d0+d1

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def diff_calculate(list forwardArray, list backwardArray,list forwardArray2, backwardArray2, footprint_sizes, offset_positions,float threshold):
    cdef unsigned int asize = len(forwardArray)
    cdef unsigned int *farr  = <unsigned int *> malloc(asize * sizeof(unsigned int))
    cdef unsigned int *rarr  = <unsigned int *> malloc(asize * sizeof(unsigned int))
    cdef unsigned int *farr2  = <unsigned int *> malloc(asize * sizeof(unsigned int))
    cdef unsigned int *rarr2  = <unsigned int *> malloc(asize * sizeof(unsigned int))
    cdef unsigned int *fpsize  = <unsigned int *> malloc(len(footprint_sizes) * sizeof(unsigned int))
    cdef unsigned int *offsets  = <unsigned int *> malloc(len(offset_positions) * sizeof(unsigned int))
    for k in range(asize):
        farr[k] = forwardArray[k]
        rarr[k] = backwardArray[k]
        farr2[k] = forwardArray2[k]
        rarr2[k] = backwardArray2[k]

    for index, size in enumerate(footprint_sizes):
        fpsize[index] = size
    for index, size in enumerate(offset_positions):
        offsets[index] = size
    #diff_wellington(unsigned int * f,  unsigned int * r, unsigned int * f2,  unsigned int * r2, unsigned int length, unsigned int * offsets, unsigned int * widths, unsigned int num_offsets)
   # tuple2 * test = diff_wellington(farr,  rarr, farr2,  rarr2, asize, offsets, fpsize, len(footprint_sizes))
    cdef tuple2 * test = diff_wellington(farr,  rarr, farr2,  rarr2, asize, offsets, fpsize, len(footprint_sizes),threshold)
    #Push the values from C arrays back into Python lists
    cdef list m,f
    m, f  = [], []
    cdef unsigned int j
    for j in range(asize):
        f.append(test.fpscores[j])
        m.append(int(test.mles[j]))

    #Clean everything up
    free(test.mles)
    free(test.fpscores)
    free(test)
    free(farr)
    free(rarr)
    free(farr2)
    free(rarr2)
    free(fpsize)
    free(offsets)

    #Return the Scores and the FP parameters
    return(f,m)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def calculate(FDR, list forwardArray_in, list backwardArray_in, footprint_sizes, shoulder_sizes, bonferroni_factor):

    #If we're looking for an FDR, shuffle the data here.

    if FDR:
        forwardArray, backwardArray = forwardArray_in[:], backwardArray_in[:]
        random.shuffle(forwardArray)
        random.shuffle(backwardArray)
    else:
        forwardArray, backwardArray = forwardArray_in, backwardArray_in

    cdef unsigned int asize = len(forwardArray)

    #Copy the forward and reverse reads and the shoulder and fp sizes into the C arrays
    cdef unsigned int * farr  = <unsigned int *> malloc(asize * sizeof(unsigned int))
    cdef unsigned int * rarr  = <unsigned int *> malloc(asize * sizeof(unsigned int))
    cdef unsigned int * shsize  = <unsigned int *> malloc(len(shoulder_sizes) * sizeof(unsigned int))
    cdef unsigned int * fpsize  = <unsigned int *> malloc(len(footprint_sizes) * sizeof(unsigned int))
    for k in range(asize):
        farr[k] = forwardArray[k]
        rarr[k] = backwardArray[k]

    for index, size in enumerate(shoulder_sizes):
        shsize[index] = size
    for index, size in enumerate(footprint_sizes):
        fpsize[index] = size

    #Farm off the Computation to C
    cdef tuple2 * test = wellington(farr,rarr,asize,shsize,len(shoulder_sizes),fpsize,len(footprint_sizes))

    #Push the values from C arrays back into Python lists
    cdef list m,f
    m, f  = [], []
    cdef unsigned int j
    cdef int minlimit = -1000
    for j in range(asize):
        if bonferroni_factor:
            f.append(min(0,test.fpscores[j] - bonferroni_factor))
        else:
            f.append(max(test.fpscores[j],minlimit))
        m.append(int(test.mles[j]))

    #Clean everything up
    free(test.mles)
    free(test.fpscores)
    free(test)
    free(farr)
    free(rarr)
    free(shsize)
    free(fpsize)

    #Return the Scores and the FP parameters
    return(f,m)
