import math
cimport cython
from libc.math cimport ceil
from libc.stdlib cimport malloc, free
import random
import matplotlib.pyplot as plt
from scipy.stats import percentileofscore as pct
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
def markus_diff_calculate(list forwardArray, list backwardArray,list forwardArray2, list backwardArray2, list footprint_sizes, list offset_positions,float threshold,plotting =0):
    cdef unsigned int asize = len(forwardArray)
    fp_scores = [101] * asize
    fp_mles = [0] * asize
    cdef unsigned int shoulder_size = 35
    cdef unsigned int centre,fp_size,mle,offset

    for offset, mle in zip(offset_positions,footprint_sizes):
        # I've commented out the expanding and wobbling for now!
        for fp_size in [mle]:# range(mle-2,mle+4,2): #[mle]:
            halffpround = int((fp_size)/2)
            for centre in [offset]:# range(offset-2,offset+4,2):

                mini_array_f  = forwardArray[centre-halffpround-shoulder_size:centre+halffpround+shoulder_size]
                mini_array_r  = backwardArray[centre-halffpround-shoulder_size:centre+halffpround+shoulder_size]

                mini_array_f2 = forwardArray2[centre-halffpround-shoulder_size:centre+halffpround+shoulder_size]
                mini_array_r2 = backwardArray2[centre-halffpround-shoulder_size:centre+halffpround+shoulder_size]

                Forward   = sum(mini_array_f[:shoulder_size])
                cForward  = sum(mini_array_f[shoulder_size:shoulder_size+fp_size])
                cBackward = sum(mini_array_r[shoulder_size:shoulder_size+fp_size])
                Backward  = sum(mini_array_r[shoulder_size+fp_size:shoulder_size+shoulder_size+fp_size])

                #This calculates the FOS score
                #t_score = (((cForward+1.0)/(fp_size *1.0)) / (Forward+1.0/(shoulder_size *1.0))) + ((cBackward+1.0/(fp_size *1.0))/(Backward+1.0/(shoulder_size *1.0)))

                #Note: we're going with Markus' idea here. Not sure how good this will be :/
                t_score = (Forward + Backward) / float((Forward + Backward + cForward + cBackward))

                scores = []
                #Now let's randomly swap shit around
                for bootstrap_iteration in range(10000):
                    decision_vector = [random.randrange(2) for i in mini_array_f]

                    newvec_f = []
                    newvec_r = []
                    newvec_c = []
                    for pos, value in enumerate(decision_vector):
                        if value == 0:
                            newvec_f.append(mini_array_f[pos])
                            newvec_r.append(mini_array_r[pos])
                            newvec_c.append(mini_array_f[pos] + mini_array_r[pos])
                        else:
                            newvec_f.append(mini_array_f2[pos])
                            newvec_r.append(mini_array_r2[pos])
                            newvec_c.append(mini_array_f2[pos] + mini_array_r2[pos])

                    #Now let's calculate the footprint score

                    Forward2  =  sum(newvec_f[:shoulder_size])
                    cForward2  = sum(newvec_f[shoulder_size:shoulder_size+fp_size])
                    cBackward2 = sum(newvec_r[shoulder_size:shoulder_size+fp_size])

                    cCombined2 =  sum(newvec_c[shoulder_size:shoulder_size+fp_size])

                    Backward2 =  sum(newvec_r[shoulder_size+fp_size:shoulder_size+shoulder_size+fp_size])
                    #Calculate the FOS score
                    #t_score2 = (((cForward2+1.0)/(fp_size *1.0)) / (Forward2+1.0/(shoulder_size *1.0))) + ((cBackward2+1.0/(fp_size *1.0))/(Backward2+1.0/(shoulder_size *1.0)))
                    t_score2 = (Forward2 + Backward2) / float((Forward2 + Backward2 + cCombined2))
                    scores.append(t_score2)
                if plotting:
                    plt.clf()
                    plt.axvline(t_score)
                    plt.hist(scores)
                    plt.title("Score : {}".format(pct(scores,t_score)))
                    plt.show()
                score = pct(scores,t_score)
                if score < fp_scores[offset]:
                    fp_scores[offset] = score
                    fp_mles[mle] = mle


    print fp_scores
    for pos, value in enumerate(fp_scores):
        if value == 101:
            fp_scores[pos] = 0

    return fp_scores, fp_mles


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
