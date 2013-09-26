cdef extern from "fastbinom.h":
    double bdtrc(int, int, float)
    double bdtr(int, int, float)

def logsf(a,b,c):
    return bdtrc(a, b, c)
def logcdf(a,b,c):
    return bdtr(a,b,c)