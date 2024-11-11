#! /usr/bin/python3
#####################
# Imports
###
import sys
import math
import numpy as np
from scipy.constants import *

#####################
# Constants
###
kb=0.0019872041 #kcal mol-1 K-1

#####################
# Functions
###

# Mean and stddev
def mean_stddev(values):
    s   = 0
    ssq = 0
    n   = float(len(values))
    for i in values:
        s   =   s + i
        ssq = ssq + i*i
    mean   = s/n
    stddev = math.sqrt(ssq/n - mean*mean)
    return mean, stddev

# Read a clean *x.xvg file with two columns
# Output: two vectors corresponding to the two columns
def read_xxvg(f):
    x = []
    y = []
    with open(f) as fi:
        for line in fi:
            if line[0] != "@" and line [0] != "#":
                vals = line.split()
                x.append(float(vals[0]))
                y.append(float(vals[1]))
    return x,y

# Read a PMF file with two columns
# Output: two vectors corresponding to the two columns (only >= 0 values are added)
def read_pmf(f):
    x = []
    y = []
    with open(f) as fi:
        for line in fi:
            if line[0] != "@" and line [0] != "#":
                vals = line.split()
                if float(vals[0]) >= 0:
                    x.append(float(vals[0]))
                    y.append(float(vals[1]))
    return x,y

# The input/output is equal to timecorr (A. M. Baptista)
# corr = Czz from equation 1 in 10.1021/jacs.6b11215
def timecorr(x, y):
    nint = len(x)
    n    = float(len(x))
    nchk = float(len(y))
    if n != nchk:
        sys.exit("(Error in timecorr function) The two arrays must have equal size.")
    #AVXs and SDs
    avx1,sd1 = mean_stddev(x)
    avx2,sd2 = mean_stddev(y)
    #deltaA eq.2.105 and 2.44 Allen and Tildesley
    xmod=[]
    ymod=[]
    for i in range(nint):
        xmod.append(x[i]-avx1)
        ymod.append(y[i]-avx2)
    #Time column
    time=[i for i in range(nint)]
    # non-normalized correlation function
    corr = []
    for i in range(nint):
        s = 0
        aaa = 0
        for j in range(0,nint-i):
            s = s + xmod[j+i]*ymod[j]
        corr.append(s/(n-i))
    # normalization and integration
    normf    = sd1*sd2
    normcorr = []
    tcorr    = []
    tcorrpar = 0
    for i in range(nint):
        normcorr.append(corr[i]/normf)
        tcorrpar = tcorrpar + corr[i]/normf
        tcorr.append(tcorrpar)
    return time,normcorr,corr,tcorr

#####################
# Functions related to equations 1-3 from 10.1021/jacs.6b11215
###
def integral_Czz(corr, var, dt):
    s   = 0
    chk = 0
    for i in range(len(corr)):
        if corr[i]>0.01*var and chk == 0:
            s = s + corr[i]*dt
        elif (corr[i]<0.01*var):
            chk=1
    return s

# The input is in ps and nm
# Output in cm2 s-1
def Dz(times, positions, b, jump):
    dt = times[1] - times[0]
    avx,stddev = mean_stddev(positions)
    var = stddev * stddev
    data = []
    autocorrs = []
    for i in range(len(times)):
        if times[i] > b:
            data.append(positions[i])
            if times[i]%jump == 0:
                avx,stddev = mean_stddev (data)
                var = stddev * stddev
                time, normcorr, corr, tcorr = timecorr(data, data)
                autocorrs.append(integral_Czz(corr, var, dt))
                data=[]
    out_data = mean_stddev(autocorrs)
    return 0.01 * var * var / out_data[0]

# The input is
# dG in kcal mol-1
# beta mol kcal-1
# The output will have the inverse units of Dz: usually cm2 s-1 -> s cm-2
def Rz(pmf, dzs, beta):
    rzs = []
    for i in range(len(pmf)):
        rzs.append(math.exp(beta * pmf[i]) / dzs[i])
    return rzs

# The input is
# dG in kcal mol-1
# beta mol kcal-1
# The output will have the inverse units of Dz: usually cm2 s-1 -> s cm-2
def Peff(UScenters, rzs):
    if len(UScenters) != len(rzs):
        sys.exit("(Error in Peff function) The two arrays must have equal size.")
    s = 0
    for i in range(1,len(UScenters)):
        s = s + ( ( UScenters[i] - UScenters[i-1] ) * 10**-7 ) * ( ( rzs[i] + rzs[i-1] ) / 2 )
    return 1/s
        

#####################
# Main
###
#if __name__ == "__main__":
    ###
    # Inputs
    #
tol        = 0.01 # ps
T          = 310 # K
xxvgs      = "./umb_{}x.xvg"
maxdist    = 18 # Angstrom
pmf        = "./aux-bsResult"
equilib    = 0 # ps
block_size = 10000 # 10000 # ps

beta = kb * T
UScenters, pmf = read_pmf(pmf)

dzs = []
for i in range(maxdist+1):
    val = '{:3.1f}'.format(float(i)/5)
    times, positions = read_xxvg(xxvgs.format(val)) 
    dzs.append(Dz(times, positions, equilib, block_size))
    
rzs = Rz(pmf, dzs, beta)
peff = Peff(UScenters, rzs)
print ('#Peff = {} cm s-1'.format(2*peff))
for i in range(len(pmf)):
    print(UScenters[i], pmf[i], dzs[i], rzs[i])
