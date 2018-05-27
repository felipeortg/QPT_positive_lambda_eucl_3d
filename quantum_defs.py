#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2018-04-11 14:04:32
# @Author  : Felipe Ortega (felipeortegagama@gmail.com)
# @Title   : Quantum definitions

import operator
from config import *

# Complex variable
I = np.complex(0,1)

def heaviside(x):
    if x<0:
        return 0
    else:
        return 1

# Product function
def prod(table):
    return reduce(operator.mul, table)

# Memoization decorators from:
# http://code.activestate.com/recipes/578231-probably-the-fastest-memoization-decorator-in-the-/

def memoize_sing(f):
    """ Memoization decorator for a function taking a single argument """
    class memodict(dict):
        def __missing__(self, key):
            ret = self[key] = f(key)
            return ret 
    return memodict().__getitem__


def memoize_mult(f):
    """ Memoization decorator for a function taking one or more arguments. """
    class memodict1(dict):
        def __getitem__(self, *key):
            return dict.__getitem__(self, key)

        def __missing__(self, key):
            ret = self[key] = f(*key)
            return ret

    return memodict1().__getitem__


# Quantum definitions: 
# -----------------
# Note that we use labels (a,b,c,d,...) multiplied by two to make them integer
# Then the functions like vi, tria, F-sym need extra factors of 1/2
# -----------------
# The definitions are memoized to increase calculation speed.

# Define coupling coefficient 

@memoize_mult
def Delt(a, b, c, K):
    return (((a+b+c)%2==0)
     * heaviside(a + b - c) 
     * heaviside(b + c - a)
     * heaviside(c + a - b)
     * heaviside(2*K - (a+b+c)))

# Define q-numbers
# Technically q(0,K)=1 should be added, but since this is actually never used, it doesn't matter.
@memoize_mult
def q(n, K):
    return ((np.exp(np.pi*n*I / (K+2)) - np.exp(-np.pi*n*I / (K+2)))
    /(np.exp(np.pi*I / (K+2)) - np.exp(-np.pi*I / (K+2))))

 
# Define q-factorial 
# This one actually works for qf(0,K)=1 because the product function gives 1 
# when the maximal value for the iterator is less than the minimum).
@memoize_mult
def qf(n, K):
    if n!=0:
        # Remember: xrange needs to receive an integer type
        return prod([q(nn, K) for nn in xrange(1,n+1)])
    else:
        return 1

 
# Define square root of signed quantum dimension (note that it should be 0 if a>K).
@memoize_mult
def vi(a, K):
    return np.sqrt(np.complex(-1,0))**a * np.sqrt(q(a+1, K)) * heaviside(K-a)

 
# Define total quantum dimension.
# Maybe could be redefine to simply sqrt(sum(vi**4))
@memoize_sing
def DT(K) :
    if K == 0:
        return 1.0
    if onlyintegers:
        # Total qdimension with only integers
        if K == 3 or K == 2:
            return np.sqrt(sum([q(2*nn+1,K)**2 for nn in xrange(2)]))
        else:
            raise ValueError('Only defined for K=2,3 case')

    else:
        return np.sqrt((K + 2)/2.) * 1/np.sin(np.pi / (K+2))

def DT_int(K):
    if K == 3 or K == 2:
        return np.sqrt(sum([q(2*nn+1,K)**2 for nn in xrange(2)]))
    else:
        raise ValueError('Only defined for K=2,3 case')

 # Define triangle function.

@memoize_mult
def tria(a, b, c, K):
    return (Delt(a, b, c, K)
    * np.sqrt(
        qf((a+b-c)/2, K) * qf((a-b+c)/2, K) * qf((-a+b+c)/2, K)
        / qf((a+b+c)/2+1, K)
        )) # Leave as integers since they are admissible spins


# Define F-symbol (memory maximum M~(10^7)elements, K<9).
@memoize_mult
def FS(a, b, e, d, c, f, K):
    if (a <= K and b <= K and e <= K 
    and d <= K and c <= K and f <= K # FS defined only for allowed spins
    and Delt(a, b, e, K) == 1        # Reality conditions
    and Delt(c, d, e, K) == 1 
    and Delt(a, c, f, K) == 1 
    and Delt(b, d, f, K) == 1):      # FS calculation
        return ((-1)**((a+b+c+d)/2.) 
        * np.sqrt(q(e+1, K) * q(f+1, K))
        * tria(a, b, e, K) * tria(c, d, e, K)
        * tria(a, c, f, K) * tria(b, d, f, K)
        * sum(
            [(-1)**nn * qf(nn+1, K) 
            * (qf((a+b+c+d)/2. - nn, K)
                * qf((a+d+e+f)/2. - nn, K)
                * qf((b+c+e+f)/2. - nn, K)
                * qf(nn - (a+b+e)/2., K) * qf(nn - (a+c+f)/2., K) 
                * qf(nn - (b+d+f)/2., K) * qf(nn - (d+c+e)/2., K)
                )**(-1)
            for nn in xrange(
                max(a + b + e, a + c + f, b + d + f, d + c + e)/2,
                min(a + b + c + d, a + d + e + f, b + c + e + f)/2+1
                )]# Leave as integers since they are admissible spins
            ))
    else:
        return 0
