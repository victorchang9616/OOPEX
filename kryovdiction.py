
import numpy as np
from scipy.linalg import solve_triangular
import csv
from joblib import Parallel, delayed
import time
import math
def add(a,b):
    c = {}
    nz_idx = set(a.keys()).union(set(b.keys())) 
    for i in nz_idx:
        c[i] = 0
        if i in a:
            c[i] += a[i]
        if i in b:
            c[i] += b[i]
    return c

def innerprod(a,b):
    c = 0.0
    nz_idx = set(a.keys()).intersection(set(b.keys())) 
    for  i in nz_idx:
        c += a[i]*b[i]
    return c
def norm(a):
    sqa=math.sqrt(innerprod(a,a))
    return sqa
def cleandict(a):
    result={}
    for x, y in a.items(): 
        if(abs(y) > 1e-7):result.update({x:y})
    return result