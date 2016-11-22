#!/usr/bin/env ipython3

import numpy as np

def do_sim():
    a = np.random.rand(10,1)
    b = np.random.rand(100,1)
    c = (a[1]+b[1]) * 1000
    return { 'a': a, 'b': b, 'c':c }
