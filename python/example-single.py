#!/usr/bin/evn python3

import numpy as np

# struct1=struct('a',zeros(10,1),'b',zeros(100,1),'c',uint32(0));
# struct_array1(1:100,1)=struct1;
#
# In matlab, `zeros` creates a matrix.  In numpy, matrices are primarily used
# for facilitating linear algebra computations.  As this are one dimentional
# matrices, the ndarray will work just fine

# `struct_array1` is a simple list of dictionary objects, each of which contains
# three elements, two ndarrays and one long value
struct_array1 = []

#parfor i=1:100
for i in range(100):
#    a=rand(10,1);
    a = np.random.rand(10,1)
#    b=rand(100,1);
    b = np.random.rand(100,1)
#    c=uint32((a(1)+b(1))*1000);
    c = (a[0]+b[0]) * 1000
#    struct_array1(i).a=a;
#    struct_array1(i).b=b;
#    struct_array1(i).c=c;
    struct_array1.append(
        {
            'a': a,
            'b': b,
            'c':c
        }
    )
