#!/usr/bin/env ipython3

from ipyparallel import Client

# See "ipyparallel" directory for sample profile
rc = Client(profile='mpi')

# `lv` is the load balanced view, will distribute tasks according
# to availability rather than round-robin (direct view)
lv = rc.load_balanced_view()
lv.block = True

# Make the numpy library available everywhere
with rc[:].sync_imports():
    import numpy

# The `parallel` decorator adds mojo to make this function
# parallel on our cluster
@lv.parallel()
def f(x):
    a = numpy.random.rand(10,1)
    b = numpy.random.rand(100,1)
    c = ( a[0] + b[0] ) * 1000
    return a,b,c

# `result` will be an array of tuples roughly consistent
# with the structures in the list in the original matlab
# code:
#
#    PYTHON                 MATLAB
#   -------------------    ---------------------
#   struct_array1[0][0] ~= struct_array1[1].a 
#   struct_array1[1][0] ~= struct_array1[2].a 
#   ...
#   struct_array1[n][0] ~= struct_array1[n+1].a 
#
# It would be possible to instead return an object, but adds
# lots of code without clearly adding value.  It is likely
# more valuable to have struct_array1 as an ndarray or matrix
# than object.

struct_array1 = f.map(range(100))
