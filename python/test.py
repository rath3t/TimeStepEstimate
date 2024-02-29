
import sys
sys.path.append('/build-Debug/python')
import mpi4py
import dune.common
import timestepestimate
import numpy

import pydoc
help_result_string = pydoc.render_doc(timestepestimate)

print(help_result_string)

mat = timestepestimate.Matrix(3, 4)

print(pydoc.render_doc(mat))

print(mat[0,0])
mat[0,0]=7
print(mat[0,0])
n=4
referenceMidSurface= numpy.ndarray(shape=(3,4))
referenceDirectors= numpy.ndarray(shape=(3,4))
displacements= numpy.ndarray(shape=(3,4))
directors= numpy.ndarray(shape=(3,4))


timestepestimate.computeStiffnessMatrix(referenceMidSurface,referenceDirectors,displacements,directors)