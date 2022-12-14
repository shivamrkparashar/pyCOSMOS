cimport cython
import numpy as np
from cython.parallel import prange
from libc.math cimport sqrt, fabs, cos
import config

cdef double Lx = config.Lx
cdef double Ly = config.Ly
cdef double Lz = config.Lz
cdef double alpha = config.alpha
cdef double beta = config.beta
cdef double gamma = config.gamma

@cython.boundscheck(False)  # deactivate bounds checking
@cython.wraparound(False)   # deactivate negative indexing.
cdef double distance(double[:] v1, double[:] v2) nogil:
    """
    Calculated distance between two vectors v1 and v2.
    
    Parameters
    ------------------------
    v1: numpy array
    v2: numpy array
    
    Returns
    ------------------------
    distance: float 
    """
    #v1, v2 = v

    cdef double dx, dy, dz
    cdef double squared_distance, dist

    dx = v1[0] - v2[0]
    dy = v1[1] - v2[1]
    dz = v1[2] - v2[2]
    
    if fabs(dx) > Lx/2:
        dx = Lx - fabs(dx)
    
    if fabs(dy) > Ly/2:
        dy = Ly - fabs(dy)
        
    if fabs(dz) > Lz/2:
        dz = Lz - fabs(dz)        
    
    #squared_distance = dx*dx + dy*dy + dz*dz

    squared_distance = dx*dx + dy*dy + dz*dz + 2*dx*dy*cos(gamma) \
                       + 2*dy*dz*cos(alpha) + 2*dz*dx*cos(beta)

    dist = sqrt(squared_distance)

    return dist

@cython.boundscheck(False)  # deactivate bounds checking
@cython.wraparound(False)   # deactivate negative indexing.
def distance_matrix_periodic(double [:]x_arr, double [:] y_arr, double [:] z_arr):
    
    cdef int Npoints = len(x_arr)

    periodic_distance_matrix = np.zeros((Npoints, Npoints))
    cdef double [:, :] periodic_distance_matrix_view = periodic_distance_matrix


    # Parallizing the computation using multiprocessing ironically slows down the speed
    cdef Py_ssize_t i, j
    cdef double [:] vi = np.zeros(3)
    cdef double [:] vj = np.zeros(3)


    #for i in range(Npoints):
    for i in prange(Npoints, nogil =True):
        vi[0] = x_arr[i]
        vi[1] = y_arr[i]
        vi[2] = z_arr[i]


        for j in range(Npoints):
            if (i >= j):
                vj[0] = x_arr[j]
                vj[1] = y_arr[j]
                vj[2] = z_arr[j]
                periodic_distance_matrix_view[i, j] = distance(vi, vj) 

    for i in prange(Npoints, nogil =True):
        for j in range(Npoints):
            if (i < j):
                periodic_distance_matrix_view[i, j] = periodic_distance_matrix_view[j, i]
 
    return periodic_distance_matrix
