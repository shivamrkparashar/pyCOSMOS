#from math import sqrt

cimport cython
import numpy as np
from cython.parallel import prange
from libc.math cimport sqrt, fabs
from config import Box
#Box = 25.832
#Lx = Box
#Ly = Box
#Lz = Box

#cdef double Box = 25.832
cdef double Lx = Box
cdef double Ly = Box
cdef double Lz = Box

#import multiprocessing


#cdef distance(double[:] v1, double[:] v2) nogil:
#def distance(double[:] v1, double[:] v2):
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
    
    squared_distance = dx*dx + dy*dy + dz*dz
    dist = sqrt(squared_distance)

    return dist

@cython.boundscheck(False)  # deactivate bounds checking
@cython.wraparound(False)   # deactivate negative indexing.
#def distance_matrix_periodic(x_arr, y_arr, z_arr):
def distance_matrix_periodic(double [:]x_arr, double [:] y_arr, double [:] z_arr):
    
    #print(len(x_arr))
    # print('3')
    cdef int Npoints = len(x_arr)
    #cdef int Npoints = x_arr.shape[0]

    periodic_distance_matrix = np.zeros((Npoints, Npoints))
    cdef double [:, :] periodic_distance_matrix_view = periodic_distance_matrix


    # Parallizing the computation using multiprocessing ironically slows down the speed
    """
    for i in range(Npoints):
        vi = np.array([x_arr[i], y_arr[i], z_arr[i]])
        vj_list = [np.array([x_arr[j], y_arr[j], z_arr[j]]) for j in range(Npoints) if i >= j]

        vi_list = [vi]* len(vj_list)


        a_pool = multiprocessing.Pool()
        input_vector_list = zip(vi_list, vj_list)
        #input_vector_list = [ (vi, np.array([x_arr[j], y_arr[j], z_arr[j]])) for j in range(Npoints) if i >= j]
        
        results = a_pool.map(distance, input_vector_list)

        for j in range(Npoints):
            if i >= j:
                periodic_distance_matrix[i][j] = results[j]
        
        a_pool.close()
    """
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
