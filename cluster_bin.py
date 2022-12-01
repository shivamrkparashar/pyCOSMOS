import config
import time
import numpy as np
from clustering_functions import dbscan, approx_equal, update_pore_type_matrix
from periodic_distance import distance_matrix_periodic as cython_periodic
import sys
sys.path.append('/home/shivam/Mypymodules')
from shivmod import PrinttoFile

element_list = ['Ac', 'Ag', 'Am', 'As', 'At', 'Au', 'B', 'Ba', 'Be', 'X']
def cluster_bin(x, y, z, bin_index_of_points, boi, compute_pbc_matrix):
    print('Currently at nbin = %d' % boi)

    # x, y, z coordinates of points in that bin
    if compute_pbc_matrix:
        xk, yk, zk, diak = [], [], [], []
        for i, index in enumerate(bin_index_of_points):
            # Extract points within a bin with bin_index = boi
            if index == boi:
                xk.append(x[i])
                yk.append(y[i])
                zk.append(z[i])
                #diak.append(diameter[i])

        # Compute distances between each point using periodic boundary conditions
        xk, yk, zk = np.array(xk), np.array(yk), np.array(zk)
        pbc_matrix = cython_periodic(xk, yk, zk)


    else: # If not computer_pbc_matrix, load fom saved numpy arrays
        start = time.time()
        pbc_matrix = np.load('saved_np_array/pbc_matrix_%d.npy' %boi)
        xk = np.load('saved_np_array/x_coordinates_%d.npy' %boi)
        yk = np.load('saved_np_array/y_coordinates_%d.npy' %boi)
        zk = np.load('saved_np_array/z_coordinates_%d.npy' %boi)
        end = time.time()
        print('Time to load np matrices: ', end - start)

    # DBSCAN to cluster points
    cluster_labels, cluster_center_list, cluster_diameter_list, fraction_of_noisy_points = \
        dbscan(xk, yk, zk, pbc_matrix, eps=3.0, min_samples=10)

    # Number of clusters and theri centers
    Ncluster = max(np.unique(cluster_labels)) + 1
    print('Number of clusters identified = %d ' % Ncluster)

    # append the centers of the clusters to the global list
    if not config.all_cluster_center_list:  # if all_cluster_center_list is empty
        print("First pore type identified")
        config.pore_type_count += 1
        config.all_cluster_center_list.append(cluster_center_list)
        config.all_cluster_diameter_list.append(cluster_diameter_list)
        config.all_cluster_pore_type_labels.append(Ncluster * [config.pore_type_count])  # pore type labels of this bin
        # Highest histogram peak always correspond to the largest pore
        #config.Ncluster_list.append(Ncluster)
        #print('List of centers of clusters', config.all_cluster_center_list)
        print('Largest pore centers identified')

    else: # For 1st, 2nd or 3rd pore type
        # Check if the cluster centers are already in the all_cluster_center_list
        # If the recently found cluster center is close to the centers found in the previous bin
        """
        if approx_equal(config.all_cluster_center_list[-2], config.all_cluster_center_list[-1]):
            # Do not change the pore_type_count
            # This pore type was already identified in the previous bin
            None
        """
        if fraction_of_noisy_points > 0.5:
            # This bin cannot be regarded as a new pore type
            print('bin contains noisy points')
            print('Starting with another bin')
            return 0
        else:
            # New pore type found
            config.pore_type_count += 1
            print('New pore type identified, pore type count = ', config.pore_type_count)
            config.all_cluster_center_list.append(cluster_center_list)
            config.all_cluster_diameter_list.append(cluster_diameter_list)
            config.all_cluster_pore_type_labels.append(Ncluster * [config.pore_type_count])

            #config.Ncluster_list.append(Ncluster)
            #print('List of centers of clusters', config.all_cluster_center_list)

    # What about the points that were classified as noise by DBSCAN? We are just ignoring them for the time being
    # Add the points to the grid which are not classified as noise

    xgrid, ygrid, zgrid = [], [], []
    pore_type_labels = []
    for i, li in enumerate(cluster_labels):
        if li != -1:  # if not noise
            pore_type_labels.append(config.pore_type_count)
            xgrid.append(xk[i])
            ygrid.append(yk[i])
            zgrid.append(zk[i])

    update_pore_type_matrix(xgrid, ygrid, zgrid, pore_type_labels)
    save_as_xyz(boi, xk, yk, zk, cluster_labels)

    return 1
    #return xk, yk, zk, pore_type_labelsL


def save_as_xyz(boi, x_arr, y_arr, z_arr, cluster_labels):
    """
    Saves the x, y, z arrays as an xyz file for visualization in vmd
    :param boi: bin number
    :param x_arr:
    :param y_arr:
    :param z_arr:
    :return:
    """
    xyz_file = 'bin_%d.xyz' %boi
    Natoms = len(x_arr)
    #element_list = ['Ar'] * Natoms

    with open(xyz_file, 'w') as out:
        out.write('%d \n' %Natoms)
        out.write('bin number = %d \n' %boi)
        for ai, xi, yi, zi in zip(cluster_labels, x_arr, y_arr, z_arr):
            out.write("%s %1.3f %1.3f %1.3f\n" %(element_list[ai], xi, yi, zi))

