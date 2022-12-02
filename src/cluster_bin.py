import config
import numpy as np
import pandas as pd
from clustering_functions import dbscan, update_pore_type_matrix
import plotly.express as px
from periodic_distance import distance_matrix_periodic as cython_periodic
import sys
sys.path.append('/home/shivam/Mypymodules')

element_list = ['Ac', 'Ag', 'Am', 'As', 'At', 'Au', 'B', 'Ba', 'Be', 'X']*10
def cluster_bin(x, y, z, bin_index_of_points, boi):
    """

    :param x:
    :param y:
    :param z:
    :param bin_index_of_points:
    :param boi:
    :param compute_pbc_matrix:
    :return:
    """
    print('Currently at nbin = %d' % boi)

    # xk, yk, zk are the points belonging to the bin boi
    xk, yk, zk, diak = [], [], [], []
    for i, index in enumerate(bin_index_of_points):
        # Extract points within a bin with bin_index = boi
        if index == boi:
            xk.append(x[i])
            yk.append(y[i])
            zk.append(z[i])
            #diak.append(diameter[i])

    # Store distances between every point and every other points within a bin in pbc_matrix variable
    # Minimum image convention is used
    xk, yk, zk = np.array(xk), np.array(yk), np.array(zk)
    pbc_matrix = cython_periodic(xk, yk, zk)

    # Using DBSCAN to cluster points
    cluster_labels, cluster_center_list, cluster_diameter_list, fraction_of_noisy_points = \
        dbscan(xk, yk, zk, pbc_matrix, eps=3.0, min_samples=10)

    # Number of clusters and their centers
    #Ncluster = max(np.unique(cluster_labels)) + 1
    Ncluster = max(cluster_labels) + 1
    print('Number of clusters identified = %d ' % Ncluster)

    # append the centers of the clusters to the global list
    if not config.all_cluster_center_list:  # if all_cluster_center_list is empty
        print("First pore type identified")
        config.pore_type_count += 1
        config.all_cluster_center_list.append(cluster_center_list)
        config.all_cluster_diameter_list.append(cluster_diameter_list)
        config.all_cluster_pore_type_labels.append(Ncluster * [config.pore_type_count])  # pore type labels of this bin

    else: # For 2nd or 3rd pore type
        # The points within a bin will be a new pore type if the
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


    # plot x, y, z coordinates in 3d
    d = {'x': xk, 'y': yk, 'z': zk, 'color': cluster_labels}
    df = pd.DataFrame(data=d)
    df["color"] = df["color"].astype(str)
    fig = px.scatter_3d(df, x='x', y='y', z='z', color="color")
    fig.update_traces(marker=dict(size=4,
                                  line=dict(width=0.8,
                                            color='Black')),
                      selector=dict(mode='markers'))

    fig.write_html("geometric_points_with_cluster_labels_for_pore_type_%d.html" %config.pore_type_count)
    fig.show()

    # What about the points that were classified as noise by DBSCAN?
    # We are just ignoring them for the time being

    # Add the points to the grid which are not classified as noise
    # And update the pore type matrix
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

