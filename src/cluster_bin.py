import config
import numpy as np
import pandas as pd
# from calculate_dbscan_parameters import calculate_dbscan_parameters
from clustering_functions import dbscan, update_pore_type_matrix, abc_to_xyz
import plotly.express as px
import mayavi.mlab as mlab
from identify_cluster_shape import identify_cluster_shape, classify_cluster_shape
import plotly.graph_objs as go
from periodic_distance import distance_matrix_periodic as cython_periodic
import sys

sys.path.append('/home/shivam/Mypymodules')

# element_list = ['Ac', 'Ag', 'Am', 'As', 'At', 'Au', 'B', 'Ba', 'Be', 'X']*30
element_list = range(0, 100, 1)


def extract_points_in_a_bin(x, y, z, geom_dia, bin_address_of_points, boi):
    """

    For a given bin boi,
    Extracts all points within the bin and store in xk, yk, zk

    :param x:
    :param y:
    :param z:
    :param geom_dia:
    :param bin_address_of_points:
    :param boi:
    :return:
    """
    # xk, yk, zk are the points belonging to the bin boi
    xk, yk, zk, diak = [], [], [], []
    for i, index in enumerate(bin_address_of_points):
        # Extract points within a bin with bin_index = boi
        if index == boi:
            xk.append(x[i])
            yk.append(y[i])
            zk.append(z[i])
            diak.append(geom_dia[i])

    print('\n================================================================')
    print('Currently at nbin = %d\n' % (boi))
    # If a bin contains no points
    if len(xk) == 0:
        print("%d does not contain any points" % boi)
    else:
        print("Average geometric diameter of bin: %1.1f Ang" % np.average(diak))
        # return 0

    xk, yk, zk = np.array(xk), np.array(yk), np.array(zk)
    return xk, yk, zk, diak


def cluster_points_within_a_bin(ak, bk, ck, bin_number):
    """
    For a given bin bin_number,
    1. Extracts all points within the bin and store in xk, yk, zk
    2. Use DBSCAN to cluster these points and calculate number of clusters, cluster centers, and size

    :param ak: numpy array, len = No. of geometric points
    :param bk: numpy array, len = No. of geometric points
    :param ck: numpy array, len = No. of geometric points
    :param bin_number: bin number (0, 1, ...20)
    """

    # Store distances between every point and every other points within a bin in pbc_matrix variable
    # Minimum image convention is used
    pbc_matrix = cython_periodic(ak, bk, ck)

    # Using DBSCAN to cluster points
    # If the value of epsilon chosen is too small then a higher number of clusters will be created,
    # and more data points will be taken as noise.
    #
    # min_samples: The number of samples (or total weight) in a neighborhood for a point to be considered as a core point.
    # This includes the point itself.

    # default values 3, and 10 respectively
    e, m = config.eps, config.Nmin
    # e, m = calculate_dbscan_parameters(np.average(diak))

    print("DBSCAN parameters, %1.3f, %1.1f \n" % (e, m))

    # Use DBSCAN
    cluster_labels, cluster_center_list = dbscan(ak, bk, ck, pbc_matrix, eps=e, min_samples=m)

    # Number of clusters and their centers
    Ncluster = max(cluster_labels) + 1
    print('Number of clusters identified = %d ' % Ncluster)

    return cluster_labels, cluster_center_list


def calculate_shape_and_size_of_cluster_within_bin(cluster_labels, bin_number, ak, bk, ck, cluster_center_list):
    Ncluster = max(cluster_labels) + 1
    # Determine the shape of each cluster
    cluster_size_list = []
    cluster_shape_list = []
    cluster_length_list = []
    cluster_orientation_list = []

    for cluster_number in range(Ncluster):
        mask = cluster_labels == cluster_number
        shape, size, length, orientation = classify_cluster_shape(bin_number, cluster_number, ak[mask], bk[mask],
                                                                  ck[mask], cluster_center_list[cluster_number])

        cluster_shape_list.append(shape)
        cluster_size_list.append(size)
        cluster_length_list.append(length)
        cluster_orientation_list.append(orientation)

    return cluster_shape_list, cluster_size_list, cluster_length_list, cluster_orientation_list


def classify_bin(cluster_shape_list, diak, cluster_labels, boi):
    """
    Determines whether a bin is primary or secondary
    :return:
    """

    Ncluster = max(cluster_labels) + 1

    fraction_of_noisy_points = np.count_nonzero(cluster_labels == -1) / len(cluster_labels)
    if fraction_of_noisy_points > 0.2:
        # This bin cannot be regarded as a new pore type
        # print('bin contains noisy points')
        # print('Starting with another bin')
        print(f"fraction of noisy points classified by DBSCAN: {fraction_of_noisy_points}")
        print('nbin = %d is a Secondary bin because fraction of noisy points is larger' % boi)
        return 0

    # calculate number of points in each cluster
    if Ncluster != 0:
        Npoints_per_cluster = (1 - fraction_of_noisy_points) * len(cluster_labels) / Ncluster
    else:
        Npoints_per_cluster = 0

    print("Number of points per cluster = %1.0f" % Npoints_per_cluster)

    if cluster_shape_list[0] == 'sphere':
        Nc = config.Nprobe / config.Volume_of_uc * np.pi / 6 * np.average(diak) ** 3
        print("Points per cluster required for a spherical primary bin = %1.0f " % Nc)

        if Npoints_per_cluster < 0.1 * Nc:
            print('nbin = %d is a Secondary bin because of too many small noisy clusters' % boi)
            return 0

        else:
            return 1

    if cluster_shape_list[0] == 'channel':
        # New pore type found
        return 1

def plot_xyz(xk, yk, zk, cluster_labels):
    """
    plot geometric probe particles using plotly express
    :param xk: numpy arrays
    :param yk: numpy arrays
    :param zk: numpy arrays
    :param cluster_labels: numpy arrays
    :return:
    """
    d = {'x': xk, 'y': yk, 'z': zk, 'color': cluster_labels}
    df = pd.DataFrame(data=d)
    df["color"] = df["color"].astype(str)
    fig = px.scatter_3d(df, x='x', y='y', z='z', color="color")
    fig.update_traces(marker=dict(size=4,
                                  line=dict(width=0.8,
                                            color='Black')),
                      selector=dict(mode='markers'))

    fig.write_html("geometric_points_with_cluster_labels_for_pore_type_%d.html" % config.pore_type_count)
    fig.show()


def save_as_xyz(boi, x_arr, y_arr, z_arr, cluster_labels):
    """
    Saves the x, y, z arrays as an xyz file for visualization in vmd
    A different element is chosen for each clusters

    :param boi: bin number
    :param x_arr:
    :param y_arr:
    :param z_arr:
    :return:
    """
    xyz_file = 'bin_%d.xyz' % boi
    Natoms = len(x_arr)
    # element_list = ['Ar'] * Natoms

    with open(xyz_file, 'w') as out:
        out.write('%d \n' % Natoms)
        out.write('bin number = %d \n' % boi)
        for ai, xi, yi, zi in zip(cluster_labels, x_arr, y_arr, z_arr):
            out.write("%s %1.3f %1.3f %1.3f\n" % (element_list[ai], xi, yi, zi))
