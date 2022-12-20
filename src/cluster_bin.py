import config
import numpy as np
import pandas as pd
#from calculate_dbscan_parameters import calculate_dbscan_parameters
from clustering_functions import dbscan, update_pore_type_matrix, abc_to_xyz
import plotly.express as px
from periodic_distance import distance_matrix_periodic as cython_periodic
import sys
sys.path.append('/home/shivam/Mypymodules')

element_list = ['Ac', 'Ag', 'Am', 'As', 'At', 'Au', 'B', 'Ba', 'Be', 'X']*10


def extract_points_in_a_bin(x, y, z, geom_dia, bin_index_of_points, boi):
    """

    For a given bin boi,
    Extracts all points within the bin and store in xk, yk, zk

    :param x:
    :param y:
    :param z:
    :param geom_dia:
    :param bin_index_of_points:
    :param boi:
    :return:
    """
    # xk, yk, zk are the points belonging to the bin boi
    xk, yk, zk, diak = [], [], [], []
    for i, index in enumerate(bin_index_of_points):
        # Extract points within a bin with bin_index = boi
        if index == boi:
            xk.append(x[i])
            yk.append(y[i])
            zk.append(z[i])
            diak.append(geom_dia[i])

    print('\n Currently at nbin = %d, & D = %1.1f Ang' % (boi, np.average(diak)))
    # If a bin contains no points
    if len(xk) == 0:
        print("%d does not contain any points" % boi)
        # return 0

    xk, yk, zk = np.array(xk), np.array(yk), np.array(zk)
    return xk, yk, zk, diak


def cluster_points_within_a_bin(xk, yk, zk):
    """
    For a given bin boi,
    1. Extracts all points within the bin and store in xk, yk, zk
    2. Use DBSCAN to cluster these points and calculate number of clusters, cluster centers, and size

    :param x: numpy array, len = No. of geometric points
    :param y: numpy array, len = No. of geometric points
    :param z: numpy array, len = No. of geometric points
    :param bin_index_of_points: numpy array of len = No. of geometric points
    Tells in which bin does each of the point goes
    :param boi: bin number (0, 1, ...20)
    """

    # Store distances between every point and every other points within a bin in pbc_matrix variable
    # Minimum image convention is used
    pbc_matrix = cython_periodic(xk, yk, zk)

    # Using DBSCAN to cluster points
    # If the value of epsilon chosen is too small then a higher number of clusters will be created,
    # and more data points will be taken as noise.
    #
    # min_samples: The number of samples (or total weight) in a neighborhood for a point to be considered as a core point.
    # This includes the point itself.

    # e and m are different for each bin and each mof
    # default values 3, and 10 respectively
    #e, m = 2., 40
    e, m = 2.5, 25
    #e, m = calculate_dbscan_parameters(np.average(diak))

    print("DBSCAN parameters, %1.3f, %1.1f " %(e, m))

    cluster_labels, cluster_center_list, cluster_diameter_list, fraction_of_noisy_points = \
        dbscan(xk, yk, zk, pbc_matrix, eps=e, min_samples=m)

    # Number of clusters and their centers
    #Ncluster = max(np.unique(cluster_labels)) + 1
    Ncluster = max(cluster_labels) + 1
    print('Number of clusters identified = %d ' % Ncluster)

    # calculate number of points in each cluster
    if Ncluster != 0:
        Npoints_per_cluster = (1-fraction_of_noisy_points)*len(xk)/Ncluster
    else:
        Npoints_per_cluster = 0

    return Ncluster, Npoints_per_cluster, cluster_labels, cluster_center_list, cluster_diameter_list, fraction_of_noisy_points


def classify_bin(diak, Npoints_per_cluster, fraction_of_noisy_points, boi):
    """
    Determines whether a bin is primary or secondary


    :return:
    """

    print("Number of points per cluster = " , Npoints_per_cluster)
    Nc = 50000/config.Volume_of_uc *np.pi/6* np.average(diak) **3
    print("Points per cluster required for a spherical primary bin = " , Nc)

    if Npoints_per_cluster < 0.5*Nc:
        print('nbin = %d is a Secondary bin because of too many small noisy clusters' % boi)
        return 0

    if Npoints_per_cluster > Nc:
        # print('nbin = %d is a Secondary bin because of big noisy clusters' % boi)
        print('Reduce the DBSCAN parameter epsilon' % boi)
        return 0

    if fraction_of_noisy_points > 0.2:
        # This bin cannot be regarded as a new pore type
        # print('bin contains noisy points')
        # print('Starting with another bin')
        print('nbin = %d is a Secondary bin because fraction of noisy points is larger' % boi)
        return 0
    else:
        # New pore type found
        return 1


def cluster_bin(xk, yk, zk, diak, boi):
    """
    For a given bin boi,
    1. Extracts all points within the bin and store in xk, yk, zk
    2. Use DBSCAN to cluster these points and calculate number of clusters, cluster centers, and size
    3. Displays all clusters of points in html file on browser
    4. Update the pore type matrix with the new pore type label

    :param x: numpy array, len = No. of geometric points
    :param y: numpy array, len = No. of geometric points
    :param z: numpy array, len = No. of geometric points
    :param bin_index_of_points: numpy array of len = No. of geometric points
    Tells in which bin does each of the point goes
    :param boi: bin number (0, 1, ...20)

    :returns:
     0 if the bin boi is a secondary bin
     1 if the bin boi is a primary bin

    """



    # Store distances between every point and every other points within a bin in pbc_matrix variable
    # Minimum image convention is used
    xk, yk, zk = np.array(xk), np.array(yk), np.array(zk)
    pbc_matrix = cython_periodic(xk, yk, zk)

    # Using DBSCAN to cluster points
    # If the value of epsilon chosen is too small then a higher number of clusters will be created,
    # and more data points will be taken as noise.
    #
    # min_samples: The number of samples (or total weight) in a neighborhood for a point to be considered as a core point.
    # This includes the point itself.

    # e and m are different for each bin and each mof
    # default values 3, and 10 respectively
    #e, m = 2., 40
    e, m = 2.5, 25
    #e, m = calculate_dbscan_parameters(np.average(diak))

    print("DBSCAN parameters, %1.3f, %1.1f " %(e, m))

    cluster_labels, cluster_center_list, cluster_diameter_list, fraction_of_noisy_points = \
        dbscan(xk, yk, zk, pbc_matrix, eps=e, min_samples=m)

    # Number of clusters and their centers
    #Ncluster = max(np.unique(cluster_labels)) + 1
    Ncluster = max(cluster_labels) + 1
    print('Number of clusters identified = %d ' % Ncluster)

    # calculate number of points in each cluster
    if Ncluster != 0:
        Npoints_per_cluster = (1-fraction_of_noisy_points)*len(xk)/Ncluster
    else:
        Npoints_per_cluster = 0

    print("Number of points per cluster = " , Npoints_per_cluster)
    Nc = 50000/config.Volume_of_uc *np.pi/6* np.average(diak) **3
    print("Points per cluster required for a spherical primary bin = " , Nc)
    # compare it with


    """
    if Npoints_per_cluster < 0.5*Nc:
        print('nbin = %d is a Secondary bin because of too many small noisy clusters' % boi)
        return 0
    """

    if Npoints_per_cluster >  Nc:
        #print('nbin = %d is a Secondary bin because of big noisy clusters' % boi)
        print('Reduce the DBSCAN parameter epsilon' % boi)
        return 0

    if fraction_of_noisy_points > 0.2:
        # This bin cannot be regarded as a new pore type
        # print('bin contains noisy points')
        # print('Starting with another bin')
        print('nbin = %d is a Secondary bin because fraction of noisy points is larger' % boi)
        return 0
    else:
        # New pore type found
        config.pore_type_count += 1
        print('New pore type identified, pore type count = ', config.pore_type_count)
        config.all_cluster_center_list.append(cluster_center_list)
        config.all_cluster_diameter_list.append(cluster_diameter_list)
        config.all_cluster_pore_type_labels.append(Ncluster * [config.pore_type_count])

    # What about the points that were classified as noise by DBSCAN?
    # We are just ignoring them for the time being
    update_pore_type_matrix(xk, yk, zk, cluster_labels)


    # plot x, y, z coordinates in 3d
    xxk, yyk, zzk = abc_to_xyz(xk, yk, zk)

    save_as_xyz(boi, xxk, yyk, zzk, cluster_labels)

    d = {'x': xxk, 'y': yyk, 'z': zzk, 'color': cluster_labels}
    df = pd.DataFrame(data=d)
    df["color"] = df["color"].astype(str)
    fig = px.scatter_3d(df, x='x', y='y', z='z', color="color")
    fig.update_traces(marker=dict(size=4,
                                  line=dict(width=0.8,
                                            color='Black')),
                      selector=dict(mode='markers'))

    fig.write_html("geometric_points_with_cluster_labels_for_pore_type_%d.html" %config.pore_type_count)
    fig.show()


    print('nbin = %d is a Primary bin' % boi)
    return 1

def plot_xyz(xk, yk, zk, cluster_labels):
    """

    :param xk:
    :param yk:
    :param zk:
    :param cluster_labels:
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

    fig.write_html("geometric_points_with_cluster_labels_for_pore_type_%d.html" %config.pore_type_count)
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
    xyz_file = 'bin_%d.xyz' %boi
    Natoms = len(x_arr)
    #element_list = ['Ar'] * Natoms

    with open(xyz_file, 'w') as out:
        out.write('%d \n' %Natoms)
        out.write('bin number = %d \n' %boi)
        for ai, xi, yi, zi in zip(cluster_labels, x_arr, y_arr, z_arr):
            out.write("%s %1.3f %1.3f %1.3f\n" %(element_list[ai], xi, yi, zi))

