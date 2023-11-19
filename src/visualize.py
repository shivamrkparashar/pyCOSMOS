import config
import numpy as np
import pandas as pd
import mayavi.mlab as mlab
from clustering_functions import dbscan, update_pore_type_matrix, abc_to_xyz



def draw_unit_cell_mayavi():
    """
    draws a unit cell box in mayavi
    :return:
    """
    black = (0, 0, 0)
    white = (1, 1, 1)
    mlab.figure(bgcolor=white, size=(800, 800))

    def draw_line(A, B):
        """
        """
        mlab.plot3d([A[0], B[0]], [A[1], B[1]], [A[2], B[2]], color=black, tube_radius=None)

    va = [0.00, config.Lx]
    vb = [0.00, config.Ly]
    vc = [0.00, config.Lz]

    cord = []  # list of xyz coordinates of vertex of unit cell
    for a in va:
        for b in vb:
            for c in vc:
                cord.append(abc_to_xyz(a, b, c))

    # indexes of the coord list which needs to be connected to form the unit cell
    edges_index = [[0, 1], [0, 2], [0, 4], [5, 4], [5, 1], [5, 7], \
                   [3, 1], [3, 2], [3, 7], [6, 4], [6, 7], [6, 2]]

    for index in edges_index:
        i, j = index
        draw_line(cord[i], cord[j])

def plot_pore_centers_mayavi():
    def plot_sphere(x0, y0, z0, r, color):
        """
        Draws a sphere with (x0, y0, z0) as center and radius = r
        :param x0: x coordinate of the sphere center
        :param y0: y coordinate of the sphere center
        :param z0: z coordinate of the sphere center
        :param r: radius of the sphere
        :param color: color in rgb format. For example (1,0,0) is red
        :return: mayavi mesh
        """
        # phi, theta = (0, 2pi), (0, pi)
        [phi, theta] = np.mgrid[0:2 * np.pi:30j, 0:np.pi:30j]
        x = np.cos(phi) * np.sin(theta)
        y = np.sin(phi) * np.sin(theta)
        z = np.cos(theta)
        return mlab.mesh(r * x + x0, r * y + y0, r * z + z0, color=color)

    def plot_cylinder(x1, y1, z1, x2, y2, z2, r, color):
        """

        draws a cylinder from (x1, y1, z1) to (x2, y2, z2)
        :param r: radius of the cylinder
        :param color: color of the cylinder
        :return:
        """
        return mlab.plot3d([x1, x2], [y1, y2], [z1, z2], tube_radius=r, tube_sides=30, color=color)

    all_cluster_pore_type_labels_flatten = [j for sub in config.all_cluster_pore_type_labels for j in sub]
    all_cluster_center_list_flatten = [j for sub in config.all_cluster_center_list for j in sub]
    all_cluster_diameter_list_flatten = [j for sub in config.all_cluster_diameter_list for j in sub]
    all_cluster_length_list_flatten = [j for sub in config.all_cluster_length_list for j in sub]
    all_cluster_shape_list_flatten = [j for sub in config.all_cluster_shape_list for j in sub]
    all_cluster_orientation_list_flatten = [j for sub in config.all_cluster_orientation_list for j in sub]

    # calculate cluster center in x y z coordinates
    ac, bc, cc = [], [], []
    for cluster_center in all_cluster_center_list_flatten:
        ac.append(cluster_center[0])
        bc.append(cluster_center[1])
        cc.append(cluster_center[2])

    ac, bc, cc = np.array(ac), np.array(bc), np.array(cc)
    xc, yc, zc = abc_to_xyz(ac, bc, cc)

    # red, blue, green, yellow
    shape_color = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0)]

    draw_unit_cell_mayavi()

    for i, center_cord in enumerate(all_cluster_center_list_flatten):

        color = shape_color[all_cluster_pore_type_labels_flatten[i]-1]
        if all_cluster_shape_list_flatten[i] == 'sphere':

            plot_sphere(xc[i], yc[i], zc[i], all_cluster_diameter_list_flatten[i]/2,
                        color=color)

        if all_cluster_shape_list_flatten[i] == 'channel':
            # v1 and v2 are vectors on the end of the cylinder
            v1 = center_cord - all_cluster_orientation_list_flatten[i] * all_cluster_length_list_flatten[i] / 2.
            v2 = center_cord + all_cluster_orientation_list_flatten[i] * all_cluster_length_list_flatten[i] / 2.

            plot_cylinder(v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], all_cluster_diameter_list_flatten[i]/2, color=color)

    #for xi, yi, zi, di, label in zip(xc, yc, zc, all_cluster_diameter_list_flatten,
    #                                 all_cluster_pore_type_labels_flatten):
    #    plot_sphere(xi, yi, zi, di / 2, color=shape_color[label - 1])

    mlab.savefig('fig_pore_centers.jpg')
    # mlab.savefig('pore_centers.obj')
