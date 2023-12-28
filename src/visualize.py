import config
import numpy as np
import pandas as pd
from matplotlib import colors
import mayavi.mlab as mlab
import os, glob
from clustering_functions import dbscan, update_pore_type_matrix, abc_to_xyz, put_point_in_box

def draw_unit_cell_mayavi():
    """
    draws a unit cell box in mayavi
    :return:
    """
    black = (0, 0, 0)
    white = (1, 1, 1)
    mlab.figure(bgcolor=white, size=(1024, 1024))

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

def plot_pore_centers_mayavi():

    all_cluster_pore_type_labels_flatten = [j for sub in config.all_cluster_pore_type_labels for j in sub]
    all_cluster_center_list_flatten = [j for sub in config.all_cluster_center_list for j in sub]
    all_cluster_diameter_list_flatten = [j for sub in config.all_cluster_diameter_list for j in sub]
    all_cluster_length_list_flatten = [j for sub in config.all_cluster_length_list for j in sub]
    all_cluster_shape_list_flatten = [j for sub in config.all_cluster_shape_list for j in sub]
    all_cluster_orientation_list_flatten = [j for sub in config.all_cluster_orientation_list for j in sub]

    # convert color to (r,g,b) between 0 and 1
    shape_color_rgb = list(map(colors.to_rgb, config.color_list))

    draw_unit_cell_mayavi()

    for i, center_cord in enumerate(all_cluster_center_list_flatten):
        ac, bc, cc = center_cord
        xc, yc, zc = abc_to_xyz(ac, bc, cc)

        color = shape_color_rgb[all_cluster_pore_type_labels_flatten[i]-1]
        if all_cluster_shape_list_flatten[i] == 'sphere':

            #plot_sphere(xc[i], yc[i], zc[i], all_cluster_diameter_list_flatten[i]/2,
            plot_sphere(xc, yc, zc, all_cluster_diameter_list_flatten[i]/2,
                        color=color)

        if all_cluster_shape_list_flatten[i] == 'channel':
            # v1 and v2 are vectors on the end of the cylinder
            v1 = center_cord - all_cluster_orientation_list_flatten[i] * all_cluster_length_list_flatten[i] / 2.
            v2 = center_cord + all_cluster_orientation_list_flatten[i] * all_cluster_length_list_flatten[i] / 2.

            v1 = np.array(abc_to_xyz(v1[0], v1[1], v1[2]))
            v2 = np.array(abc_to_xyz(v2[0], v2[1], v2[2]))
            plot_cylinder(v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], all_cluster_diameter_list_flatten[i]/2, color=color)

    #for xi, yi, zi, di, label in zip(xc, yc, zc, all_cluster_diameter_list_flatten,
    #                                 all_cluster_pore_type_labels_flatten):
    #    plot_sphere(xi, yi, zi, di / 2, color=shape_color[label - 1])

    mlab.savefig('fig_pore_centers.jpg')

    ciffile = glob.glob('*.cif')[0]
    pdbfile = cif_to_pdb(ciffile)
    visualize_pdb(pdbfile)
    # mlab.savefig('pore_centers.obj')


def cif_to_pdb(ciffile):

    # Convert cif file into pdb using openbabel
    pdbfile = ciffile.replace('.cif', '.pdb')
    if not os.path.exists(pdbfile):
        os.system('obabel -icif %s -opdb -O%s' %(ciffile, pdbfile))

    return pdbfile


def read_pdb(pdbfile):

    name, X, Y, Z = [], [], [], []
    with open(pdbfile, 'r') as inp:
        for line in inp:
            if 'HETATM' in line:
                s = line.split()
                n, x, y, z = s[-1], float(s[5]), float(s[6]), float(s[7])
                name.append(n)
                X.append(x)
                Y.append(y)
                Z.append(z)

    return name, X, Y, Z

def visualize_pdb(pdbfile):
    """
    visualize pdb file in mayavi
    draws atoms as sphere
    bonds are not plotted
    :param pdbfile: path to pdb file
    :return: None
    """
    name, X, Y, Z = read_pdb(pdbfile)
    radius = {'C': 0.8, 'H': 0.5, 'N': 0.8, 'O': 0.8, 'Zn': 1.2, 'Zr': 1.2}
    color = {'C': 'grey', 'H': 'snow', 'N': 'cyan', 'O': 'darkorange', 'Zn': 'purple', 'Zr':
             'purple'}
    color_rgb = {k: colors.to_rgb(v) for k, v in color.items()}

    for ni, x, y, z in zip(name, X, Y, Z):
        if ni not in radius.keys():
            radius[ni] = 1.2
            color_rgb[ni] = colors.to_rgb('purple')
        plot_sphere(x, y, z, radius[ni], color_rgb[ni])

    # set angle and save images
    scene = mlab.gcf()
    scene.scene.parallel_projection = True

    scene.scene.x_minus_view()
    scene.scene.save('img_x_minus.png')

    scene.scene.y_plus_view()
    scene.scene.save('img_y_plus.png')

    scene.scene.z_minus_view()
    scene.scene.save('img_z_minus.png')

    scene.scene.isometric_view()
    scene.scene.save('img_isometric.png')

    
    view_parameter_list = [(0, 30, 50), (0, 45, 50), (30, 0, 50), (45, 0, 50), (30, 30, 50), (30, 45, 50), (45, 45, 50)]
        
    for (azimuth, elevation, distance) in view_parameter_list:
        mlab.view(azimuth=azimuth, elevation=elevation, distance=distance)
        scene.scene.save(f'img_azimuth_{azimuth}_elevation_{elevation}_distance_{distance}.png')
    """
    # custom 1
    scene.scene.camera.position = [52.05004265999461, 54.272623083262744, -124.47365033936694]
    scene.scene.camera.focal_point = [14.314302832764572, 13.977785034102457, 9.314656441821107]
    scene.scene.camera.view_angle = 30.0
    scene.scene.camera.view_up = [-0.07972495641319911, 0.9604526701211306, 0.26678568136634456]
    scene.scene.camera.clipping_range = [75.6048828398952, 232.09492615776236]
    scene.scene.camera.compute_view_plane_normal()
    scene.scene.save('img_custom1.png')

    # custom 2
    scene.scene.camera.position = [-82.59548057900444, 96.79764098478564, 77.84447772048276]
    scene.scene.camera.focal_point = [14.314302832764572, 13.977785034102457, 9.314656441821107]
    scene.scene.camera.view_angle = 30.0
    scene.scene.camera.view_up = [0.48825015022348367, 0.8195292312442335, -0.2999727153308802]
    scene.scene.camera.clipping_range = [71.27508842595759, 237.54524525670973]
    scene.scene.camera.compute_view_plane_normal()
    scene.scene.save('img_custom2.png')

    # custom 3
    scene.scene.camera.position = [65.23359320485599, -59.71152353549938, 124.11622146482044]
    scene.scene.camera.focal_point = [9.481975127823052, 13.720633808393874, 12.55257857562868]
    scene.scene.camera.view_angle = 30.0
    scene.scene.camera.view_up = [0.9219161580594645, 0.17443063037889395, -0.3458967370393918]
    scene.scene.camera.clipping_range = [70.48147970504198, 240.27387800193512]
    scene.scene.camera.compute_view_plane_normal()
    scene.scene.save('img_custom3.png')
    """


