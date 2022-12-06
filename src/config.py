from math import ceil
import sys
import numpy as np

# default parameters
INPUT = 'IRMOF-1.vpsdpts'
Box = 25.832
Npores = 2

INPUT = sys.argv[1]
Box = float(sys.argv[2])
Npores = int(sys.argv[3])
File_pore_type_matrix = sys.argv[4]

print('INPUT is ', INPUT)
print('Number of pore types = ', Npores)

Lx, Ly, Lz = Box, Box, Box
pore_type_matrix = -1 * np.ones((int(ceil(Lx)),int(ceil(Ly)),int(ceil(Lz)))) # 3D pore type matrix
pore_type_matrix_2 = -1 * np.ones((int(ceil(Lx)),int(ceil(Ly)),int(ceil(Lz)))) # 3D pore type matrix
pore_type_matrix_3 = -1 * np.ones((int(ceil(Lx)),int(ceil(Ly)),int(ceil(Lz)))) # 3D pore type matrix
all_cluster_center_list = []  # [ [c1, c2], [c3, c4, c5, c6] ......... ] where ci = np.array([xi, yi, zi])
#Ncluster_list = [] # Number of clusters identified by DBSCAN in each bin
all_cluster_pore_type_labels = [] # [[0, 0,], [1, 1, 1, 1] ]
all_cluster_diameter_list = []
x, y, z, diameter = [], [], [], []
bin_index_of_points = []
pore_type_count = 0