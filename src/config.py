from math import ceil, radians
import sys
import numpy as np

INPUT = sys.argv[1] # vpsdpts file
Lx = float(sys.argv[2])
Ly = float(sys.argv[3])
Lz = float(sys.argv[4])
alpha_degree = float(sys.argv[5])
beta_degree = float(sys.argv[6])
gamma_degree = float(sys.argv[7])
NumberOfPoreTypes = int(sys.argv[8])
# output file
File_pore_type_matrix = "pore_type_matrix.csv"

# convert angles to radians
alpha, beta, gamma = radians(alpha_degree), radians(beta_degree), radians(gamma_degree)
# calculate volume of unit cell
Volume_of_uc = Lx*Ly*Lz*np.sqrt(1-np.cos(alpha)**2 - np.cos(beta)**2
                -np.cos(gamma)**2 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))

print('INPUT is ', INPUT)
print("Unit cell: %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f" %(Lx, Ly, Lz, alpha_degree, beta_degree, gamma_degree))
print("Number of Pore types: ", NumberOfPoreTypes)
print("Volume of uc [A^3]: %1.3f "%Volume_of_uc)

pore_type_matrix_with_pore_type_labels = -1 * np.ones((int(ceil(Lx)), int(ceil(Ly)), int(ceil(Lz)))) # 3D pore type matrix
pore_type_matrix_with_cluster_labels = -1 * np.ones((int(ceil(Lx)), int(ceil(Ly)), int(ceil(Lz)))) # 3D pore type matrix
all_cluster_center_list = []  # [ [c1, c2], [c3, c4, c5, c6] ......... ] where ci = np.array([xi, yi, zi])
all_cluster_pore_type_labels = [] # [[0, 0,], [1, 1, 1, 1] ]
all_cluster_diameter_list = []
bin_index_of_points = []
pore_type_count = 0
# density of geometric points (1/Ang^3)
rho = 50000/Volume_of_uc
