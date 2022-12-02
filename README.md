# Automatic compartmentalization
This algorithm divides the unit cell of a MOF into pore compartments


## How to use?
To build your cython file:
`python3 setup_periodic_distance.py build_ext --inplace`


## Required: Perform Zeo++ pore size distribution calculation
Zeo++ link:


## Inputs
1. *vpspts file from Zeo++ pore size distribution calculation
2. Number of pore types

# Outputs
1. pore_type_matrix_with_cluster_labels.csv
2. pore_type_matrix_with_pore_type_labels.csv
3. pore_type_matrix_with_cluster_center_labels.html
4. pore_type_matrix_with_pore_type_labels.html
5. geometric_points_with_cluster_labels_for_pore_type_*Npores.html 

