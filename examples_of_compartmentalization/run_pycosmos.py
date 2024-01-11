import os
import sys

ciffile = sys.argv[1]
Npore_type = int(sys.argv[2])
eps = float(sys.argv[3])
Nmin = int(sys.argv[4])

input_file = 'input_pycosmos.dat'
print('ciffile: %s, Npore types: %d' %(ciffile, Npore_type))
print('DBSCAN parameters: eps = %1.2f, Nmin = %d' %(eps, Nmin))

pycosmos = '/home/shivam/Desktop/Automatic_compartmentalization/Program/pyCOSMOS/src/main.py'
folder_cif = '/home/shivam/Desktop/Automatic_compartmentalization/Program/pyCOSMOS/new_example/cif_files/'
path_cif = folder_cif + ciffile

def extract_cell_parameters(path_cif):
    """
    Returns Lx, Ly, Lz, alpha, beta, gamma from cif file

    """
    with open(path_cif, 'r') as inp:
        for line in inp:
            if '_cell_length_a' in line:
                Lx = float(line.split()[1])

            if '_cell_length_b' in line:
                Ly = float(line.split()[1])

            if '_cell_length_c' in line:
                Lz = float(line.split()[1])

            if '_cell_angle_alpha' in line:
                alpha = float(line.split()[1])

            if '_cell_angle_beta' in line:
                beta = float(line.split()[1])

            if '_cell_angle_gamma' in line:
                gamma = float(line.split()[1])

    return Lx, Ly, Lz, alpha, beta, gamma


def create_cif_folder(ciffile):
    new_folder = ciffile.replace('.cif','')
    try:
        os.mkdir(new_folder)
    except:
        print(new_folder, 'exists, overwritting files')
    os.chdir(new_folder)

#os.system("python %s %s %1.4f %1.4f %1.4f %1.3f %1.3f %1.3f %d | tee summary.txt" %(pycosmos, vpsdpts_file, Lx, Ly, Lz, alpha, beta, gamma, Npore_type))

def create_input_file():

    with open(input_file, 'w') as out:
        out.write('vpsdpts\t %s\n'  %vpsdpts_file)
        out.write('lx\t %1.4f\n' %Lx)
        out.write('ly\t %1.4f\n' % Ly)
        out.write('lz\t %1.4f\n' %Lz)
        out.write('alpha\t %1.4f\n' %alpha)
        out.write('beta\t %1.4f\n' %beta)
        out.write('gamma\t %1.4f\n' %gamma)
        out.write('npore\t %d\n' %Npore_type)
        out.write('eps\t %1.3f \n' %eps)
        out.write('nmin\t %d \n' %Nmin)

    os.system("cp %s . " %path_cif)
    os.system("python %s %s | tee summary.txt" %(pycosmos, input_file))

if __name__ == '__main__':

    Lx, Ly, Lz, alpha, beta, gamma = extract_cell_parameters(path_cif)
    vpsdpts_file = '/home/shivam/Desktop/Automatic_compartmentalization/Program/pyCOSMOS/new_example/vpsdpts/' + ciffile.replace('.cif','.vpsdpts') 

    create_cif_folder(ciffile)
    create_input_file()
