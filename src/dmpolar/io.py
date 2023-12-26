import numpy as np


# read & write QM input and output
def read_xyz(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    tot_lines = len(lines)
    natoms = int(lines[0])
    tot_coords = []
    ptr = 0
    while ptr < tot_lines:
        atoms = []
        coords = []
        for line in lines[ptr+2:ptr+2+natoms]:
            atom, x, y, z = line.split()
            atoms.append(atom)
            coords.append([float(x), float(y), float(z)])
        tot_coords.append(coords)
        ptr += natoms + 2
    return atoms, np.array(tot_coords)