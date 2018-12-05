'''
Loading the positions from a traj file. May add additional support for different file types.
Carter Francis
csfrancis@wisc.edu
'''
import numpy as np

def load_traj(file_name):
    with open (file_name) as traj:
        file_list = [line.strip(" \n").split(" ") for line in traj]
    timestep = []
    num_atoms = []
    box_bounds = []
    atoms_pos = []
#This shouldn't work and I don't know why it does but whatever....
    for line in file_list:
        if "ITEM:" in line:
            holder = []
            id = line[1]
            if id is not None:
                if "TIMESTEP" == id:
                    timestep.append(holder)
                elif "NUMBER" == id:
                    num_atoms.append(holder)
                elif "BOX" == id:
                    box_bounds.append(holder)
                elif "ATOMS" == id:
                    atoms_pos.append(holder)
        else:
            holder.append(line)
    position = []
    atom_type = []
    for positions in atoms_pos:
        position_holder = []
        atom_holder = []
        for pos in positions:
            atom_holder.append(pos[1])
            position_holder.append(pos[2:5])
        atom_type.append(atom_holder)
        position.append(position_holder)
    position=np.float_(position)
    box_bounds = np.float_(box_bounds)
    atom_type = np.float_(atom_type)
    return timestep, num_atoms, box_bounds, atom_type,position
