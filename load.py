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

# There is probably a better way to do this. This just seaches for the heading and prints the values below that...
def read_lammps_out(out_file, heading):
    with open(out_file) as f:
        thermo=[]
        read = False
        for line in f:
            if "Loop time" in line:
                read = False
            if read:
                # How many times can you cast one thing....
                thermo.append(list(np.array(line.rstrip("\n").split(), dtype=float)))
            if heading in line:
                read = True
    headings = heading.split()
    print(headings)
    thermo = [x for x in zip(*thermo)]
    return thermo, headings


def reduce(filename,lammps_out,datapoints):
    '''
    Notice that the timestep needs to be a multiple of the traj output...
    :param filename:
    :param lammps_out:
    :param timestep:
    :return:
    '''
    timestep, num_atoms, box_bounds, atom_type, position = load_traj(filename)
    stable_len = 154
    seperation = len(position[stable_len:])/datapoints
    indexes = [int(stable_len+np.floor(i*seperation))for i in range(0,datapoints)]
    indexes = indexes
    p = position[indexes]
    na = np.reshape(num_atoms,len(num_atoms))[indexes]
    bb = box_bounds[indexes]
    at = atom_type[indexes]
    print(np.shape(timestep))
    print(indexes)
    ts = np.reshape(timestep,len(timestep))[indexes]
    thermo, headings = read_lammps_out(lammps_out,"Step Temp Press TotEng PotEng Enthalpy")
    ind=list(np.array(np.multiply(indexes,5),dtype=int)) #  pythp
    print(ind)
    temp = np.reshape(thermo[1], len(thermo[1]))[ind] # python is broken... it hates reslicing...
    return p,na, bb, at, ts, temp


'''
reduce('/home/carter/Documents/Classes/760/FinalProject/1E12Cool/traj.lammpstrj',
       '/home/carter/Documents/Classes/760/FinalProject/1E12Cool/out_1E12.lammps',100)
 '''