'''
Using tess created by https://github.com/wackywendell
Carter Francis
csfrancis@wisc.edu
'''
import numpy as np

from tess import Container
from load import load_traj

import timeit
# input position files...
TwoE12Cooling = '/Users/shaw/Shaw/MSE760/FinalProject/2E12Cool/traj.lammpstrj'
OneE12Cooling = '/home/carter/Documents/Classes/760/FinalProject/1E12Cool/traj.lammpstrj'
FiveE11Cooling = '/Users/shaw/Shaw/MSE760/FinalProject760(ExternalCluster)/5E11Cool/traj.lammpstrj'
TwoE11Cooling = '/Users/shaw/Shaw/MSE760/FinalProject760(ExternalCluster)/2E11Cool_2/traj.lammpstrj'


def get_radii(atom_list, radii):
    radi_list = atom_list
    for num,r in enumerate(radii):
        radi_list[atom_list==num+1]=r
    return radi_list


def make_container(boundingbox, atom_type, atom_pos, radii=[0.23, 0.14], central_atom=2):
    is_cental_atom = atom_type == central_atom
    limit = (boundingbox[0][1] - boundingbox[0][0], boundingbox[1][1] - boundingbox[1][0],
             boundingbox[2][1] - boundingbox[2][0])
    radii = get_radii(atom_type, radii)
    cntr = Container(atom_pos, limits=limit, periodic=True, radii=radii)
    cntr = [c for is_c, c in zip(is_cental_atom, cntr) if is_c]
    return cntr


def compute_freq(container,verbose=False):
    face_frequency = [v.face_freq_table() for v in container]
    l = [len(le) for le in face_frequency]
    max_size = max(l)
    face_frequency_padded = [f+[0]*(max_size-le) for f,le in zip(face_frequency,l)]
    index,freq = np.unique(face_frequency_padded,axis=0, return_counts=True)
    top_ten = sorted(range(len(freq)), key=lambda i: freq[i])[-10:]
    top_ten_freq = np.divide(freq[top_ten], len(container))
    if verbose:
        print("The top 10 Vornoi indexes are:")
        [print(i,f)for i, f in zip(index[top_ten],top_ten_freq)]
        print("The percent sum is:", np.sum(top_ten_freq))

    return index, freq, top_ten


def determine_surface_area(container):
    surface_area = [v.surface_area() for v in container]
    print(np.shape(surface_area))
    print("The max is:", max(surface_area))
    print("The min is:", min(surface_area))
    print("The Average is:", np.average(surface_area))
    return surface_area


def determine_volume(container):
    volume = [v.volume() for v in container]
    print("The max is:", max(volume))
    print("The min is:", min(volume))
    print("The Average is:", np.average(volume))
    return volume


def determine_index_surface_area(indices,container):
    containers = []

    for index in indices:
        index_container = [c for c in container if list(c.face_freq_table()) == list(np.trim_zeros(index, 'b'))]
        containers.append(index_container)


    for cont in containers:
        determine_surface_area(cont)
        determine_volume(cont)

#Testing
'''
t, n, bb, at, ap = load_traj(OneE12Cooling)
f_time = len(t)-1
cntr = make_container(boundingbox=bb[f_time],atom_type=at[f_time],atom_pos=ap[f_time])
print(cntr)
determine_surface_area(cntr)
i,f,top = compute_freq(cntr,verbose=True)
determine_index_surface_area(i[top],container=cntr)
'''