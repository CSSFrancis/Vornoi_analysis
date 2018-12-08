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

def get_face_freq(container,edge_thresh):
    vcell_list = [v.vertices() for v in container] # all of the indicies of a face
    vlist_list = [v.face_vertices()for v in container] # of the positions of the indiceis
    r_thresh2 = edge_thresh*edge_thresh
    indexes = []
    for vcell,vlist in zip(vcell_list,vlist_list):
        edges = np.zeros(15)
        for face in vlist:
            n_edge = 0
            for i,index in enumerate(face):
                a = index
                b = face[(i+1)%len(face)] # wrapping to complete face
                dx = vcell[a][0] - vcell[b][0]
                dy = vcell[a][1] - vcell[b][1]
                dz = vcell[a][2] - vcell[b][2]
                r2 = dx*dx+dy*dy+dz*dz
                if r2 > r_thresh2:
                    n_edge +=1
            edges[n_edge-1] = edges[n_edge-1]+1
        indexes.append(list(np.array(edges,dtype=int)))
    print(indexes)
    return indexes


def compute_freq(container,verbose=False):
    face_frequency = [v.face_freq_table() for v in container]
    l = [len(le) for le in face_frequency]
    max_size = 15       # much easier than trying to deal with trailing
    face_frequency_padded = [f+[0]*(max_size-le) for f,le in zip(face_frequency,l)]
    index,freq = np.unique(face_frequency_padded,axis=0, return_counts=True)
    top_ten = sorted(range(len(freq)), key=lambda i: freq[i])[-10:]
    top_ten_freq = np.divide(freq[top_ten], len(container))
    if verbose:
        print("The top 10 Vornoi indexes are:")
        [print(i,f)for i, f in zip(index[top_ten],top_ten_freq)]
        print("The percent sum is:", np.sum(top_ten_freq))

    return index, freq, index[top_ten], top_ten_freq


def determine_surface_area(container, verbose=False):
    surface_area = [v.surface_area() for v in container]
    if verbose:
        print(np.shape(surface_area))
        print("The max is:", max(surface_area))
        print("The min is:", min(surface_area))
        print("The Average is:", np.average(surface_area))
    return surface_area, np.average(surface_area), np.std(surface_area)


def determine_volume(container,verbose=False):
    volume = [v.volume() for v in container]
    if verbose:
        print("The max is:", max(volume))
        print("The min is:", min(volume))
        print("The Average is:", np.average(volume))
    return volume, np.average(volume),np.std(volume)


def characterize_index(indices,container):
    containers = []
    average_areas = []
    area_stds = []
    average_volumes = []
    volume_stds = []
    for index in indices:
        index_container = [c for c in container if list(c.face_freq_table()) == list(np.trim_zeros(index, 'b'))]
        containers.append(index_container)
    for cont in containers:
        _, average_area, area_std = determine_surface_area(cont)
        _, average_volume, volume_std = determine_volume(cont)
        average_areas.append(average_area)
        area_stds.append(area_std)
        average_volumes.append(average_volume)
        volume_stds.append(volume_std)
    return average_areas, area_stds, average_volumes, volume_stds


#Testing

t, n, bb, at, ap = load_traj(TwoE12Cooling)
f_time = len(t)-1
cntr = make_container(bb[f_time],at[f_time],ap[f_time])
#get_freq()
#get_index()
#get_volume()
get_face_freq(cntr,0.1)
print([v.face_freq_table() for v in cntr])

'''
#Eventaully I should make extend the contianer class to just add functionality...
class extended_container(Container):
    def __init__(self,boundingbox, atom_type, atom_pos, radii=[0.23, 0.14], central_atom=2):
        is_cental_atom = atom_type == central_atom
        limit = (boundingbox[0][1] - boundingbox[0][0], boundingbox[1][1] - boundingbox[1][0],
                 boundingbox[2][1] - boundingbox[2][0])
        radi_list = atom_type
        for num, r in enumerate(radii):
            radi_list[atom_type == num + 1] = r
        Container.__init__(self,atom_pos, limits=limit, periodic=True, radii=radi_list)
        self.container = [c for is_c, c in zip(is_cental_atom, self) if is_c]
        _, _, self.indices, self.ind_freq = self.get_freq()

    def get_freq(self, verbose=False):
        face_frequency = [v.face_freq_table() for v in self.container]
        l = [len(le) for le in face_frequency]
        max_size = 15  # much easier than trying to deal with trailing
        face_frequency_padded = [f + [0] * (max_size - le) for f, le in zip(face_frequency, l)]
        index, freq = np.unique(face_frequency_padded, axis=0, return_counts=True)
        top_ten = sorted(range(len(freq)), key=lambda i: freq[i])[-10:]
        top_ten_freq = np.divide(freq[top_ten], len(self.container))
        if verbose:
            print("The top 10 Vornoi indexes are:")
            [print(i, f) for i, f in zip(index[top_ten], top_ten_freq)]
            print("The percent sum is:", np.sum(top_ten_freq))
        return index, freq, index[top_ten], top_ten_freq

    def get_surface_area(self, verbose=False):
        surface_area = [v.surface_area() for v in self.container]
        if verbose:
            print(np.shape(surface_area))
            print("The max is:", max(surface_area))
            print("The min is:", min(surface_area))
            print("The Average is:", np.average(surface_area))
        return surface_area, np.average(surface_area), np.std(surface_area)

    def get_volume(self, verbose=False):
        volume = [v.volume() for v in self.container]
        if verbose:
            print("The max is:", max(volume))
            print("The min is:", min(volume))
            print("The Average is:", np.average(volume))
        return volume, np.average(volume), np.std(volume)

    def get_index(self):
        containers = []
        average_areas = []
        area_stds = []
        average_volumes = []
        volume_stds = []
        for index in self.indices:
            index_container = [c for c in self.container if list(c.face_freq_table()) == list(np.trim_zeros(index, 'b'))]
            containers.append(index_container)
        for cont in containers:
            _, average_area, area_std = cont.get_surface_area()
            _, average_volume, volume_std = determine_volume(cont)
            average_areas.append(average_area)
            area_stds.append(area_std)
            average_volumes.append(average_volume)
            volume_stds.append(volume_std)
        return average_areas, area_stds, average_volumes, volume_stds

'''