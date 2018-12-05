'''
Using tess created by https://github.com/wackywendell
Carter Francis
csfrancis@wisc.edu
'''
import numpy as np

from tess import Container
from load import load_traj

# input position files...
TwoE12Cooling = '/home/carter/Documents/Classes/760/FinalProject/2E12Cool/traj.lammpstrj'
OneE12Cooling = '/home/carter/Documents/Classes/760/FinalProject/1E12Cool/traj.lammpstrj'
t, n, bb, at, ap = load_traj(OneE12Cooling)

print(len(t))
Time = 744
limit = (bb[Time][0][1]-bb[Time][0][0],bb[Time][1][1]-bb[Time][1][0],bb[Time][2][1]-bb[Time][2][0])
#must have the right limits from the bounding box....
cntr = Container(ap[Time], limits=limit, periodic=True)
is_copper = at[Time] == 2
con = [c for is_c,c in zip(is_copper,cntr) if is_c]
def compute_freq(container):
    face_frequency = [v.face_freq_table() for v in container]
    l = [len(le) for le in face_frequency]
    max_size = max(l)
    face_frequency_padded = [f+[0]*(max_size-le) for f,le in zip(face_frequency,l)]
    index,freq = np.unique(face_frequency_padded,axis=0, return_counts=True)
    return index, freq
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

determine_surface_area(con)
i,f = compute_freq(con)

print(len(con))
print(i[0])
print(f)
top_ten = sorted(range(len(f)), key=lambda i: f[i])[-10:]
print(i[top_ten])
print(np.sum(np.divide(f[top_ten], len(con))))
print(top_ten)
determine_index_surface_area(i[top_ten],container=con)