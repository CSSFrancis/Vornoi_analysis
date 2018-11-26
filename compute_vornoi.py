'''
Using tess created by https://github.com/wackywendell
Carter Francis
csfrancis@wisc.edu
'''
import numpy as np

from tess import Container
from load import load_traj


t,n, bb, at,ap =load_traj('/Users/shaw/Shaw/MSE760/FinalProject/try2/traj.lammpstrj')


Time = 27
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
i,f = compute_freq(con)
print(i[0])
print(f)
top_ten = sorted(range(len(f)), key=lambda i: f[i])[-15:]
print(i[top_ten])
print(np.divide(f[top_ten], len(con)))

