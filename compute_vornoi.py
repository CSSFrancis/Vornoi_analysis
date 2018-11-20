'''
Using tess created by https://github.com/wackywendell
Carter Francis
csfrancis@wisc.edu
'''

from tess import Container
from load import load_traj

t,n, bb, at,ap =load_traj('/home/carter/Documents/Classes/760/FinalProject/FinalProject/try2/traj.lammpstrj')

cntr = Container(ap[1], limits=(100,100,100), periodic=True)


area = [v.normals for v in cntr]
print(area)