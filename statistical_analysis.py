from compute_vornoi import  *
from load import load_traj, reduce
import numpy as np
import matplotlib.pyplot as plt
import timeit


# maybe this would be better organized as a class.
# need to do something to reduce the memory constraints...
# Rethink how this is done..
class Timeseries:
    def __init__(self, traj_file,lammps_file, sampling=10):
        self.sampling = sampling
        pos, natom, bobo, atype, tstep, self.temp = reduce(traj_file, lammps_file, self.sampling)
        self.containers = [make_container(bb, at, ap) for bb, at, ap in zip(bobo, atype, pos)]
        self.all_average_area =[]
        self.all_area_std = []
        self.all_average_volume = []
        self.all_volume_std = []
        self.all_top_ten = []
        self.all_top_ten_freq =[]
        for container in self.containers:
            _, _, top_ten, top_ten_freq = compute_freq(container)
            average_area, area_std, average_volume, volume_std = characterize_index(top_ten, container)
            self.all_average_area.append(average_area)
            self.all_area_std.append(area_std)
            self.all_average_volume.append(average_volume)
            self.all_volume_std.append(volume_std)
            self.all_top_ten.append(top_ten)
            self.all_top_ten_freq.append(top_ten_freq)
        # note this is kind of weird implementation because it is unwrapped to help with the calculations.
        # it just helps by giving all the possible Vornoi indices an index... (I should maybe do that myself.
        self.indexes, self.positions, self.frequency = np.unique(np.reshape(self.all_top_ten, (10*sampling, 15)), axis=0,
                                                                 return_inverse=True,return_counts=True)

    def get_top(self):
        freq = np.reshape(self.all_top_ten_freq, (10*self.sampling)) # unwrap
        resampled_freq = []
        for index in range(1, len(self.indexes)):
            where = np.where(index == self.positions)[0]  # once again the pythons weird arrays strike again...
            is_in = np.floor_divide(where,10)
            index_freq = []
            count = 0
            for i in range(0,10):
                if i in is_in:
                    index_freq.append(freq[where[count]])
                    count += 1
                else:
                    index_freq.append(0)
            resampled_freq.append(index_freq)
        return resampled_freq

    def plot_top(self):
        rf = self.get_top()
        print(self.temp)
        for r in rf:
            plt.scatter(self.temp,r)
        plt.show()


def analyze_timeseries(traj_file,lammps_file):
    pos, natom, bobo, atype, tstep, temp = reduce(traj_file,lammps_file, 10)

    containers = [make_container(bb, at, ap) for bb, at, ap in zip(bobo, atype, pos)]
    all_average_area =[]
    all_area_std = []
    all_average_volume = []
    all_volume_std = []
    for container in containers:
        _, _, top_ten, top_ten_freq = compute_freq(container)
        average_area, area_std, average_volume, volume_std = characterize_index(top_ten, container)
        all_average_area.append(average_area)
        all_area_std.append(area_std)
        all_average_volume.append(average_volume)
        all_volume_std.append(volume_std)
    return



    print(np.shape(containers))



timeSeries  = Timeseries('/home/carter/Documents/Classes/760/FinalProject/1E12Cool/traj.lammpstrj',
                   '/home/carter/Documents/Classes/760/FinalProject/1E12Cool/out_1E12.lammps')
print(timeSeries.get_top())
timeSeries.plot_top()