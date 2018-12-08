from compute_vornoi import  *
from load import load_traj, reduce
import numpy as np
import matplotlib.pyplot as plt
import timeit


# maybe this would be better organized as a class.
# need to do something to reduce the memory constraints...
# Rethink how this is done..
class Timeseries:
    def __init__(self, ts, temp, aa, astd, av, vstd, vi, vf,aas, astds, vas, vstds):
        self.timestep = ts
        self.sampling = 100
        self.temp = temp
        self.surface_area = aa
        self.area_std = astd
        self.average_volume = av
        self.volume_std = vstd
        # determining the unique vornoi indices and giving them an index...
        self.indexes, self.positions = np.unique(np.reshape(vi, (10*self.sampling, 15)),
                                                                 axis=0, return_inverse=True)
        self.all_average_area, self.all_area_std = self.get_areas(aas, astds)
        self.all_average_volume, self.all_volume_std = self.get_volumes(vas, vstds)

        self.all_top_ten_freq = self.get_top(vf)
        # note this is kind of weird implementation because it is unwrapped to help with the calculations.
        # it just helps by giving all the possible Vornoi indices an index... (I should maybe do that myself.

    def get_index_stats(self, property):
        unwrapped = np.reshape(property, (10 * self.sampling))  # unwrap
        resampled = []
        for index in range(1, len(self.indexes)):
            where = np.where(index == self.positions)[0]  # once again the pythons weird arrays strike again...
            is_in = np.floor_divide(where, 10)
            index_prop = []
            count = 0
            for i in range(0, self.sampling):
                if i in is_in:
                    index_prop.append(float(unwrapped[where[count]]))
                    count += 1
                else:
                    index_prop.append(0)
            resampled.append(index_prop)
        return resampled

    def get_top(self,freq,plot=False):
        resampled_freq =self.get_index_stats(freq)
        print(np.sum(resampled_freq, axis =1))
        return resampled_freq

    def get_areas(self, aa, astd, plot=False):
        area = self.get_index_stats(aa)
        astd = self.get_index_stats(astd)
        if plot:
            for a in area:
                plt.scatter(self.temp, a)
            plt.show()
        return area, astd

    def get_volumes(self,vas, vstds, plot=False):
        volumes = self.get_index_stats(vas)
        vstd = self.get_index_stats(vstds)
        if plot:
            for v in volumes:
                plt.scatter(self.temp, v)
            plt.show()
        return volumes, vstd

    def plot_top(self):
        for r in self.all_top_ten_freq:
            plt.scatter(self.temp,r)
        plt.show()

    def get_info_from_index(self, index):
        # python is dumb as hell sometimes
        pos = [num for num,ind in enumerate(self.indexes) if all(ind == index )][0]
        return self.all_average_area[pos], self.all_area_std[pos], self.all_average_volume[pos],\
               self.all_volume_std[pos],self.all_top_ten_freq[pos]
class Sample:
    def __init__(self,timeseries, quench_rates):
        self.timeseries = timeseries  # list of timeseries objects
        self.qenchrates = quench_rates

    def compare_index(self, index):
        avg_area,area_std, avg_vol, vol_std, freq =[],[],[],[],[]
        for quench in self.timeseries:
            holder = quench.get_info_from_index(index)
            avg_area.append(holder[0])
            area_std.append(holder[1])
            avg_vol.append(holder[2])
            vol_std.append(holder[3])
            freq.append(holder[4])
        return avg_area, area_std, avg_vol, vol_std, freq
