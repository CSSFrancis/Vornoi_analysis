# this program is designed to preform the analysis and output the data in a way

from load import load_traj, reduce
from compute_vornoi import *
from statistical_analysis import Timeseries, Sample
import matplotlib.pyplot as plt
def output_analysis(traj_file, themo_file, data_points, outfile):
    pos, num_atom, boundbox, atom_type, timestep, temp = reduce(traj_file, themo_file, data_points)
    with open(outfile, "w") as f:
        for bb, at, p,t,ts in zip(boundbox, atom_type, pos,temp, timestep):
            cont = make_container(bb, at, p)  # make sure to overwrite the container so that there isn't memory leakage

            _, _, indices, top_ten_freq = compute_freq(cont)
            _, average_sarea, std_sarea = determine_surface_area(cont)
            _, average_vol, std_vol = determine_volume(cont)
            average_areas, area_stds, average_volumes, volume_stds = characterize_index(indices, cont)
            f.write("Timestep:" + ts+"\n")
            f.write("Temperature: %g\n" %t)
            f.write("AverageArea: %g\n" % average_sarea)
            f.write("AreaStd: %g\n" % std_sarea)
            f.write("AverageVolume: %g\n" %average_vol)
            f.write("VolumeStd: %g\n" %std_vol)
            f.write("VornoiIndicies:")
            for index in indices:
                for i in index:
                    f.write("%g " % i)
                f.write(",")
            f.write("\n")
            write_list("VornoiFreq:", f, top_ten_freq)
            write_list("AverageAreas:",f, average_areas)
            write_list("AreaStds:", f, area_stds)
            write_list("AverageVolumes:", f, average_volumes)
            write_list("VolumeStds:", f, volume_stds)

    return
def write_list(name, file, list):
    file.write(name)
    [file.write("%g "%l) for l in list]
    file.write("\n")
    return


def read_analysis(input_file):
    '''
    Returns a timeseries object which can be manipulated for further analysis.
    '''
    ts, temp, aa, astd, av, vstd, vi, vf,aas, astds,vas,vstds = [], [], [], [], [], [], [], [], [], [], [], []

    with open(input_file) as f:
        for line in f:
            if "Timestep:" in line:
                ts.append(int(line.strip("\n").split(":")[1]))
            elif "Temperature:" in line:
                temp.append(float(line.strip("\n").split(":")[1]))
            elif "AverageArea:" in line:
                aa.append(float(line.strip("\n").split(":")[1]))
            elif "AreaStd:" in line:
                astd.append(float(line.strip("\n").split(":")[1]))
            elif "AverageVolume:" in line:
                av.append(float(line.strip("\n").split(":")[1]))
            elif "VolumeStd:" in line:
                vstd.append(float(line.strip("\n").split(":")[1]))
            elif "VornoiIndicies:" in line:
                vi.append([line.split()for line in line.strip("\n").split(":")[1].split(",")][0:10])
            elif "VornoiFreq:" in line:
                vf.append(list(np.array(line.strip("\n").split(":")[1].split())))
            elif"AverageAreas:"in line:
                aas.append(list(np.array(line.strip("\n").split(":")[1].split())))
            elif "AreaStds:" in line:
                astds.append(list(np.array(line.strip("\n").split(":")[1].split())))
            elif "AverageVolumes:"in line:
                vas.append(list(np.array(line.strip("\n").split(":")[1].split())))
            elif "VolumeStds:"in line:
                vstds.append(list(np.array(line.strip("\n").split(":")[1].split(),dtype=float)))
        vi = [[[int(x) for x in y] for y in v]for v in vi]
        timeSeries = Timeseries(ts, temp, aa, astd, av, vstd, vi, vf,aas, astds,vas,vstds)

    return timeSeries


OneE12Cooling_lammps = '/home/carter/Documents/Classes/760/Organized/1E12Cooling_5000Atoms.lammps'
OneE12Cooling_traj = '/home/carter/Documents/Classes/760/Organized/1E12Cooling_5000Atoms.lammpstrj'
TwoE12Cooling_lammps ='/home/carter/Documents/Classes/760/Organized/2E12Cooling_5000Atoms.lammps'
TwoE12Cooling_traj = '/home/carter/Documents/Classes/760/Organized/2E12Cooling_5000Atoms.lammpstrj'
FiveE11Cooling_lammps ='/home/carter/Documents/Classes/760/Organized/5E11Cooling_5000Atoms.lammps'
FiveE11Cooling_traj = '/home/carter/Documents/Classes/760/Organized/5E11Cooling_5000Atoms.lammpstrj'
FiveE11Cooling_lammps_2 ='/home/carter/Documents/Classes/760/Organized/5E11Cooling_5000Atoms_1.lammps'
FiveE11Cooling_traj_2 ='/home/carter/Documents/Classes/760/Organized/5E11Cooling_5000Atoms_1.lammpstrj'
TwoE11Cooling_lammps ='/home/carter/Documents/Classes/760/Organized/2E11Cooling_5000Atoms.lammps'
TwoE11Cooling_traj ='/home/carter/Documents/Classes/760/Organized/2E11Cooling_5000Atoms.lammpstrj'
TwoE11Cooling_traj_Higher_Q ='/home/carter/Documents/Classes/760/Organized/2E11Cooling_10000AtomsHigherQ.lammpstrj'
TwoE11Cooling_lammps_Higher_Q ='/home/carter/Documents/Classes/760/Organized/2E11Cooling_10000AtomsHigherQ.lammps'

timeseries_list = ["TwoE12Cooling.out", "OneE12Cooling.out", "FiveE11Cooling.out", "TwoE11Cooling_HQ.out"]
ts_list = []
#t_st = read_analysis("TwoE11Cooling_HQ.out")
#avg_area, area_std, avg_vol, vol_std, freq = t_st.get_info_from_index(index=[0, 0, 0, 0, 0, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0])
#print("The average area", avg_area)
for t in timeseries_list:
    ts_list.append(read_analysis(t))
s=Sample(ts_list, timeseries_list)
avg_area, area_std, avg_vol, vol_std, freq = s.compare_index(index=[0, 0, 0, 0, 0, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0])
[plt.scatter(range(0,100),a) for a in vol_std]
plt.show()

#ts = read_analysis("out.out")
#ts.plot_top()
#print(ts.all_volume_std)
#ts.get_info_from_index([0,0,0,0,0,12,0,0,0,0,0,0,0,0,0])
#output_analysis(FiveE11Cooling_traj_2,FiveE11Cooling_lammps_2, 100, 'FiveE11Cooling_2.out')
