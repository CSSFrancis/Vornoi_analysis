'''
This program is designed to preform the analysis and output the data in a way that is more easily read by Origin/excel
Carter Francis
csfrancis@wisc.edu
'''

from load import load_traj, reduce
from compute_vornoi import *
from statistical_analysis import Timeseries, Sample
import matplotlib.pyplot as plt


def output_analysis(traj_file, themo_file, data_points, outfile):
    pos, num_atom, boundbox, atom_type, timestep, temp = reduce(traj_file, themo_file, data_points)
    with open(outfile, "w") as f:
        for bb, at, p, t, ts, na in zip(boundbox, atom_type, pos,temp, timestep, num_atom):
            cont = make_container(bb, at, p)  # make sure to overwrite the container so that there isn't memory leakage
            density = (float(na)/.6022)*(91.224*.36+63.546*.64)/(float((bb[0][1]-bb[0][0])**3))
            _, _, indices, top_ten_freq = compute_freq(cont)
            _, average_sarea, std_sarea = determine_surface_area(cont)
            _, average_vol, std_vol = determine_volume(cont)
            average_areas, area_stds, average_volumes, volume_stds = characterize_index(indices, cont)
            f.write("Timestep:" + ts+"\n")
            f.write("Temperature: %g\n" %t)
            f.write("Density: %g\n" % density)
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
    ts, temp, aa, astd, av, vstd, vi, vf,aas, astds,vas,vstds,den = [], [], [], [], [], [], [], [], [], [], [], [], []

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
            elif "Density:" in line:
                den.append(float(line.strip("\n").split(":")[1]))
        vi = [[[int(x) for x in y] for y in v]for v in vi]
        timeSeries = Timeseries(ts, temp, aa, astd, av, vstd, vi, vf,aas, astds,vas,vstds,den)

    return timeSeries

def print_index_details(index,out_file,sample):
    avg_area, area_std, avg_vol, vol_std, freq, temp = sample.compare_index(index=index)
    with open(str(index)+out_file+"Average_Area.csv", "w") as f:
        f.write("2E12Cooling,1E12Cooling,5E11Cooling,2E11Cooling,2E11Cooling_1000,2E11Cooling_HQ\n")
        np.savetxt(f, list(zip(*avg_area)), delimiter = ',')
    with open(str(index)+out_file + "Area_std.csv", "w") as f:
        f.write("2E12Cooling,1E12Cooling,5E11Cooling,2E11Cooling,2E11Cooling_1000,2E11Cooling_HQ\n")
        np.savetxt(f, list(zip(*area_std)), delimiter=',')
    with open(str(index)+out_file + "Average_Vol.csv", "w") as f:
        f.write("2E12Cooling,1E12Cooling,5E11Cooling,2E11Cooling,2E11Cooling_1000,2E11Cooling_HQ\n")
        np.savetxt(f, list(zip(*avg_vol)), delimiter=',')
    with open(str(index)+out_file + "Vol_std.csv", "w") as f:
        f.write("2E12Cooling,1E12Cooling,5E11Cooling,2E11Cooling,2E11Cooling_1000,2E11Cooling_HQ\n")
        np.savetxt(f, list(zip(*vol_std)), delimiter=',')


def print_density(outfile,sample):
    ts = sample.timeseries
    with open(outfile, "w") as f:
        density = [t.density for t in ts]
        f.write("2E12Cooling,1E12Cooling,5E11Cooling,2E11Cooling,2E11Cooling_1000,2E11Cooling_HQ\n")
        np.savetxt(f, list(zip(*density)), delimiter=',')


def print_average_values(outfile,sample):
    ts = sample.timeseries
    with open(outfile+"Average_Area.csv", "w") as f:
        aa = [t.surface_area for t in ts]
        f.write("2E12Cooling,1E12Cooling,5E11Cooling,2E11Cooling,2E11Cooling_1000,2E11Cooling_HQ\n")
        np.savetxt(f, list(zip(*aa)), delimiter=',')
    with open(outfile + "Area_STD.csv", "w") as f:
        ast = [t.area_std for t in ts]
        f.write("2E12Cooling,1E12Cooling,5E11Cooling,2E11Cooling,2E11Cooling_1000,2E11Cooling_HQ\n")
        np.savetxt(f, list(zip(*ast)), delimiter=',')
    with open(outfile + "Average_vol.csv", "w") as f:
        av = [t.average_volume for t in ts]
        f.write("2E12Cooling,1E12Cooling,5E11Cooling,2E11Cooling,2E11Cooling_1000,2E11Cooling_HQ\n")
        np.savetxt(f, list(zip(*av)), delimiter=',')
    with open(outfile + "vol_STD.csv", "w") as f:
        vst = [t.volume_std for t in ts]
        f.write("2E12Cooling,1E12Cooling,5E11Cooling,2E11Cooling,2E11Cooling_1000,2E11Cooling_HQ\n")
        np.savetxt(f, list(zip(*vst)), delimiter=',')
    with open(outfile + "temp.csv", "w") as f:
        vst = [t.temp for t in ts]
        f.write("2E12Cooling,1E12Cooling,5E11Cooling,2E11Cooling,2E11Cooling_1000,2E11Cooling_HQ\n")
        np.savetxt(f, list(zip(*vst)), delimiter=',')

def print_freq(out_file,sample):
    ts =sample.timeseries
    with open(out_file,"w") as f:
        for t in ts:
            print(list(zip(*t.all_top_ten_freq))[99])
            print(t.all_top_ten_freq)
'''
# These are just the generalized input files from my computer...

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
TwoE11Cooling_lammps_10000 ='/home/carter/Documents/Classes/760/Organized/2E11Cooling_10000Atoms.lammps'
TwoE11Cooling_traj_10000 ='/home/carter/Documents/Classes/760/Organized/2E11Cooling_10000Atoms.lammpstrj'
TwoE11Cooling_traj_Higher_Q ='/home/carter/Documents/Classes/760/Organized/2E11Cooling_10000AtomsHigherQ.lammpstrj'
TwoE11Cooling_lammps_Higher_Q ='/home/carter/Documents/Classes/760/Organized/2E11Cooling_10000AtomsHigherQ.lammps'
OneE11Cooling_lammps = '/home/carter/Documents/Classes/760/Organized/1E11Cooling_5000Atoms.lammps'
OneE11Cooling_traj = '/home/carter/Documents/Classes/760/Organized/1E11Cooling_5000Atoms.lammpstrj'

# for outputting analysis
output_analysis(TwoE11Cooling_traj, TwoE11Cooling_lammps, 100, 'OutputFiles/TwoE11Cooling.out')
output_analysis(TwoE11Cooling_traj, TwoE11Cooling_lammps, 100, 'OutputFiles/OneE11Cooling.out')
output_analysis(FiveE11Cooling_traj, FiveE11Cooling_lammps, 100, 'OutputFiles/FiveE11Cooling.out')
output_analysis(TwoE12Cooling_traj, TwoE12Cooling_lammps, 100, 'OutputFiles/TwoE12Cooling.out')
output_analysis(OneE12Cooling_traj, OneE12Cooling_lammps, 100, 'OutputFiles/OneE12Cooling.out')
output_analysis(TwoE11Cooling_traj_10000, TwoE11Cooling_lammps_10000, 100, 'OutputFiles/TwoE11Cooling_10000.out')
output_analysis(TwoE11Cooling_traj_Higher_Q, TwoE11Cooling_lammps_Higher_Q, 100, 'OutputFiles/TwoE11Cooling_HQ.out')

'''
# list of timeseries outputs
timeseries_list = ["OutputFiles/TwoE12Cooling.out", "OutputFiles/OneE12Cooling.out", "OutputFiles/FiveE11Cooling.out",
                   "OutputFiles/TwoE11Cooling.out","OutputFiles/TwoE11Cooling_10000.out",
                   "OutputFiles/TwoE11Cooling_HQ.out"]
ts_list = []

for t in timeseries_list:
    ts_list.append(read_analysis(t))
s = Sample(ts_list, timeseries_list)  # creating a sample class which is basically multiple time series
avg_area, area_std, avg_vol, vol_std, freq, temp = s.compare_index(index=[0, 0, 0, 0, 2, 8, 1, 0, 0, 0, 0, 0, 0, 0, 0])

# out putting a bunch of files
print_index_details([0, 0, 0, 0, 1, 10, 2, 0, 0, 0, 0, 0, 0, 0, 0],'out',s)
print_index_details([0, 0, 0, 0, 2, 8, 1, 0, 0, 0, 0, 0, 0, 0, 0],'out',s)
print_index_details([0, 0, 0, 0, 2, 8, 2, 0, 0, 0, 0, 0, 0, 0, 0],'out',s)
print_index_details([0, 0, 0, 0, 3, 6, 3, 0, 0, 0, 0, 0, 0, 0, 0],'out',s)
print_density("density.csv",s)
print_average_values("SystemVarible",s)
print_freq("freq.csv",s)
[plt.scatter(t, a) for a, t in zip(avg_area,temp)]
plt.xlim(2000, 0)
plt.show()

