# this program is designed to preform the analysis and output the data in a way

from load import load_traj, reduce
from compute_vornoi import *
from statistical_analysis import Timeseries

def output_analysis(traj_file, themo_file, data_points, outfile):
    pos, num_atom, boundbox, atom_type, timestep, temp = reduce(traj_file,themo_file,data_points)
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
        timeSeries = Timeseries(ts, temp, aa, astd, av, vstd, vi, vf,aas, astds,vas,vstds)

    return timeSeries



ts = read_analysis("out.out")
ts.plot_top()
#output_analysis('/home/carter/Documents/Classes/760/FinalProject/1E12Cool/traj.lammpstrj',
#               '/home/carter/Documents/Classes/760/FinalProject/1E12Cool/out_1E12.lammps',100, 'out.out')