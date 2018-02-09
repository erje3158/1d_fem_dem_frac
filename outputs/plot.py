import matplotlib.pyplot as plt
import numpy as np
import subprocess as sp
import os

def plot_disp(path_to_data, num_of_fig, conv_factors):

    print("----> Displacement Plot")

    base_dir = os.getcwd()
    os.chdir(path_to_data)

    with open('d.txt') as f:
        lines = f.readlines()
        dstr = [line.split()[0] for line in lines]

    #d = np.zeros(len(dstr))
    d = []
    for i in dstr:
        d.append(float(i))
    #d = d * conv_factors[0]
    #print(d[:10])

    with open('t.txt') as f:
        lines = f.readlines()
        tstr = [line.split()[0] for line in lines]

    t = []
    for i in tstr:
        t.append(float(i))

    plt.figure(num_of_fig)
    plt.plot(t,d)

    plt.xlabel('Time (ms)')
    plt.ylabel('Displacement (mm)')

    plt.savefig("disp_" + str(num_of_fig) + ".pdf")

    os.chdir(base_dir)





def main():
    print("\n\n\n")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("                  1D FEM-DEM PLOTS                    ")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("\n\n\n")

    dispfactor   = 10**-3 #mm
    stressfactor = 10**-6 #MPa
    timefactor   = 10**-3 #ms

    factors = [dispfactor, stressfactor, timefactor]

    plot_disp("./debug",1, factors)


if __name__ == "__main__": 
    main()
