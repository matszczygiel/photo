#!/usr/bin/python3

import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import scipy.integrate

kvals = 10
include = -1

###############################

def get_fit(infile):
    file = open(infile, 'r')
    inp = file.readlines()
    file.close()

    stride = 7 + kvals
    phind = range(1, len(inp), stride)
    bnind = range(0, len(inp), stride)

    phot = []
    for i in phind:
        words = inp[i].split()
        phot.append(float(words[-1]))

    basis_name = []
    for i in bnind:
        words = inp[i].split()
        basis_name.append(str(words[1]))

    csind = range(6, len(inp), stride)
    crss = []
    kthet = []
    for i in csind:
        cv = []
        kv = []
        for k in range(kvals):
            words = inp[i + k].split()
            kv.append(float(words[0]))
            cv.append(float(words[include]))
        kthet.append(kv)
        crss.append(cv)
        
    data = []
    for a in range(len(phot)):
        kv = []
        cv = []
        for i in range(kvals):
            kv.append(kthet[a][i])
            cv.append(crss[a][i])
        for i in range(2, kvals+1):
            kv.append(np.pi - kthet[a][kvals - i])
            cv.append(crss[a][kvals - i])
        for i in range(1, kvals):
            kv.append(np.pi + kthet[a][i])
            cv.append(crss[a][i])
        for i in range(2, kvals+1):
            kv.append(2 * np.pi - kthet[a][kvals - i])
            cv.append(crss[a][kvals - i])

        data.append( [phot[a], kv, cv] )

    cs = []
    for i in range(len(phot)):
        cs.append(scipy.interpolate.CubicSpline(data[i][1], data[i][2], bc_type='periodic'))

    return basis_name, phot, cs

def integrate(infile, outfile_txt):
    basis, phot, fit = get_fit(infile)
    xs = np.arange(0.0, np.pi / 2, 0.005)

    res = np.empty((len(phot), 2), dtype=float)
    res[:,0] = phot
    for i in range(len(phot)):
        res[i, 1] = 4 * np.pi * scipy.integrate.simps(fit[i](xs) * np.sin(xs), xs)

    res = res[res[:,0].argsort()]

    print(res)
    np.savetxt(outfile_txt, res, header=basis[0], fmt=['%3.1f', '%5.3f'])

##################################################

file_dc = "tot_dip_ci_R0.dat"
file_vc = "tot_vel_ci_R0.dat"
file_dh = "tot_dip_hf_R0.dat"
file_vh = "tot_vel_hf_R0.dat"


integrate("parallel_dip_ci_R0.out", file_dc)
integrate("parallel_vel_ci_R0.out", file_vc)
integrate("parallel_dip_hf_R0.out", file_dh)
integrate("parallel_vel_hf_R0.out", file_vh)



exp_file = "/home/mateusz/Documents/cpp/photo/data/h2/exp_results/total.dat"
exp_data = np.loadtxt(exp_file, skiprows=2, ).transpose()

basis = open(file_dc, 'r').readline().split()[1]



data_dc = np.loadtxt(file_dc).transpose()
data_vc = np.loadtxt(file_vc).transpose()
data_dh = np.loadtxt(file_dh).transpose()
data_vh = np.loadtxt(file_vh).transpose()

xs = np.arange(data_dc[0, 0], data_dc[0, -1], 0.01)

cs_dc = scipy.interpolate.CubicSpline(data_dc[0], data_dc[1], bc_type='natural')
cs_vc = scipy.interpolate.CubicSpline(data_vc[0], data_vc[1], bc_type='natural')
cs_dh = scipy.interpolate.CubicSpline(data_dh[0], data_dh[1], bc_type='natural')
cs_vh = scipy.interpolate.CubicSpline(data_vh[0], data_vh[1], bc_type='natural')


plt.clf()
plt.xlim((xs[0], xs[-1]))
plt.ylim((0, 7.5))

plt.plot(data_dc[0], data_dc[1], 'o', markersize=3)
plt.plot(data_vc[0], data_vc[1], 'o', markersize=3)
plt.plot(data_dh[0], data_dh[1], 'o', markersize=3)
plt.plot(data_vh[0], data_vh[1], 'o', markersize=3)

plt.plot(xs, cs_dc(xs), label='dip CI')
plt.plot(xs, cs_vc(xs), label='vel CI')
plt.plot(xs, cs_dh(xs), label='dip HF')
plt.plot(xs, cs_vh(xs), label='vel HF')


plt.errorbar(exp_data[0], exp_data[1], yerr=0.03*exp_data[1],fmt='o',capsize=2, capthick=1, markersize=3, label='exp')
plt.title(basis)
plt.legend()

plt.savefig("total.png")





