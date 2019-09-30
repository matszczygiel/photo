#!/usr/bin/python3

import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate

include = -1


def plot_from_file(infile, llabel):
    file = open(infile, 'r')
    inp = file.readlines()
    file.close()

    stride = 7
    phind = range(1, len(inp), stride)

    phot = []
    for i in phind:
        words = inp[i].split()
        phot.append(float(words[-1]))

    csind = range(6, len(inp), stride)
    crss = []
    for i in csind:
        words = inp[i].split()
        crss.append(float(words[include]))

    xs = np.arange(phot[0], phot[-1], 0.01)
    cs = scipy.interpolate.CubicSpline(phot, crss)

    plt.xlim((phot[0], phot[-1]))
    plt.plot(phot, crss, 'o', markersize=3)
    plt.plot(xs, cs(xs), label=llabel)


exp_file = "/home/mateusz/Documents/photo/he/he_exp_data.dat"
exp_data = np.loadtxt(exp_file, skiprows=2, ).transpose()


def plot_exp():
    plt.errorbar(exp_data[0], exp_data[1], yerr=0.03*exp_data[1],
                 fmt='o', capsize=2, capthick=1, markersize=3, label='exp')


basis = open("cNbN_dip_ci.out", 'r').readline().split()[1]

plt.clf()
plt.ylim((1, 7.5))
plot_from_file("cNbN_dip_ci.out", llabel='dip CI')
plot_from_file("cNbN_vel_ci.out", llabel='vel CI')
plot_from_file("cNbN_dip_hf.out", llabel='dip HF')
plot_from_file("cNbN_vel_hf.out", llabel='vel HF')
plot_exp()


plt.title(basis, va='bottom')
plt.legend()
plt.savefig('total.png')
