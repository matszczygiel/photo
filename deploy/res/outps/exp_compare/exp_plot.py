#!/usr/bin/python3

import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate

kvals = 10
include = -2


infiles = sys.argv[1:]
exppath = "/home/mateusz/Documents/cpp/photo/data/h2/exp_results/"
exp_phot = [21.2, 23.1, 26.9, 40.8]

##############################
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
######################

for infile in infiles:
    savefile = infile[:-4]

    basis, phot, cs = get_fit(infile)
    xs = np.arange(0.0, 2*np.pi, 0.01)

    for i in range(len(phot)):
        if phot[i] not in exp_phot:
            continue
            
        exp_file = exppath
        if infile.find("parall") != -1:
            exp_file += "parallel_" + str(phot[i]) + "eV.dat"
        elif infile.find("perpend") != -1:
            exp_file += "perpendicular_" + str(phot[i]) + "eV.dat"

        exp = np.transpose(np.loadtxt(exp_file))

        plt.clf()
        plt.figure(figsize=(13, 9))

        ax = plt.subplot(111, projection='polar')
#        ax.plot(xs, cs[i](xs), label=str(phot[i]))
        cs_max = max(cs[i](xs))
        ax.plot(xs, cs[i](xs)/cs_max, label=str(phot[i]))
        if len(exp) != 0: 
#            ax.errorbar(exp[0] * 2* np.pi / 360, exp[1]/100, yerr=exp[2]/100,fmt='o',capsize=2, capthick=1, markersize=3, label='exp')
            max_exp = np.average(np.sort(exp[1])[-4])
            ax.errorbar(exp[0] * 2* np.pi / 360, exp[1]/max_exp, yerr=exp[2]/max_exp,fmt='o',capsize=2, capthick=1, markersize=3, label='exp')
        ax.set_rlabel_position(-22.5)
        ax.set_thetagrids(range(0,360, 15))
        ax.set_rmin(0.)
        ax.set_rorigin(0.)
        ax.grid(True)
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        ax.set_title(basis[i] + "         " + savefile, va='bottom')
        plt.savefig(savefile + "_" + str(phot[i]) + ".png")













