# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 10:54:12 2021

@author: uqxzho15
"""

################################################################################

import math, xlrd, xlwt, os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

################################################################################

l = 1800          # area size is l*l.
r = 20            # bipolar axon arbor radius.
mode = 2          # 1 for regular, 2 for random.
tess = 2          # 1 for sqaure, 2 for regular hexagon.
regmindis = 10    # the minimum distance between neighbouring axon arbor in regular mode.
pct = 1.00        # percentage of area that bipolar axon arbors occupy.
rep = 1        # how many repetitions for each annulus setting. 0< rep <= 1000.
ovd = 3
bipolar_writing_structure = 7

# read previously generated bipolar arrays.

bipolarfile = "'C:/Users/uqxzho15/Desktop/py code/bipolar sim 3 sustained off alpha/200 arrays 1800.xls'"
bipolarinput = bipolarfile[1:-1]

bipolarwb = xlrd.open_workbook(bipolarinput)

for n in range(len(bipolarwb.sheet_names())):
    if bipolarwb.sheet_names()[n] == 'Sheet 1':
        bps = bipolarwb.sheets()[n]
        bpline = []
        for row in range(bps.nrows):
            col_value = []
            for col in range(bps.ncols):
                col_value.append(bps.cell(row,col).value)
            bpline.append(col_value)
            
bipolarlist = []
for i in range(len(bpline)):
    if bpline[i][0] == 'array':
        bipolar = []
        for n in range(bipolar_writing_structure):
            for t in range(len(bpline[i+2+2*n])):
                if bpline[i+2+2*n][t] != '':
                    bipolar.append((bpline[i+2+2*n][t], bpline[i+2+2*n+1][t]))
        bipolarlist.append(bipolar)
bipolarselect = bipolarlist[0:rep]

################################################################################


for bipolar in bipolarselect:

    fig, ax = plt.subplots()
    ax = plt.gca(aspect='equal')
    ax.cla()  
    ax.set_xlim((-900, 900))
    ax.set_ylim((-900, 900))
                
    output_list = []
    for cell in bipolar:
        output_list.append(plt.Circle((cell[0], cell[1]), r, color='grey', fill=False, linewidth=0.1))

    for cell in output_list:
        ax.add_artist(cell)
                
    ax.axis('off')
    figname = 'bipolar array illustration - array 1.svg'
    fig.savefig(figname, format = 'svg', dpi = 1200)
            










