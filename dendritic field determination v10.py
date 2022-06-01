# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 14:16:25 2020

@author: uqxzho15

system log:
v.01
1. develop a way to define the contour of dendritic field inside out.
v.02
1. apply bipolar array
2. plot contour, BPs within dend field, and connected BPs anainst unconnected.
v.03
1. save all figures as vector plot
v.04
1. use a sliding window of 4degree sliding at 0.5degree.
v.05
1. interate through all morph files in a given folder.
2. interate 200 times using  new 1800x1800 arrays.
3. generate an excel file to store all coverage data.
v.06
1. use the enriched morphology to resolve a higher resolution of the contour.
v.08
1. improve dend field determination algorithm to avoid dend segments ignored at distal part.
2. enrich function problem: for i in range(0) or for i in range(-1) didnt work.
3. add intersections between dend segment and sector edge. not done.
v.09
1. talked to Stephen. Dont worry too much about tiny bits left outside.
2. more concerned about the sharp dive in parts. should join up in a more smooth contour.
3. try plot all end point first.
4. explore the parameter space and decide to use all points sliding 10 window 30.
5. then use intersection iteration to include all left-out points till none is left out.
6. pack dend field determination process into function.
7. add intersection interation function.
8. add check intersection function.
v.09a
1. intersection iteration in v09 has made improvements(added a few points on contour) but not perfect.
2. adjust intersection iteration function, use all points instead of furthest point to find intersections.
v.09b
1. adjust the point-in-angle condition statement. now using math.atan2(y,x).
2. adjust the intersection slope condition. now using x-in and y-in.
3. adjust the add-up condition. now checking for pre-existing and break loop after.
4. adjust plot feature. now adjusting the width of dend by 0.57528.
v.09d
1. adapted from v09b. attempts in v09c failed.
v.10_output
1. mannual version.
2. add function to output coordinates for cell and dend contour.


"""

################################################################################

import math, xlrd, xlwt, os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

################################################################################

# functions.

def cover(array, ar1, ar2, r):
    covered = []
    for cell in array:
        d = math.sqrt(math.pow(cell[0],2) + math.pow(cell[1],2))
        if not ((d+r<ar1) or (d-r>ar2)):
            covered.append(cell)    
    return covered
    
def calcdis(array):
    minsep = []
    for i in array:
        n = 0
        for j in array:
            if i != j:
                if n == 0:
                    mindis = math.sqrt(math.pow(i[0]-j[0],2) + math.pow(i[1]-j[1],2)) - 2 * r
                    n = 1
                if n == 1:
                    dis = math.sqrt(math.pow(i[0]-j[0],2) + math.pow(i[1]-j[1],2)) - 2 * r
                    if dis < mindis:
                        mindis = dis
        minsep.append(mindis)
    avr = sum(minsep) / len(minsep)
    return avr

def overlapf(covered, DScell):
    overlap = []
    for cell in covered:
        for point in DScell:
            A = point[3] - point[1]
            B = point[0] - point[2]
            C = point[2] * point[1] - point[0] * point[3]
            
            if (math.pow(A,2)+math.pow(B,2)) == 0:
                continue
            
            D = math.sqrt(math.pow(A*cell[0]+B*cell[1]+C,2)/(math.pow(A,2)+math.pow(B,2)))
            PA = math.sqrt(math.pow(cell[0]-point[0],2)+math.pow(cell[1]-point[1],2))
            PB = math.sqrt(math.pow(cell[0]-point[2],2)+math.pow(cell[1]-point[3],2))
            S = math.sqrt(math.pow(point[0]-point[2],2)+math.pow(point[1]-point[3],2))
            if D < r + 0.5*point[4]:
                if ( PA + PB ) <= ( D + math.sqrt(math.pow(D,2)+math.pow(S,2)) ):
                    overlap.append(cell)
                    break
                elif ( PA <= (r+0.5*point[4]) ) or ( PB <= (r+0.5*point[4]) ):
                    overlap.append(cell)
                    break
    return overlap

def dendrolen(DScell, ar1, ar2list):
    dendro = []
    for i in range(10):
        IR = ar1[i]
        THIN = caldendrolen(DScell, ar1[i], ar2list[0][i])
        MIDDLE = caldendrolen(DScell, ar1[i], ar2list[1][i])
        THICK = caldendrolen(DScell, ar1[i], ar2list[2][i])
        dendro.append([IR, THIN, MIDDLE, THICK])
    return dendro

def caldendrolen(DScell, ar1, ar2):
    dendrodis = 0
    for cell in DScell:
        x1 = cell[0]
        y1 = cell[1]
        d1 = math.sqrt(math.pow(x1,2)+math.pow(y1,2))
        x2 = cell[2]
        y2 = cell[3]
        d2 = math.sqrt(math.pow(x2,2)+math.pow(y2,2))
        l = math.sqrt(math.pow(x1-x2,2)+math.pow(y1-y2,2))
        A = math.pow(x2-x1,2)+math.pow(y2-y1,2)
        B = 2*((x2-x1)*(x1)+(y2-y1)*(y1))
        Car1 = math.pow(x1,2)+math.pow(y1,2)-math.pow(ar1,2)
        Car2 = math.pow(x1,2)+math.pow(y1,2)-math.pow(ar2,2)
        
        if l == 0:
            continue
        
        if ar1<=d1<=ar2 and ar1<=d2<=ar2:
            dendrodis += l
        elif (d1<ar1 and d2<ar1) or (d1>ar2 and d2>ar2):
            dendrodis += 0
        elif d1<ar1 and d2>ar2:
            u1ar1 = (-B + math.sqrt(math.pow(B,2)-4*A*Car1))/(2*A)
            u2ar1 = (-B - math.sqrt(math.pow(B,2)-4*A*Car1))/(2*A)
            if 0<u1ar1<1:
                uar1 = u1ar1
            elif 0<u2ar1<1:
                uar1 = u2ar1
            xp1 = x1 + uar1*(x2-x1)
            yp1 = y1 + uar1*(y2-y1)
            u1ar2 = (-B + math.sqrt(math.pow(B,2)-4*A*Car2))/(2*A)
            u2ar2 = (-B - math.sqrt(math.pow(B,2)-4*A*Car2))/(2*A)
            if 0<u1ar2<1:
                uar2 = u1ar2
            elif 0<u2ar2<1:
                uar2 = u2ar2
            xp2 = x1 + uar2*(x2-x1)
            yp2 = y1 + uar2*(y2-y1)
            dpp = math.sqrt(math.pow(xp1-xp2,2)+math.pow(yp1-yp2,2))
            dendrodis += dpp
        elif (d1-ar1)*(d2-ar1)<0:
            u1 = (-B + math.sqrt(math.pow(B,2)-4*A*Car1))/(2*A)
            u2 = (-B - math.sqrt(math.pow(B,2)-4*A*Car1))/(2*A)
            if 0<u1<1:
                u = u1
            elif 0<u2<1:
                u = u2
            xp = x1 + u*(x2-x1)
            yp = y1 + u*(y2-y1)
            if d1>ar1:
                xi = x1
                yi = y1
            elif d2>ar1:
                xi = x2
                yi = y2
            dpi = math.sqrt(math.pow(xp-xi,2)+math.pow(yp-yi,2))
            dendrodis += dpi
        elif (d1-ar2)*(d2-ar2)<0:
            u1 = (-B + math.sqrt(math.pow(B,2)-4*A*Car2))/(2*A)
            u2 = (-B - math.sqrt(math.pow(B,2)-4*A*Car2))/(2*A)
            if 0<u1<1:
                u = u1
            elif 0<u2<1:
                u = u2
            xp = x1 + u*(x2-x1)
            yp = y1 + u*(y2-y1)
            if d1<ar2:
                xi = x1
                yi = y1
            elif d2<ar2:
                xi = x2
                yi = y2
            dpi = math.sqrt(math.pow(xp-xi,2)+math.pow(yp-yi,2))
            dendrodis += dpi
    return dendrodis

def enriching_cell(DScell, resolution):
    
    newcontour = []
    for sec in DScell:
        xs = sec[0]
        ys = sec[1]
        xe = sec[2]
        ye = sec[3]
        dis = math.sqrt((xs-xe)**2+(ys-ye)**2)
        num_point = dis // resolution 
        
        for i in range(int(num_point) - 1):
            xi = xs + i * (1/num_point) * (xe-xs)
            yi = ys + i * (1/num_point) * (ye-ys)
            xj = xs + (i+1) * (1/num_point) * (xe-xs)
            yj = ys + (i+1) * (1/num_point) * (ye-ys)
            
            newcontour.append([xi, yi, xj, yj, sec[4], sec[5]])
                    
    return newcontour

def get_dendfield(DScell):

    sliding = 10
    n_slide = int(360 / sliding + 1)
    anglelist = list(np.linspace(0,2*math.pi,n_slide))
    #print(anglelist)
    dendfield = []
    sliding_range = 15 / 180 * math.pi
    
    anglemaxlist = []
    for i in range(len(anglelist)-1):
        seclist = []
        dislist = []
        low = anglelist[i] - sliding_range
        high = anglelist[i] + sliding_range
        
        # use all pooints
        for cell in DScell:
            x = cell[2]
            y = cell[3]
            dis = math.sqrt(x**2 + y**2)
            
            cell_degree = math.atan2(y,x) / math.pi *180
            if cell_degree < 0:
                cell_degree += 360
            
            low_degree = low / math.pi *180
            high_degree = high / math.pi *180
            
            if low_degree <= cell_degree <= high_degree:
                seclist.append(cell)
                dislist.append(dis)
            

        """
        # use only end points
        for cell in DScell:
            if cell[5] == 'EP':
                x = cell[2]
                y = cell[3]
                dis = math.sqrt(x**2 + y**2)
                sin = y / dis
                cos = x / dis
                
                con1 = (math.sin(low)-sin)*(math.sin(high)-sin)
                con2 = (math.cos(low)-cos)*(math.cos(high)-cos)
                    
                if con1 < 0 and con2 < 0:
                    seclist.append(cell)
                    dislist.append(dis)
        """
        
        for cell in seclist:
            dis = math.sqrt(cell[2]**2 + cell[3]**2)
            if dis == max(dislist) and cell not in anglemaxlist:
                anglemaxlist.append(cell)
                dendfield.append([0,0,cell[2],cell[3],1])
        
    print('the contour contains', len(dendfield), 'points')

    return dendfield, anglemaxlist

def intersection_iteration(DScell, dendfield, anglemaxlist):

    dendfield_iter = dendfield
    dendfield_iter.append(dendfield[0])
    anglemaxlist_iter = anglemaxlist
    anglemaxlist_iter.append(anglemaxlist[0])
    print(len(dendfield_iter))

    outside_checker = 1
    while outside_checker != 0:
        
        dendfield_fixed = []
        dendfield_fixed.append(dendfield_iter[0])
        anglemaxlist_fixed = []
        anglemaxlist_fixed.append(anglemaxlist_iter[0])
            
        for i in range(len(dendfield_iter)-1):
            
            outside_checker = 0
            
            low = dendfield_iter[i]
            high = dendfield_iter[i+1]
            
            #cell_low = anglemaxlist_iter[i]
            cell_high = anglemaxlist_iter[i+1]
            
            xlow = low[2]
            ylow = low[3]
            xhigh = high[2]
            yhigh = high[3]
            
            sector_seclist = []
            #sector_dislist = []
            
            for cell in DScell:
                x = cell[2]
                y = cell[3]
                dis = math.sqrt(x**2 + y**2)
                    
                cell_degree = math.atan2(y,x) / math.pi *180
                if cell_degree < 0:
                    cell_degree += 360
                    
                low_degree = math.atan2(ylow,xlow) / math.pi *180
                high_degree = math.atan2(yhigh,xhigh) / math.pi *180
                
                if low_degree < 0:
                    low_degree += 360
                if high_degree < 0:
                    high_degree += 360
                    
                if low_degree > high_degree:
                    high_degree += 360
                
                if low_degree <= cell_degree <= high_degree:
                    sector_seclist.append(cell)
                    #sector_dislist.append(dis)
                
                intersect_seclist = []
                intersect_dislist = []
            
            for cell in sector_seclist:
                
                xcell = cell[2]
                ycell = cell[3]
                
                xa1 = 0
                ya1 = 0
                xa2 = xcell
                ya2 = ycell
                xb1 = xlow
                yb1 = ylow
                xb2 = xhigh
                yb2 = yhigh
                
                if check_intersection(xa1,ya1,xa2,ya2,xb1,yb1,xb2,yb2):
                    dis = math.sqrt(xcell**2 + ycell**2)
                    intersect_seclist.append(cell)
                    intersect_dislist.append(dis)
            
            if len(intersect_seclist) != 0:
            
                for cell in intersect_seclist:
                    dis = math.sqrt(cell[2]**2 + cell[3]**2)
                    if dis == max(intersect_dislist) and cell not in anglemaxlist_fixed:
                        dendfield_fixed.append([0,0,cell[2],cell[3],1])
                        anglemaxlist_fixed.append(cell)
                        dendfield_fixed.append(high)
                        anglemaxlist_fixed.append(cell_high)
                        outside_checker = 1
                        break
                        
            else:
                dendfield_fixed.append(high)
                anglemaxlist_fixed.append(cell_high)
                
        dendfield_iter = dendfield_fixed
        anglemaxlist_iter = anglemaxlist_fixed
        print(len(dendfield_iter))
        
    del dendfield_iter[-1]
    del anglemaxlist_iter[-1]
    
    return dendfield_iter, anglemaxlist_iter

def check_intersection(xa1,ya1,xa2,ya2,xb1,yb1,xb2,yb2):
    
    if xa1 != xa2 and xb1 != xb2:
        ka = (ya2 - ya1) / (xa2 - xa1)
        ba = ya1 - ka * xa1
        kb = (yb2 - yb1) / (xb2 - xb1)
        bb = yb1 - kb * xb1
        
        if ka != kb:
            
            xinter = - (ba - bb) / (ka - kb)
            #yinter = ka * xinter - ba
            
            if (xinter-xa1)*(xinter-xa2)<0 and (xinter-xb1)*(xinter-xb2)<0:
                truth = True
            else:
                truth = False
        else:
            truth = False
    else:
        truth = False
        
    return truth

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

filepath = 'C:/Users/uqxzho15/Desktop/py code/bipolar sim 3 sustained off alpha/raw data/data/concentric/'
#filepath = 'C:/Users/uqxzho15/Desktop/py code/bipolar sim 3 sustained off alpha/raw data/data/test/'
#filepath = 'C:/Users/uqxzho15/Desktop/py code/bipolar sim 3 sustained off alpha/raw data/unrotated test'
#filepath = 'C:/Users/uqxzho15/Desktop/py code/bipolar sim 3 sustained off alpha/raw data/newly added 210212'
#filepath = 'C:/Users/uqxzho15/Desktop/py code/bipolar sim 3 sustained off alpha/raw data/newly added 210215'

filelist = []
for file in os.listdir(filepath):
    file_path = os.path.join(filepath, file)
    if os.path.splitext(file_path)[1] == '.xlsx':
        filelist.append(file_path)

# input ganglion cell and analyse.

data = []
avg_data = []

for inputfile in filelist:
    
    # input ganglion cell.
    
    cellname = ''
    for i in range(len(inputfile)):
        if inputfile[i]=='c' and inputfile[i+1]=='e' and inputfile[i+2]=='l' and inputfile[i+3]=='l':
            cellname = inputfile[i-9:i-1] + ' cell ' + inputfile[i+5:i+7]
    print('now running cell: ', cellname)
    
    wb = xlrd.open_workbook(inputfile)
    
    for n in range(len(wb.sheet_names())):
        if wb.sheet_names()[n] == 'Segment Points - Dendrites':
            s = wb.sheets()[n]
            line = []
            for row in range(s.nrows):
                col_value = []
                for col in range(s.ncols):
                    col_value.append(s.cell(row,col).value)
                line.append(col_value)

    DScell = []
    n = 0
    for cell in line:
        if n == 0:
            n = 1
        else:
            DScell.append((float(cell[2]), float(cell[3]), float(cell[5]), float(cell[6]), float(cell[13]), str(cell[8])))
    # cell(startX, startY, endX, endY, averaged diameter, point type)
    # cell[   0       1      2     3           4               5    ]


    #DScell_enriched = enriching_cell(DScell, 0.5)

    ################################################################################

    dendfield = []
    anglemaxlist = []
    dendfield, anglemaxlist = get_dendfield(DScell)

    ################################################################################
    
    dendfield_trans = dendfield
    anglemaxlist_trans = anglemaxlist
    
    dendfield, anglemaxlist = intersection_iteration(DScell, dendfield_trans, anglemaxlist_trans)
    
    ################################################################################

    n_plot = 0
    n_rep = 0
    data_line = []
    data_line.append(cellname)
    
    for bipolar in bipolarselect:
        
        dendfield_within = overlapf(bipolar, dendfield)
        dendfield_connected = overlapf(dendfield_within, DScell)
        
        coverage = len(dendfield_connected) / len(dendfield_within)
        data_line.append(coverage)
        #print('coverage is:', coverage)
    
        if n_rep % 10 == 0:
            print(n_rep + 1, '/', rep, 'done...')
        n_rep += 1

    ################################################################################
    
        if n_plot == 0:
            n_plot = 1
    
            fig, ax = plt.subplots()
            ax = plt.gca(aspect='equal')
            ax.cla()  
            ax.set_xlim((-900, 900))
            ax.set_ylim((-900, 900))
            
            for cell in DScell:
                line1 = [(cell[0], cell[1]), (cell[2], cell[3])]
                (line1_xs, line1_ys) = zip(*line1)
                line = Line2D(line1_xs, line1_ys, linewidth=cell[4]/7.27*1.5*0.57528, color='black')
                line.set_zorder(0)
                ax.add_line(line)
            
            for i in range(len(anglemaxlist)-1):
                line1 = [(anglemaxlist[i][2], anglemaxlist[i][3]), (anglemaxlist[i+1][2], anglemaxlist[i+1][3])]
                (line1_xs, line1_ys) = zip(*line1)
                line = Line2D(line1_xs, line1_ys, linewidth=0.3, color='red')
                line.set_zorder(0)
                ax.add_line(line)
        
            line1 = [(anglemaxlist[-1][2], anglemaxlist[-1][3]), (anglemaxlist[0][2], anglemaxlist[0][3])]
            (line1_xs, line1_ys) = zip(*line1)
            line = Line2D(line1_xs, line1_ys, linewidth=0.3, color='red')
            line.set_zorder(0)
            ax.add_line(line)
            
            ax.axis('off')
            figname = cellname + ' 1 dend field contour 27.svg'
            fig.savefig(figname, format = 'svg', dpi = 1200)
        
        ################################################################################
            """        
            fig, ax = plt.subplots()
            ax = plt.gca(aspect='equal')
            ax.cla()  
            ax.set_xlim((-900, 900))
            ax.set_ylim((-900, 900))
            
            output_list = []
            for cell in bipolar:
                output_list.append(plt.Circle((cell[0], cell[1]), r, color='grey', fill=False, linewidth=0.1))
            for cell in dendfield_within:
                output_list.append(plt.Circle((cell[0], cell[1]), r, facecolor='red', edgecolor='none', alpha=0.7))
            
            for cell in output_list:    
               ax.add_artist(cell)
            
            for cell in DScell:
                line1 = [(cell[0], cell[1]), (cell[2], cell[3])]
                (line1_xs, line1_ys) = zip(*line1)
                line = Line2D(line1_xs, line1_ys, linewidth=cell[4]/7.27*1.5*0.57528, color='black')
                line.set_zorder(0)
                ax.add_line(line)
            
            ax.axis('off')
            figname = cellname + ' 2 bps within dend field.svg'
            fig.savefig(figname, format = 'svg', dpi = 1200)
            
        ################################################################################
           
            fig, ax = plt.subplots()
            ax = plt.gca(aspect='equal')
            ax.cla()  
            ax.set_xlim((-900, 900))
            ax.set_ylim((-900, 900))
            
            output_list = []
            for cell in dendfield_within:
                if cell not in dendfield_connected:
                    output_list.append(plt.Circle((cell[0], cell[1]), r, facecolor='green', edgecolor='none', alpha=0.7))
            for cell in dendfield_connected:
                output_list.append(plt.Circle((cell[0], cell[1]), r, facecolor='red', edgecolor='none', alpha=0.7))
            
            for cell in output_list:    
               ax.add_artist(cell)
            
            for cell in DScell:
                line1 = [(cell[0], cell[1]), (cell[2], cell[3])]
                (line1_xs, line1_ys) = zip(*line1)
                line = Line2D(line1_xs, line1_ys, linewidth=cell[4]/7.27*1.5*0.57528, color='black')
                line.set_zorder(0)
                ax.add_line(line)
            
            ax.axis('off')
            figname = cellname + ' 3 bps covered and not covered.svg'
            fig.savefig(figname, format = 'svg', dpi = 1200)

        ################################################################################
        
            fig, ax = plt.subplots()
            ax = plt.gca(aspect='equal')
            ax.cla()  
            ax.set_xlim((-900, 900))
            ax.set_ylim((-900, 900))
            
            output_list = []
            for cell in bipolar:
                output_list.append(plt.Circle((cell[0], cell[1]), r, color='grey', fill=False, linewidth=0.1))
            for cell in dendfield_within:
                if cell not in dendfield_connected:
                    output_list.append(plt.Circle((cell[0], cell[1]), r, facecolor='green', edgecolor='none', alpha=0.7))
            for cell in dendfield_connected:
                output_list.append(plt.Circle((cell[0], cell[1]), r, facecolor='red', edgecolor='none', alpha=0.7))
            
            for cell in output_list:    
               ax.add_artist(cell)
            
            for cell in DScell:
                line1 = [(cell[0], cell[1]), (cell[2], cell[3])]
                (line1_xs, line1_ys) = zip(*line1)
                line = Line2D(line1_xs, line1_ys, linewidth=cell[4]/7.27*1.5*0.57528, color='black')
                line.set_zorder(0)
                ax.add_line(line)

            for i in range(len(anglemaxlist)-1):
                line1 = [(anglemaxlist[i][2], anglemaxlist[i][3]), (anglemaxlist[i+1][2], anglemaxlist[i+1][3])]
                (line1_xs, line1_ys) = zip(*line1)
                line = Line2D(line1_xs, line1_ys, linewidth=0.5, color='blue')
                line.set_zorder(0)
                ax.add_line(line)
        
            line1 = [(anglemaxlist[-1][2], anglemaxlist[-1][3]), (anglemaxlist[0][2], anglemaxlist[0][3])]
            (line1_xs, line1_ys) = zip(*line1)
            line = Line2D(line1_xs, line1_ys, linewidth=0.5, color='blue')
            line.set_zorder(0)
            ax.add_line(line)
            
            ax.axis('off')
            figname = cellname + ' 4 full elements.svg'
            fig.savefig(figname, format = 'svg', dpi = 1200)
            """        
    ################################################################################
        
    print(cellname, 'done.\n')
    data.append(data_line)
    avg_data.append([cellname, np.mean(data_line[1:])])

data.append([' '])
data.append(['avg data'])
for line in avg_data:
    data.append(line)

book = xlwt.Workbook(encoding="utf-8")
sheet1 = book.add_sheet("Sheet 1")
for p in range(len(data)):
    for q in range(len(data[p])):
        sheet1.write(p,q,data[p][q])
book.save('Coverage Summary.xls')
    
    
    
    