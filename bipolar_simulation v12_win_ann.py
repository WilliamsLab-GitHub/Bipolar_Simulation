"""
Bipolar cell axon arbor distribution simulation
v.12_win_ann
1. bug in line 137, function caldendrolen fixed. for ar2 condition, where should be an 'and' was an 'or'.
2. ar2 thin 195.0 should be 195.5. fixed.
3. change the conditions in caldendrolen function for better accuracy. now dpp condition first.

v.11_win_ann
1. modified for reading 1000 previously generated bipolar arrays
2. modified for batch anaylsis. now can automatically analyse all files in one folder.
3. delete module for calculating normalised results.
4. restructure output to put results on top.
5. bug fixed.
6. At least 10 times faster and reproducible. 3.6min for one array before, now 12s.

v.10_win_ann_cov test spe
1. fix the rounding problem for the annulus radius matrix.
2. put area, deltas and ovd in the main part to save time.
3. add cell density.
4. automatically go through all pct-ovd combinations.
5. overlap is now a function to deploy.
6. be able to recognise cell name and automatically name the output excel file.
7. define new functions named dendrolen and caldendrolen, to calculate dendritic length in annulus.
8. bug fixed.
9. structure re-organised.
10. this version is specialised in coverage test.
11. bug prevented in overlapf and caldendrodis when 2 end points of a segment overlap in x and y.
    (this problem occurs specially in 3D tracing.)

v.09_win_ann
1. modified for reading x1,y1,x2,y2,thickness of dendritic segments.
2. overlap algorithm updated.
3. bug fixed.
4. change area into numbers in result.

v.08_win_ann
1. modified for running in windows environment.
2. modified for annulus simmulation.
3. bug in function 'cover' fixed. now parameter 'r' has been included in the input variables.
4. script modified for reading excel files to fit in more format.
   cos Arne gave me files in all different shapes.
5. output file name format modified to include more parameter setting.

v.07a
1. modified for running in windows environment.

v.07
1. coupled with the morphology.

v.06a
1. automatically calculate the the average area for each simulation.
2. now display coverage for each simulation.
3. skip plotting to save time.
4. bug fixed.
5. modify line 100 for generating random array, saving 1/3 of the processing time.
"""

#import matplotlib.pyplot as plt
import math, xlrd, xlwt, os

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


# main body.

# specify patameters.

l = 1400          # area size is l*l.
r = 20            # bipolar axon arbor radius.
mode = 2          # 1 for regular, 2 for random.
tess = 2          # 1 for sqaure, 2 for regular hexagon.
regmindis = 10    # the minimum distance between neighbouring axon arbor in regular mode.
pct = 1.00        # percentage of area that bipolar axon arbors occupy.
rep = 200        # how many repetitions for each annulus setting. 0< rep <= 1000.
ovd = 3

ar1 =       [0.0,     25.0,     37.5,     75.0,    112.5,    150.0,    187.5,    225.0,    262.5,    300.0]
ar2thin =   [31.250,   46.157,    55.193,    86.314,    120.598,    156.250,    192.571,    229.260,    266.169,    303.221]
ar2middle = [31.250,   46.157,    55.193,    86.314,    120.598,    156.250,    192.571,    229.260,    266.169,    303.221]
ar2thick =  [31.250,   46.157,    55.193,    86.314,    120.598,    156.250,    192.571,    229.260,    266.169,    303.221]
ar2list =   [ar2thin, ar2middle, ar2thick]

# read previously generated bipolar arrays.

bipolarfile = "'C:/Users/uqxzho15/Desktop/py code/bipolar sim 3 sustained off alpha/500 arrays.xls'"
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
        for n in range(5):
            for t in range(len(bpline[i+2+2*n])):
                if bpline[i+2+2*n][t] != '':
                    bipolar.append((bpline[i+2+2*n][t], bpline[i+2+2*n+1][t]))
        bipolarlist.append(bipolar)
bipolarselect = bipolarlist[0:rep]

# read ganglion cell list from designated folder.

filepath = 'C:/Users/uqxzho15/Desktop/py code/bipolar sim 3 sustained off alpha/raw data/data/batch 2/'

filelist = []
for file in os.listdir(filepath):
    file_path = os.path.join(filepath, file)
    if os.path.splitext(file_path)[1] == '.xlsx':
        filelist.append(file_path)

# input ganglion cell and analyse.

for inputfile in filelist:
    
    # input ganglion cell.
    
    cellname = ''
    for i in range(len(inputfile)):
        if inputfile[i]=='c' and inputfile[i+1]=='e' and inputfile[i+2]=='l' and inputfile[i+3]=='l':
            cellname = inputfile[i-9:i-1] + ' c' + inputfile[i+5:i+7]
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
            for n in range(len(cell)):
                if cell[n] == 'Start X':
                    n0 = n
                    break
            n = 1
        else:
            DScell.append((float(cell[n0]), float(cell[n0+1]), float(cell[n0+3]), float(cell[n0+4]), float(cell[n0+11])))

    # form headings of the output sheet.

    sheet = []
    if mode == 1:
        sheet.append(['l','r','mode','tess','regmindis','coverage'])
        sheet.append([l,r,mode,tess])
    elif mode == 2:
        sheet.append(['l','r','mode','pct','avrmindis','coverage', 'dens/mm2'])
        sheet.append([l,r,mode,pct])
    sheet.append([' '])
    
    arrayline = ['array',' ']
    for i in range(rep):
        arrayline.append(i+1)
            
    avrmindisline = ['avrmindis',' ']    
    coverageline = ['coverage',' ']    
    densityline = ['dens/mm2',' ']
    
    data = [[] for a in range(30)]
    datan = 0
    for ar2 in ar2list:
        for i in range(len(ar1)):
            data[datan].append(ar1[i])
            data[datan].append(ar2[i])
            datan += 1
    
    # analyse all bipolar arrays.
    
    t = 0
    print('total times of simulation: ', rep)
    for bipolar in bipolarselect:
        t += 1
        print('----------', t)

        if mode == 1:
            avrmindisline.append(regmindis)
        elif mode == 2:
            avrmindisline.append(calcdis(bipolar))
        
        coverageline.append(len(bipolar)*3.1415*math.pow(r,2)/math.pow(l,2))
        densityline.append(len(bipolar)/math.pow(l,2)*1000*1000)
               
        datan = 0
        for ar2 in ar2list:
            for i in range(len(ar1)):
                covered = cover(bipolar, ar1[i], ar2[i], r)
                overlap = overlapf(covered, DScell)
                
                totalnumber = len(overlap)#*3.1415*math.pow(r,2)
                data[datan].append(totalnumber)
                datan += 1
    
    # summarise results.
    
    result = []
    for i in range(10):
        IR = ar1[i]
        THIN = sum(data[i][2:])/len(data[i][2:])
        MIDDLE = sum(data[i+10][2:])/len(data[i+10][2:])
        THICK = sum(data[i+20][2:])/len(data[i+20][2:])
        result.append([IR, THIN, MIDDLE, THICK])
    
    sheet[1].append(sum(avrmindisline[2:])/len(avrmindisline[2:]))
    sheet[1].append(sum(coverageline[2:])/len(coverageline[2:]))
    sheet[1].append(sum(densityline[2:])/len(densityline[2:]))
    
    sheet.append(['Result'])
    sheet.append(['IR','Thin','Middle','Thick'])
    for line in result:
        sheet.append(line)
    sheet.append([' '])
    
    sheet.append(['Dendro lenth'])
    sheet.append(['IR','Thin','Middle','Thick'])
    dendro = dendrolen(DScell, ar1, ar2list)
    for line in dendro:
        sheet.append(line)
    sheet.append([' '])
    
    section = [2, 252, 502, 752]
    for x in section:
        arraysection = arrayline[0:2] + arrayline[x:x+250]
        avrmindissection = avrmindisline[0:2] + avrmindisline[x:x+250]
        coveragesection = coverageline[0:2] + coverageline[x:x+250]
        densitysection = densityline[0:2] + densityline[x:x+250]
        sheet.append(arraysection)
        sheet.append(avrmindissection)
        sheet.append(coveragesection)
        sheet.append(densitysection)
        datan = 0
        for line in data:
            datasection = line[0:2] + line[x:x+250]
            datan += 1
            if datan == 1:
                sheet.append(['thin'])
                sheet.append(['IR','OR'])
            elif datan == 11:
                sheet.append(['middle'])
                sheet.append(['IR','OR'])
            elif datan == 21:
                sheet.append(['thick'])
                sheet.append(['IR','OR'])
            sheet.append(datasection)
        sheet.append(' ')
    
    finalcov = sum(coverageline[2:])/len(coverageline[2:])
    outputfile = 'bp sim '+cellname+' pct '+str(pct)+' ovd '+str(ovd)+' cov ' +str(round(finalcov,5))+'.xls'
    book = xlwt.Workbook(encoding="utf-8")
    sheet1 = book.add_sheet("Sheet 1")
    for p in range(len(sheet)):
        for q in range(len(sheet[p])):
                sheet1.write(p,q,sheet[p][q])
    book.save(outputfile)
    print('done for ', cellname)

print('\ndone\n')
            
