"""
bipolar array generating algorithm.

"""

#import matplotlib.pyplot as plt
import random, math, xlwt

def formarray(l, r, mode, tess, regmindis, pct, area, deltas):
    array = []
    if mode == 1:
        array = regulararray(l, r, tess, regmindis)
    elif mode == 2:
        array = randomarray(l, r, pct, area, deltas)
    return array

def regulararray(l, r, tess, regmindis):
    regarray = []
    rr = r + 0.5*regmindis
    if tess == 1:
        rangeX = (-0.5*l, -0.5*l+2*rr)
        rangeY = (0.5*l-2*rr, 0.5*l)
        x0 = random.randrange(*rangeX)
        y0 = random.randrange(*rangeY)
        x = x0
        y = y0
        while y > -0.5*l:
            while x < 0.5*l:
                regarray.append((x,y))
                x += 2*rr
            y -= 2*rr
            x = x0
    elif tess == 2:
        rangeX = (-0.5*l, -0.5*l+2*rr)
        rangeY = (int(0.5*l-math.sqrt(3)*rr), 0.5*l)
        x0 = random.randrange(*rangeX)
        y0 = random.randrange(*rangeY)
        x = x0
        y = y0
        i = 0
        while y > -0.5*l:
            while x < 0.5*l:
                if -0.5*l < x:
                    regarray.append((x,y))
                x += 2*rr
            i += 1
            y -= math.sqrt(3)*rr
            x = x0 - i * rr
    return regarray

def randomarray(l, r, pct, area, deltas):
    randarray = []
    num = int((pct * math.pow(l,2)) / (3.1415 * math.pow(r,2)))
    i = 0
    exclude = set()
    while i < num:
        [(x,y)] = random.sample(area,1)
        randarray.append((x,y))
        exclude = set((x+dx, y+dy) for (dx,dy) in deltas)
        area = area - exclude
        i += 1
        if len(area) == 0:
            break
    return randarray

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

# part 1. specify all the parameters.

l = 1800          # area size is l*l.
r = 20            # bipolar axon arbor radius.
mode = 2          # 1 for regular, 2 for random.
tess = 2          # 1 for sqaure, 2 for regular hexagon.
regmindis = 10    # the minimum distance between neighbouring axon arbor in regular mode.
pct = 1.00        # percentage of area that bipolar axon arbors occupy.
rep = 250         # how many repetitions for each annulus setting.
ovd = 3
bipolar_writing_structure = 7

'''
pctlist =   [0.1,   0.2,    0.3,    0.4,    0.5,    1.0,    1.0,    1.0,    1.0,    1.0]
ovdlist =   [0,     0,      0,      0,      0,      1,      3,      5,      7,      9]

ar1 =       [0.0,     50.0,     75.0,     150.0,    225.0,    300.0,    375.0,    450.0,    525.0,    600.0]
ar2thin =   [125.0,   134.5,    146.0,    195.0,    257.5,    325.0,    395.5,    467.0,    539.5,    613.0]
ar2middle = [146.0,   154.5,    164.5,    209.5,    268.5,    333.5,    402.5,    473.0,    545.0,    617.5]
ar2thick =  [218.0,   223.5,    230.5,    264.5,    313.0,    371.0,    433.5,    500.0,    568.5,    638.5]
ar2list =   [ar2thin, ar2middle, ar2thick]
'''

area = set()
for x in range(int(-0.5*l), int(0.5*l+1)):
    for y in range(int(-0.5*l), int(0.5*l+1)):
        area.add((x,y))

deltas = set()
for x in range(-2*r-ovd, 2*r-ovd+1):
    for y in range(-2*r-ovd, 2*r-ovd+1):
        if math.pow(x,2) + math.pow(y,2) <= math.pow(2*r-ovd,2):
            deltas.add((x,y))

sheet = []
for i in range(rep):
    print('----------', i+1)
    bipolar = formarray(l, r, mode, tess, regmindis, pct, area, deltas)
    print(len(bipolar))
    
    sheet.append(['array','l','r','mode','pct','ovd','avrmindis','cov', 'dens/mm2'])
    avrmindis = calcdis(bipolar)
    cov = len(bipolar)*3.1415*math.pow(r,2)/math.pow(l,2)
    dens = len(bipolar)/math.pow(l,2)*math.pow(1000,2)
    sheet.append([i+1, l, r, mode, pct, ovd, avrmindis, cov, dens])
    
    xline = [[] for t in range(bipolar_writing_structure)]
    yline = [[] for t in range(bipolar_writing_structure)]
    
    for t in range(len(xline)):
        for n in range(250):
            if (t*250+n) < len(bipolar):
                xline[t].append(bipolar[t*250+n][0])
                yline[t].append(bipolar[t*250+n][1])
    
    for x in range(len(xline)):
        sheet.append(xline[x])
        sheet.append(yline[x])
    
    print('++++++++++', i+1)

outputfile = '200 arrays 1800.xls'
book = xlwt.Workbook(encoding="utf-8")
sheet1 = book.add_sheet("Sheet 1")
for p in range(len(sheet)):
    for q in range(len(sheet[p])):
        sheet1.write(p,q,sheet[p][q])
book.save(outputfile)
print('\ndone\n')
