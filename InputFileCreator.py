import numpy as np
import csv
import glob
import re
import os

listofposlist =[]
cwd = os.getcwd()
labellist = []

for file in glob.glob(str(cwd)+'/PlanetsUniq/*.txt'):
    linepos = 0
    with open(file, 'r') as file2:
        filedata = file2.read()

    mass = filedata.find('kg')
    massend = filedata[mass:].find('$')
    mass1 = float(filedata[mass:massend+mass].translate(None, 'kg )~=-+'))

    massorder1 = filedata.find('Mass, 10^')
    massorder2 = filedata.find('Mass (10^')
    massorder3 = filedata.find('Mass Pluto (10^')
    massorder = 0

    if massorder1 >=  0:
        massorder = float(filedata[massorder1+len('Mass, 10^'):massorder1+len('Mass, 10^')+2])
    if massorder2 >= 0:
        massorder = float(filedata[massorder2+len('Mass (10^'):massorder2+len('Mass (10^')+2])
    if massorder3 >= 0:
        massorder = float(filedata[massorder3+len('Mass Pluto (10^'):massorder3+len('Mass Pluto (10^')+2])

    mass2 = mass1*10.0**float(massorder)
    
    pos = filedata.find('$$SOE')
    pos = filedata[pos+6: pos + filedata[pos+6:].find('\n')]
    poslist = pos.split(',')
    del poslist[0], poslist[0],poslist[-1],poslist[-1],poslist[-1]

    poslist.insert(0, mass2)
    listofposlist.append(poslist)
    name = str(file)
    name = name[name.find("PlanetsUniq/")+12:name.find(".txt")]
    labellist.append(name)


bodies = np.array(listofposlist)
labels = np.array(labellist)
bodiesuncorrected = np.insert(bodies, 1,labels,axis=1)

#### CENTER OF MASS CORRECTIONS ####

Num_Bodies = len(bodies[:,0].tolist())
ml = []
P = np.zeros((int(Num_Bodies),3))
v_cor = np.zeros((int(Num_Bodies),3))

for i in range(Num_Bodies):
    m = float(bodiesuncorrected[i,0])
    ml.append(m)
    vi = np.array(bodiesuncorrected[i,5:], dtype=np.float64)
    P[i] = vi*m

Ptot = np.sum(P, dtype = np.float64, axis = 0)
mtot = sum(ml)
v_CoM = Ptot/mtot

for j in range(Num_Bodies):
    vj = np.array(bodiesuncorrected[j,5:], dtype=np.float64)
    v_cor[j] = np.subtract(vj,v_CoM)

bodiesuncorrecteddel = np.delete(bodiesuncorrected, [5,6,7], axis=1)
bodiesfinal = np.concatenate((bodiesuncorrecteddel,v_cor), axis=1)
