### Importing necessary modules ###

import numpy as np
import glob
import os


### Initialising lists to be appended and info to find file path ###

listofpos_list =[]       #a list to contain the lists of body information read from text files
cwd = os.getcwd()        #gets the current working directory to know where to find planetary information
label_list = []          #a list of the name of the bodies


### Opening and scraping every Nasa Ephemeris results for each body ### 

for file in glob.glob(str(cwd)+'/PlanetsUniq/*.txt'):   #opens all nasa body information within folder
    with open(file, 'r') as file2:                      #opens and reads information within file
        filedata = file2.read()

    #finds the mass data in the file
    mass = filedata.find('kg')             #finds the position of kg symbol to locate its value after
    massend = filedata[mass:].find('$')    #finds the end of the mass value 
    mass1 = float(filedata[mass:massend+mass].translate(None, 'kg )~=-+'))  #removes all symbols which aren't numbers to isolate mass value and convert to float

    massorder1 = filedata.find('Mass, 10^')       #finds order of magniutde of mass measurement
    massorder2 = filedata.find('Mass (10^')       #accounts for variations between the ways NASA presents mass data
    massorder3 = filedata.find('Mass Pluto (10^') #accounts for the specific pluto variation :(
    massorder = 0                                 #accounts for if there is no mass order
    
    if massorder1 >=  0:
        massorder = float(filedata[massorder1+len('Mass, 10^'):massorder1+len('Mass, 10^')+2]) #converts order of magnitude to usual float for each variation
    if massorder2 >= 0:
        massorder = float(filedata[massorder2+len('Mass (10^'):massorder2+len('Mass (10^')+2])
    if massorder3 >= 0:
        massorder = float(filedata[massorder3+len('Mass Pluto (10^'):massorder3+len('Mass Pluto (10^')+2])

    mass2 = mass1*10.0**float(massorder)                                        #gives the mass of the body in required format (float)
    
    #finds the position data in the file
    position = filedata.find('$$SOE')                                            #unique string which appears before coordinate data 
    position = filedata[position+6: position + filedata[position+6:].find('\n')] #finds the end of the line => coordinates at only one time
    pos_list = position.split(',')                                               #creates a list from the raw nas coordinate string
    del pos_list[0], pos_list[0],pos_list[-1],pos_list[-1],pos_list[-1]          #deletes all elements in list that aren't x,y,z and vx,vy,vz coordinates

    pos_list.insert(0, mass2)                         #inserts the mass value into the position list at the start
    listofpos_list.append(pos_list)                   #adds position list to the list of position lists
    #finds the name of the body and adds it to a list
    body_name = str(file)                             #starts by taking name as all file path
    body_name = body_name[body_name.find("PlanetsUniq/")+12:body_name.find(".txt")] #reduces name to file path between last dir name and .txt => File name
    label_list.append(body_name)                      #adds scraped name to list 

#puts the data from the files into arrays
bodies = np.array(listofpos_list)                       #array of mass followed by coordinantes for every planet
labels = np.array(label_list)                           #array of body names
bodiesuncorrected = np.insert(bodies, 1,labels,axis=1)  #joins the body names to start of other array in style required to be read for Particle3D objects


#### CENTER OF MASS CORRECTIONS ####

Num_Bodies = len(bodies[:,0].tolist())    #Calculates an integer number of bodies in system
mass_list= []                             #creates empty list for body masses
momentum = np.zeros((int(Num_Bodies),3))  #creates empty array for momentums
v_cor = np.zeros((int(Num_Bodies),3))     #creates empty array for corrected velocities

for i in range(Num_Bodies):                 #loop to calculate momentum for every body
    mass = float(bodiesuncorrected[i,0])    #finds mass of specific body
    mass_list.append(mass)                  #creates list of body masses
    v_init = np.array(bodiesuncorrected[i,5:], dtype=np.float64)  #creates array of velocities for each body in different row
    momentum[i] = v_init*mass               #calculates momentum for each individual body

Ptot = np.sum(momentum, dtype = np.float64, axis = 0)  #Total momentum of all bodies
mtot = sum(mass_list)                                  #Total mass of all bodies
v_CoM = Ptot/mtot                                      #The center of mass velocity array

for j in range(Num_Bodies):                                  #loop to correct each body velocity relative to CoM velocity
    vj = np.array(bodiesuncorrected[j,5:], dtype=np.float64) #creates array of each body's velocity coordinantes
    v_cor[j] = np.subtract(vj,v_CoM)                         #corrects each body's velocity by subtracting CoM velocity

bodiesuncorrecteddel = np.delete(bodiesuncorrected, [5,6,7], axis=1)  #deletes old body velocity coordinates
bodiesfinal = np.concatenate((bodiesuncorrecteddel,v_cor), axis=1)    #inserts the correcred velocities where the old velocities were
