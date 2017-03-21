### Importing neccesary modules ###
import math as m
import numpy as np
import Vector as vctr
import sys
from Particle3D import Particle3D as p3D
import InputFileCreator as I
import datetime
import csv
import matplotlib.pyplot as plt


### Creating variables, arrays and files to write ###

#Getting all bodies used in simulation from input file creator
bodies = I.bodies_final

#Opening files
VMDFILE = open('VMDtrajectory.xyz', 'w')
Orbit_info = open('Orbital Information.txt', 'w')

#Reading initial conditions from CSV file dictionary
input_file = csv.DictReader(open("SimulationParameters.csv"))
for row in input_file:
    dt = int(row["TimeStep"])
    numstep = int(row["NumberSteps"])
    vmdstep = int(row["VMDSteps"])

#getting number of bodies - useful to know how many bodies to loop through in for loop
Num_Bodies = len(bodies[:,0].tolist())

#starting time at 0
time = 0
f=0.0
u=0.0
w=0.0
f_new = 0.0
u_new =0.0
period = 0.0
kinetic_energy =0.0
total_kinetic= 0.0
potential_energy =0.0
total_potential =0.0

#Getting date for progress bar in NICE format

force_array = np.zeros((int(Num_Bodies),3))
new_force_array = np.zeros((int(Num_Bodies),3))

###Defining Periapse and Apoapsis###

Periapse = np.zeros((int(Num_Bodies),1))

Theta = np.zeros((int(Num_Bodies),1))
Period = []

Planet_list = []
#make initial periapse abitrarily large so that the sunsep and moonsep will always be less to start
for x1 in xrange(Num_Bodies):
    particle0 = p3D.from_file(bodies[x1].tolist())
    Planet_list.append(particle0)

for x3 in xrange(Num_Bodies):
    particle0 = Planet_list[x3]
    if particle0.label == "SUN":
        sunnumber = x3
        pSun = particle0
        pSun_vel = pSun.velocity

for x4 in xrange(Num_Bodies):
    particle0 = Planet_list[x4]
    Period.append([])
    if x4 != sunnumber:
        Periapse[x4] = p3D.mag_sep(particle0,pSun)
    #moon
    if particle0.label == "MOON":
        moonnumber = x4
        pMoon = Planet_list[x4]
    if particle0.label == "EARTH":
        earthnumber = x4
        pEarth = Planet_list[x4]


PeriapseMoon = [p3D.mag_sep(pMoon,pEarth)]

Apoapsis = np.zeros((int(Num_Bodies),1))
ApoapseMoon = [0.0]
ThetaMoon = 0.0
PeriodMoon = []
Energy_list = []
time_list = []

###Periapse and Apoapsis are initialised###

###All variables etc required are created###


###Running through simulation for the number of steps and timestep size defined in simulation parameters###

for i in range(1,numstep):
    if i%vmdstep ==0:
        VMDFILE.write('{}\n'.format(Num_Bodies))
        VMDFILE.write('Point = {}\n'.format(i/vmdstep))
        for x1 in range(Num_Bodies):
            particle0 = Planet_list[x1]
            VMDFILE.write('{}'.format(particle0))

    if i%(numstep/10) ==0:
        print "As of " +str(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")) +" --> The program is " +str(i*100/numstep) +"% completed."

    for x1 in xrange(Num_Bodies):
        #Initialising a "main" planet
        particle0 = Planet_list[x1]
    
        ###Periapse and Apoapse Bollocks### 

        #loop to find moon periapse and apoapse + period
        if x1 == earthnumber:
            MoonSep = p3D.mag_sep(pMoon,pEarth)
            if PeriapseMoon[0] > MoonSep:
                PeriapseMoon[0] = MoonSep
            if ApoapseMoon[0] < MoonSep:
                ApoapseMoon[0] = MoonSep
            moon_ang_vel = pMoon.ang_vel(pMoon,pEarth,w)
            ThetaMoon = moon_ang_vel*dt+ThetaMoon
            if ThetaMoon > 2*m.pi:
                PeriodMoon.append(i*dt/(60.0*60.0*24.0) - sum(PeriodMoon) )
                ThetaMoon = ThetaMoon-2*m.pi

        #loop to find other planet periapse and apoapse + period
        if x1 != moonnumber and x1 != sunnumber:
            Theta[x1] = particle0.ang_vel(particle0,pSun,w)*dt+Theta[x1]
            if Theta[x1] > 2*m.pi:
                Period[x1].append(i*dt/(60.0*60.0*24.0)-sum(Period[x1]))
                Theta[x1] = Theta[x1]-2*m.pi
            if i == numstep-1 and Period[x1]==[]:
                Period[x1].append(i*dt*2*m.pi/(Theta[x1]*60.0*60.0*24.0))

            SunSep = p3D.mag_sep(particle0,pSun)
            if Periapse[x1] > SunSep:
                Periapse[x1] = SunSep
            if Apoapsis[x1] < SunSep:
                Apoapsis[x1] = SunSep

        ###End of Periapse and Apoapse Bollocks###

        kinetic_energy = (1.0/2.0)*particle0.mass*(vctr.SqMag(particle0.velocity))
        total_kinetic = total_kinetic+kinetic_energy

        ###loop to calculate all forces acting on main planet [x1] from secondary planets [x2]###
        for x2 in xrange(Num_Bodies):
            if x1 < x2:
                #initialising a "secondary" planet#
                particle1 = Planet_list[x2]
                
                #Calculating force from between secondary planet to main planet
                interaction_force = particle0.force(particle0,particle1,f)
                
                potential_energy = particle0.poten(particle0,particle1,u)
#                if particle1.label == "SUN":
#                    print particle0.label
#                    print potential_energy   # Check to see if it matches rough calculations and it does!
                total_potential = potential_energy+total_potential


                #updating force arrays using force consideration F(x1 -> x2) = -F(x2 -> x1)
                new_force_array[x1] = new_force_array[x1]+interaction_force
                new_force_array[x2] = new_force_array[x2]-interaction_force
        ###Done calculating total force from all secondary planets on main planet###

        particle0.newvel(dt,0.5*(force_array[x1]+new_force_array[x1]))        #calculates the new velocity, using verlet method, of main planet
        force_array[x1] = new_force_array[x1]
        particle0.newnewpos(dt,force_array[x1])
        
        del Planet_list[x1]
        Planet_list.insert(x1,particle0)
        new_force_array[x1] = [0,0,0]                                         #set to 0 because at each new timestep entire forces need to be recalculated

    total_energy = total_kinetic+total_potential
    Energy_list.append(total_energy)
    time_list.append((i*dt)/(60.0*60.0*24.0*365.25))
    kinetic_energy = 0.0
    total_kinetic = 0.0
    total_potential =0.0
    total_energy = 0.0
    time = time + dt

for x1 in xrange(Num_Bodies):
    particle0 = Planet_list[x1]
    if particle0.label == "MOON":
        Orbit_info.write('{} \n'.format(particle0.label))
        Orbit_info.write('The Apoapse of this body is {}km\n'.format(float(ApoapseMoon[0])))
        Orbit_info.write('The Periapse of this body is {}km\n'.format(float(PeriapseMoon[0])))
        Orbit_info.write('The Orbital Period of this body is {} days \n\n'.format( float(sum(PeriodMoon)) / float(len(PeriodMoon)) ) )
    if x1 != moonnumber and x1 != sunnumber:
        Orbit_info.write('{} \n'.format(particle0.label))
        Orbit_info.write('The Apoapse of this body is {}km\n'.format(float(Apoapsis[x1])))
        Orbit_info.write('The Periapse of this body is {}km\n'.format(float(Periapse[x1])))
        Orbit_info.write('The Orbital Period of this body is {} days '.format( float(sum(Period[x1])) / float(len(Period[x1])) ) )
        if float(sum(Period[x1])) / float(len(Period[x1])) > 365.25:
            Orbit_info.write('(or {} years)'.format( float(sum(Period[x1])) / (365.25*float(len(Period[x1]))) ) )
        Orbit_info.write('\n\n')

plt.plot(time_list,Energy_list)
plt.title('The Total Energy of the Particle as a Function of Time')
plt.xlabel('Time (Years)')
plt.ylabel('Energy (Mega Joules)')

fig = plt.gcf()
fig.set_size_inches(14, 6)
plt.savefig('Total Energy Fluctuations.png', dpi=300)

VMDFILE.close()
Orbit_info.close()

