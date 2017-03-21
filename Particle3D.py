"""
This code creates a class called Particle3D. This will let us model the 3D behaviour of a particle using two different methods, Euler and Verlet. It also lets us model the interactions between two particles.

Begin by importing the necessary modules, including the Vector module from the previous exercise
"""
import math as m
import numpy as np
import Vector as vctr
import sys

"""
Create a class called Particle3D
"""
class Particle3D(object):

    def __init__(self,mass,label,x_pos,y_pos,z_pos,x_vel,y_vel,z_vel):
	"""
	:param mass: the mass of the particle
	:param label: the particle number
	:param x_pos: the x position of the particle
	:param y_pos: the y position of the particle
	:param z_pos: the z position of the particle
	:param x_vel: the x component of the velocity of the particle
	:param y_vel: the y component of the velocity of the particle
	:param z_vel: the z component of the velocity of the particle
	:return: Lets the data be called and arranges the position and velocity as arrays
	"""
        self.mass = mass
        self.label = label
        self.position = np.array([x_pos,y_pos,z_pos])
        self.velocity = np.array([x_vel,y_vel,z_vel])

#    def __str__(self):
#        """
#        :return: Returns the particle number and a position vector, useful for printing
#        """
#        return "label = " + str(self.label) + ", r = " + str(self.position)

    #how str is defined is changed so that it is a position in a VMD file
    def __str__(self):
        self.positionvmd = str(self.position[0]) +" " +str(self.position[1]) +" " +str(self.position[2])
        return "{} {} \n".format(self.label,self.positionvmd)

    def compute_KE(self):
        """
        :return: Returns the kinetic energy of the particle
        """
        return (1.0/2.0) * self.mass * vctr.SqMag(self.velocity)

    def new_vel(self,dt,f):
        """
        :param dt: the time step
        :param f: the force acting on the particle
        return: Returns the new velocity of the particle
        """
        self.velocity = self.velocity + dt*f/self.mass
       
    def orderone_pos(self,dt):
        """
        :param dt: the time step
        :return: Returns the new position using the Euler method
        """
        self.position = self.position + dt*self.velocity

    def ordertwo_pos(self,dt,f):
        """
        :param dt: the time step
        :param f: the force acting on the particle
        return: Returns the new position using the Verlet method
        """
        self.position = self.position +dt*self.velocity + dt**2.0*f/(2.0*self.mass)

    @staticmethod
    def from_file(file_name):
        """
        :param file_name: the name of the file which contains the data
        return: Returns the data from the file as a particle, using the Particle3D class
        """
        tokens = file_name
        return Particle3D(float(tokens[0]), str(tokens[1]), float(tokens[2]), float(tokens[3]), float(tokens[4]), float(tokens[5]), float(tokens[6]), float(tokens[7]))

    @staticmethod
    def vec_sep(self,particle2):
        """
        :param particle2: the second particle included in the system
        :return: Returns the vector seperation of the two particles
        """
        vecsep = self.position - particle2.position
        return vecsep

    @staticmethod
    def mag_sep(self,particle2):
        """
        :param particle2: the second particle included in the system
        :return: Returns the vector magnitude seperation of the two particles
        """
        magsep = vctr.Mag(self.position - particle2.position)
        return magsep

    @staticmethod
    def force(self,particle2,f):
        """
        :param particle2: the second particle included in the system
        :return: Returns the vector force from particle2 acting on self
	"""
	f = -6.67408*10**(-20)*self.mass*particle2.mass*Particle3D.vecsep(self,particle2)/(vctr.SqMag(Particle3D.vecsep(self,particle2))**(3.0/2.0))
	return f

    @staticmethod
    def poten(self,particle2,u):
        """
        :param particle2: the second particle included in the system
        :return: Returns the potential energy between the two particles
	"""
	u = -6.67408*10**(-20)*self.mass*particle2.mass/Particle3D.mag_sep(self,particle2)
	return u

    @staticmethod
    def ang_vel(self,particle2,w):
        """
        :param particle2: the second particle included in the system
        :return: returns the angular velocity with respect to particle2 as the central particle
	"""
        w = vctr.Mag(vctr.cross(Particle3D.vecsep(particle2,self),(self.velocity-particle2.velocity))/(vctr.SqMag(Particle3D.vecsep(self,particle2))) )
        return w

