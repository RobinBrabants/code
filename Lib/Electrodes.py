# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 14:56:51 2017

@author: joverbee
"""
import math
from abc import ABCMeta, abstractmethod
class Electrode(object):
    #a class to represent all electrodes from which to derive the actual shaped objects later
    #should be virtual class, but not same as in c++, need to figure out how to prevent instantiating this
    def __init__(self,name,position,potential):
        self.name=name
        self.position=position
        self.potential=potential
    def printname(self): #output name and some info on the object
        print('name: ',self.name,'\n');
    @abstractmethod
    def getclosestdistancetopoint(self,x,y,z):    #get closest distance from a point to this object (empty here as the generic object is not defined, needs to be specified in each individual object)
        pass
    @abstractmethod
    def ispointinobject(self,x,y,z):
        pass
    @abstractmethod
    def electrode_type():
        pass
    
#now instantiate specific objects and give them their specific functionality by overwriting the generic function        
           
class Sphere(Electrode):
    def __init__(self,name,position,radius,V):
        Electrode.__init__(self, name,position,V)
        self.radius=radius
    def getclosestdistancetopoint(self,x,y,z): #make specific code for this object
        return 0
    def ispointinobject(self,x,y,z):       
        return math.sqrt((x-self.position[0])**2+(y-self.position[1])**2+(z-self.position[2])**2)<self.radius
    def getclosestdistancetopoint(self,x,y,z):    #get closest distance from a point to 
       return abs(math.sqrt((x-self.position[0])**2+(y-self.position[1])**2+(z-self.position[2])**2)-self.radius)
    def electrode_type(self):
        return 'sphere electrode'