#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 21:32:28 2017

@author: jv
"""

from Lib.Electrodes import Sphere

#example use of classes
mysphere = Sphere('Mysphere', [0,0,0],1,100) #creates an electrode, more specifically a sphere
mysphere.printname()
#test if point is inside the object that we made 
print(mysphere.ispointinobject(0,0,0))
print(mysphere.ispointinobject(0,2,0))
print(mysphere.ispointinobject(2,0,0))

#now find closest distance to object
print(mysphere.getclosestdistancetopoint(0,0,0))
print(mysphere.getclosestdistancetopoint(1,1,1))


#extend the electrode class to hold all possible shapes of electrodes
#make use of combinations of primitive electrodes to create more complicated ones
#avoid duplication of code 

#now make a list holding several objects of the Electrode class
#run over the list and check if a point is in any of the electrodes

#could even make a special class holding a list of electrodes and knowing how to tell a point is in one of its electrodes

#then create a timing loop to test the implementation for timing
#tip: first check can always be, assume the object is a sphere located on the center of the object with a radius equal to the furthest point from the center
#if the point with which to compare is further out then this, then no need to dig deeper
