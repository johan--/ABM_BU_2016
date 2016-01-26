# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# import modules
import numpy
from random import randint
import matplotlib.pyplot as plt

# define sea size and array
sea_size = 250
sea = numpy.zeros((sea_size,sea_size))

# sets Y:Y, X:X range where there is food and amount
sea[175:220, 175:220]=20
sea[125:150, 125:150]=20
sea[50:75, 50:75]=20

# set number of Sharks
sharks = 5

# set start positions 
xpos = [125]*sharks
ypos = [125]*sharks

# variable to hold count how much they eat
sharks_food = [0]*sharks

# initial heading
heading = [randint(0,8)]*sharks

# std
std = 1

def makefig():
    plt.clf()
    plt.imshow(sea)
    plt.scatter(ypos,xpos)
    axes = plt.gca()
    axes.set_xlim([0,sea_size])
    axes.set_ylim([0,sea_size])

plt.ion()
fig=plt.figure()

#plt.imshow(sea)

for tstep in range(800):
    for shark in range(sharks):
        
        # NO FOOD, expend less energy
        if sea[xpos[shark],ypos[shark]] == 0:
            sharks_food[shark] - 1 # energy loss
            change_direction = int(numpy.random.normal(0, std)) # creates a variable with value -1, 0, or 1
            heading[shark] = heading[shark] + change_direction # altering the heading based on the random number
            
            if heading[shark] == 1: # note double equals sign - this is the test for equality. Single equals means assignment (x = x+1 is assignment, x == x+1 is incorrect)
                ypos[shark] = ypos[shark] + 1 # move upwards
            if heading[shark] == 2:
                xpos[shark] = xpos[shark]+1
                ypos[shark] = ypos[shark]+2 # up and right
            if heading[shark] == 3:
                xpos[shark] = xpos[shark]+1 # right
            if heading[shark] == 4:
                xpos[shark] = xpos[shark] + 1
                ypos[shark] = ypos[shark] -1
            if heading[shark] == 5:
                ypos[shark] = ypos[shark] - 1
            if heading[shark] == 6:
                ypos[shark] = ypos[shark] - 1
                xpos[shark] = xpos[shark] -1
            if heading[shark] == 7:
                xpos[shark] = xpos[shark] - 1
            if heading[shark] == 8:
                xpos[shark] = xpos[shark]-1
                ypos[shark] = ypos[shark]+2
        
        # FOOD present
        if sea[xpos[shark],ypos[shark]] != 0:
            sharks_food[shark] + 1 # energy loss
            sea[xpos[shark],ypos[shark]] = sea[xpos[shark],ypos[shark]] - 2 # food decrease
            change_direction = int(numpy.random.normal(0, std*3)) # creates a variable with value -1, 0, or 1
            heading[shark] = heading[shark] + change_direction # altering the heading based on the random number
            
            if heading[shark] == 1: 
                ypos[shark] = ypos[shark] + 1 # move upwards
            if heading[shark] == 2:
                xpos[shark] = xpos[shark]+1
                ypos[shark] = ypos[shark]+1 # up and right
            if heading[shark] == 3:
                xpos[shark] = xpos[shark]+1 # right
            if heading[shark] == 4:
                xpos[shark] = xpos[shark] + 1
                ypos[shark] = ypos[shark] -1
            if heading[shark] == 5:
                ypos[shark] = ypos[shark] - 1
            if heading[shark] == 6:
                ypos[shark] = ypos[shark] - 1
                xpos[shark] = xpos[shark] -1
            if heading[shark] == 7:
                xpos[shark] = xpos[shark] - 1
            if heading[shark] == 8:
                xpos[shark] = xpos[shark]-1
                ypos[shark] = ypos[shark]+1
            
        ### this next section of code is 'error checking' I need 'heading' to stay between 1 and 8 (one of the 8 possible directions including diagonal)
        # this value between 1 and 8 is vital in the next bit of code - so I make sure it will always be in this range here
        if heading[shark] > 8:
            heading[shark] = heading[shark] - 8
        if heading[shark] < 1:
            heading[shark] = heading[shark] + 8

            
        ##### error checking - make sure it does not exceed the limits of the environment
        if xpos[shark] >= (sea_size-1):
            xpos[shark] = 1
        if xpos[shark] < 1:
            xpos[shark] = (sea_size-1)
        if ypos[shark] >= (sea_size-1):
            ypos[shark] = 1
        if ypos[shark] < 1:
            ypos[shark] = (sea_size-1)

    makefig()
    
    plt.pause(0.1)
    
        #plt.xlim(0,500)
        #plt.ylim(0,500)
            
            
        