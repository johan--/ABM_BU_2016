import numpy 

import matplotlib.pyplot as plt # functions for graphs - making dispalys easier - the 'as plt' means we only need to refer to it as plt in the future

from random import randint

def makeFig(): # a user built function
    plt.clf()
    plt.scatter(xpos,ypos)
    axes = plt.gca()
    axes.set_xlim([0,env_size])
    axes.set_ylim([0,env_size])
    
plt.ion() # enable interactivity
fig=plt.figure() # make a figure

########################################
###Parameters
number_agents = 30
env_size = 99
model_time = 100
#########################################

xpos = numpy.zeros(number_agents) 
ypos = numpy.zeros(number_agents) 
heading = numpy.zeros(number_agents) 

for x in range(number_agents):
    xpos[x] = randint(0,env_size)
    ypos[x] = randint(0,env_size)
    heading[x] = randint(1,8)


for timestep in range(model_time): # the model will run for however long is specificed by the model_time variable - initially four timesteps
    for agent in range(number_agents): # there are four agents - zero to 3
        # note that this line is now indented twice - anything here runs in the agent loop - but the agent loop itself runs every timestep
        change_direction = randint(-1,1) # creates a variable with value -1, 0, or 1
        heading[agent] = heading[agent] + change_direction # altering the heading based on the random number
    
        ### this next section of code is 'error checking' I need 'heading' to stay between 1 and 8 (one of the 8 possible directions including diagonal)
        # this value between 1 and 8 is vital in the next bit of code - so I make sure it will always be in this range here
        if heading[agent] > 8:
            heading[agent] = heading[agent] - 8
        if heading[agent] < 1:
            heading[agent] = heading[agent] + 8
        ## end of the error checking - I now have a direction - up, down, left, right and the diagonals. 
        
        if heading[agent] == 1: # note double equals sign - this is the test for equality. Single equals means assignment (x = x+1 is assignment, x == x+1 is incorrect)
            ypos[agent] = ypos[agent] + 1 # move upwards
        if heading[agent] == 2:
            xpos[agent] = xpos[agent]+1
            ypos[agent] = ypos[agent]+1 # up and right
        if heading[agent] == 3:
            xpos[agent] = xpos[agent]+1 # right
        if heading[agent] == 4:
            xpos[agent] = xpos[agent] + 1
            ypos[agent] = ypos[agent] -1
        if heading[agent] == 5:
            ypos[agent] = ypos[agent] - 1
        if heading[agent] == 6:
            ypos[agent] = ypos[agent] - 1
            xpos[agent] = xpos[agent] -1
        if heading[agent] == 7:
            xpos[agent] = xpos[agent] - 1
        if heading[agent] == 8:
            xpos[agent] = xpos[agent]-1
            ypos[agent] = ypos[agent]+1       
           
        # now error check for positions - keep between 0 and 99
           
        
        if xpos[agent] >= (env_size):
            xpos[agent] = 1
        if xpos[agent] < 1:
            xpos[agent] = (env_size)
        if ypos[agent] >= (env_size):
            ypos[agent] = 1
        if ypos[agent] < 1:
            ypos[agent] = (env_size)
        
    ### where indentation gets difficult - this is now out of the agent, but in the timestep loop
    makeFig()      

    plt.pause(0.1)
    
