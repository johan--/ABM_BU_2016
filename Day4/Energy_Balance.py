
# coding: utf-8

# This is my energy balance script

# In[1]:

# number python
import numpy 

# imports plotting library
import matplotlib.pyplot as plt 

# lets you plot in a cell
#get_ipython().magic('matplotlib inline')


# In[2]:

# initial params
temp = [290, 295, 300, 305, 310]
refT = 294
devT = 5 # days
boltz = (8.62 * (10 ** -5)) ## ** power
activation_energy = -0.62


# In[3]:

# agents
egg = [] # development of egg


# In[11]:

# define a python function for arreneus 
def arreneus(temperature): # brackets hold values to pass to the function
    x = 1/(numpy.exp(((1/refT)-(1/temperature)) * ((activation_energy / boltz)))*devT)
    return x # local variable only


# In[12]:

for t in temp: # for loop ie for each temperature in vairable temp
    print(t)
    egg.append(arreneus(t)) # append result of arreneus(t in temp) to egg 


# In[13]:

days = [1]*len(egg) # set count of days for development as 1 for number of eggs

for (c,e) in enumerate(egg): # "enerumate" this counts how many loops
    while e <= 1: # while egg development is less than equal to 1
        e = e+e # egg n equals itself + itself
        days[c]=days[c]+1 # +1 to no of days


# In[14]:

plt.scatter(temp, days, marker='x')

plt.xlabel('daily T in kelvin')
plt.ylabel('days to development')
plt.title('egg!!!')

plt.show()
#plt.savefig('path')


# In[15]:

# end


# In[ ]:



