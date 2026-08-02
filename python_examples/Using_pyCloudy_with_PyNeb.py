#!/usr/bin/env python
# coding: utf-8

# # Changing atomic data using PyNeb

# It is possible to extract from the Cloudy model the electron temperature and density and the ionic fractions to re-compute at each zone of the nebula the emissivities of the lines, using the PyNeb code.
# This is NOT coherent in the fact that changing the line emissivities change the cooling and then the electron temeprature. And only collisional effects are taken into account. But this can nevertheless helps to understand the effect of choosing one set of atomic data or another one in the analysis of a nebula.

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import pyCloudy as pc
import pyneb as pn
import os
home_dir = os.environ['HOME'] + '/'
pc.config.cloudy_exe = '/usr/local/Cloudy/c25.00_rc2/source/cloudy.exe'


# In[2]:


# We are using the model from the example 1
Mod = pc.CloudyModel('/tmp/models/model_1')


# In[3]:


# Print some data about the model
Mod.print_stats()


# In[8]:


# Print all the different atomic data avilable in Pyneb for the [OIII] lines
print(pn.atomicData.getAllAvailableFiles('O3',data_type='atom', mark_current=False))
print('--------------------------------------------------')
print(pn.atomicData.getAllAvailableFiles('O3',data_type='coll', mark_current=False))


# In[10]:


pc.log_.level=1
pn.log_.level=2
# Loops on the different As. 
i = 0
f, ax = plt.subplots()
for O3_atom in pn.atomicData.getAllAvailableFiles('O3',data_type='atom', mark_current=False):
    pn.atomicData.setDataFile(O3_atom) # Change the datafile used in PyNeb
    O3 = pn.Atom('O',3, NLevels=6)
    Mod.add_emis_from_pyneb('new_a5007_{}'.format(i), O3, wave=5007)
    Mod.add_emis_from_pyneb('new_a4363_{}'.format(i), O3, wave=4363)
    ax.plot(Mod.radius, Mod.get_emis('new_a5007_{}'.format(i))/Mod.get_emis('new_a4363_{}'.format(i)), label=O3_atom) # Plot the diagnostic ratio
    i += 1
ax.set_xlabel('Radius [cm]')
ax.set_ylabel('[OIII] 5007/4363')
ax.legend(loc=3)
ax.set_ylim((0., 600));


# In[11]:


pc.log_.level=1
pn.log_.level=2
i = 0
f, ax = plt.subplots()
# The same but changing the collision strengths
for O3_coll in pn.atomicData.getAllAvailableFiles('O3',data_type='coll', mark_current=False):
    pn.atomicData.setDataFile(O3_coll)
    O3 = pn.Atom('O',3, NLevels=6)
    Mod.add_emis_from_pyneb('new_c5007_{}'.format(i), O3, wave=5007)
    Mod.add_emis_from_pyneb('new_c4363_{}'.format(i), O3, wave=4363)
    ax.plot(Mod.radius, Mod.get_emis('new_c5007_{}'.format(i))/Mod.get_emis('new_c4363_{}'.format(i)), label=O3_coll) # Plot the diagnostic ratio
    i += 1
ax.set_xlabel('Radius [cm]')
ax.set_ylabel('[OIII] 5007/4363')
ax.legend(loc=3)
ax.set_ylim((0., 600));


# In[6]:


Mod.emis_labels


# In[12]:


pc.log_.level=1
pn.log_.level=2
# Define the data that will be used to compute Te
pn.atomicData.setDataFile('o_iii_coll_SSB14.dat')
pn.atomicData.setDataFile('o_iii_atom_FFT04.dat')
O3 = pn.Atom('O',3, NLevels=6)
i = 0
for O3_coll in pn.atomicData.getAllAvailableFiles('O3',data_type='coll', mark_current=False):
    tem_diag = Mod.get_emis_vol('new_c5007_{}'.format(i))/Mod.get_emis_vol('new_c4363_{}'.format(i))
    tem = O3.getTemDen(tem_diag, den = 1e4, wave1 = 5007, wave2 = 4363)
    print('{0:27s} [OIII]5007/4363 = {1:5.1f} Te = {2:6.1f}'.format(O3_coll, tem_diag, tem))
    i += 1
pn.atomicData.setDataFile('o_iii_coll_AK99.dat')
i = 0
for O3_atom in pn.atomicData.getAllAvailableFiles('O3',data_type='atom', mark_current=False):
    tem_diag = Mod.get_emis_vol('new_a5007_{}'.format(i))/Mod.get_emis_vol('new_a4363_{}'.format(i))
    tem = O3.getTemDen(tem_diag, den = 1e4, wave1 = 5007, wave2 = 4363)
    print('{0:27s} [OIII]5007/4363 = {1:5.1f} Te = {2:6.1f}'.format(O3_atom, tem_diag, tem))
    i += 1


# In[13]:


print(pn.atomicData.getAllAvailableFiles('S2',data_type='atom', mark_current=False))
print('--------------------------------------------------')
print(pn.atomicData.getAllAvailableFiles('S2',data_type='coll', mark_current=False))


# In[14]:


i = 0
f, ax = plt.subplots()
for S2_atom in pn.atomicData.getAllAvailableFiles('S2',data_type='atom', mark_current=False):
    pn.atomicData.setDataFile(S2_atom)
    S2 = pn.Atom('S',2, NLevels=6)
    Mod.add_emis_from_pyneb('new_a6716_{}'.format(i), S2, wave=6716)
    Mod.add_emis_from_pyneb('new_a6731_{}'.format(i), S2, wave=6731)
    ax.plot(Mod.radius, Mod.get_emis('new_a6716_{}'.format(i))/Mod.get_emis('new_a6731_{}'.format(i)), label=S2_atom)
    i += 1
ax.set_xlabel('Radius [cm]')
ax.set_ylabel('[SII] 6716/6731')
ax.legend(loc=2);


# In[15]:


i = 0
for S2_atom in pn.atomicData.getAllAvailableFiles('S2',data_type='atom', mark_current=False):
    dens_diag = Mod.get_emis_vol('new_a6716_{}'.format(i))/Mod.get_emis_vol('new_a6731_{}'.format(i))
    dens = S2.getTemDen(dens_diag, tem=1e4, wave1=6716, wave2=6731)
    print('{0:27s} [SII]6716/31 {1:5.3f}, density = {2:5.1f}'.format(S2_atom, dens_diag, dens))
    i += 1

