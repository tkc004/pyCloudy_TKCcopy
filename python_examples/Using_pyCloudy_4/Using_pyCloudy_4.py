#!/usr/bin/env python
# coding: utf-8

# # How to take account of the slit position when computing line intensities (even for a spherical nebula)

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _example_utils import find_cloudy_exe, save_fig

import numpy as np

script_dir = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(script_dir / "temp_models" / ".mplconfig"))
os.environ.setdefault("XDG_CACHE_HOME", str(script_dir / "temp_models" / ".cache"))
import matplotlib.pyplot as plt




import pyCloudy as pc
# Changing the location and version of the cloudy executable.
cloudy_exe = find_cloudy_exe(script_dir)
pc.config.cloudy_exe = str(cloudy_exe)
from pyCloudy.utils.astro import conv_arc




# The directory in which we will have the model
# You may want to change this to a different place so that the current directory
# will not receive all the Cloudy files.
temp_model_dir = script_dir / 'temp_models'
temp_model_dir.mkdir(exist_ok=True)
fig_dir = script_dir / 'figures'
fig_dir.mkdir(exist_ok=True)
dir_ = str(temp_model_dir) + '/'




# Define some parameters of the model:
model_name = 'model_4'
full_model_name = '{0}{1}'.format(dir_, model_name)
dens = 4. #log cm-3
Teff = 45000. #K
qH = 47. #s-1
r_min = 5e16 #cm
dist = 1.26 #kpc




# these are the commands common to all the models (here only one ...)
options = ('no molecules',
            'COSMIC RAY BACKGROUND',
            )




emis_tab = ['H  1  4861.32A',
            'H  1  6562.80A',
            'Ca B  5875.64A',
            'N  2  6583.45A',
            'O  1  6300.30A',
            'O  2  3726.03A',
            'O  2  3728.81A',
            'O  3  5006.84A',
            'O  3  4363.21A',
            'O 3R  4363.00A',
            'O 3C  4363.00A',
            'S  2  6716.44A',
            'S  2  6730.82A',
            'Cl 3  5517.71A',
            'Cl 3  5537.87A',
            'O  1  63.1679m',
            'O  1  145.495m',
            'C  2  157.636m']




abund = {'He' : -0.92, 'C' : 6.85 - 12, 'N' : -4.0, 'O' : -3.40, 'Ne' : -4.00, 
         'S' : -5.35, 'Ar' : -5.80, 'Fe' : -7.4, 'Cl' : -7.00}




# Defining the object that will manage the input file for Cloudy
c_input = pc.CloudyInput(full_model_name)




# Filling the object with the parameters
# Defining the ionizing SED: Effective temperature and luminosity.
# The lumi_unit is one of the Cloudy options, like "luminosity solar", "q(H)", "ionization parameter", etc... 
c_input.set_BB(Teff = Teff, lumi_unit = 'q(H)', lumi_value = qH)




# Defining the density. You may also use set_dlaw(parameters) if you have a density law defined in dense_fabden.cpp.
c_input.set_cste_density(dens)




# Defining the inner radius. A second parameter would be the outer radius (matter-bounded nebula).
c_input.set_radius(r_in=np.log10(r_min))
c_input.set_abund(ab_dict = abund, nograins = True)
c_input.set_other(options)
c_input.set_iterate() # (0) for no iteration, () for one iteration, (N) for N iterations.
c_input.set_sphere() # () or (True) : closed geometry, or (False): open geometry.
c_input.set_emis_tab(emis_tab) # better use read_emis_file(file) for long list of lines, where file is an external file.
c_input.set_distance(dist=dist, unit='kpc', linear=True) # unit can be 'kpc', 'Mpc', 'parsecs', 'cm'. If linear=False, the distance is in log.




# Writing the Cloudy inputs. to_file for writing to a file (named by full_model_name). verbose to print on the screen.
c_input.print_input(to_file = True, verbose = False)




# Running Cloudy with a timer. Here we reset it to 0.
pc.log_.timer('Starting Cloudy', quiet = True, calling = 'test1')
c_input.run_cloudy()
pc.log_.timer('Cloudy ended after seconds:', calling = 'test1')




c_output = pc.CloudyModel(full_model_name)
c_output.print_stats()




# define the size of the 3D cube and instanciate the object that manage it.
cube_size = 201
M_sphere = pc.C3D(c_output, dims=cube_size, center=True, n_dim=1)




# plot the image of the OIII emission
plt.imshow(M_sphere.get_emis('O__3_500684A').sum(0));
save_fig(plt.gcf(), fig_dir / 'o3_image.png')




# A function in form of lambda to transform size in cm into arcsec, for a distance "dist" defined above.
arcsec = lambda cm: conv_arc(dist=dist, dist_proj=cm)




def make_mask(ap_center=[0., 0.], ap_size=[1., 1.]):
    """
    This returns a mask (values between 0. and 1.) to be multiplied to the image to take the flux passing through an aperture.
    An pc.C3D object named M_sphere must exist outside theis function
    """
    x_arc = arcsec(M_sphere.cub_coord.x_vec)
    y_arc = arcsec(M_sphere.cub_coord.y_vec)
    z_arc = arcsec(M_sphere.cub_coord.z_vec)
    X, Y = np.meshgrid(y_arc, x_arc)
    bool_mask = ((X > ap_center[0] - ap_size[0]/2.) & 
            (X <= ap_center[0] + ap_size[0]/2.) & 
            (Y > ap_center[1] - ap_size[1]/2.) & 
            (Y <= ap_center[1] + ap_size[1]/2.))
    mask = np.zeros_like(X)
    mask[bool_mask] = 1.0
    return mask




# we define the mask. Can be change to see the effect of the aperture on line intensities
mask = make_mask(ap_center=[1.5, 2.3], ap_size=[50, 1.5])




# Check that the mask is not empty
print(mask.size)
print(mask.sum())




# We plot the OIII image and overplot the mask.
plt.imshow(M_sphere.get_emis('O__3_500684A').sum(0), interpolation='None')
plt.colorbar()
plt.contour(mask);
save_fig(plt.gcf(), fig_dir / 'o3_aperture.png')




# Hbeta is computed for the whole object and throught the aperture
Hb_tot = (M_sphere.get_emis('H__1_486132A')*M_sphere.cub_coord.cell_size).sum()
Hb_slit = ((M_sphere.get_emis('H__1_486132A')*M_sphere.cub_coord.cell_size).sum(1) * mask).sum()
print(Hb_tot, Hb_slit)
if Hb_tot == 0 or Hb_slit == 0:
    raise ValueError('Hbeta flux is zero for the total model or selected aperture')




# For every line, we compute the intensity for the whole object and throught the aperture.
# We also print out the difference due to the slit.
for label in M_sphere.m[0].emis_labels:
    I_tot = np.divide(M_sphere.get_emis(label).sum()*M_sphere.cub_coord.cell_size, Hb_tot)
    I_slit = np.divide((M_sphere.get_emis(label).sum(1) * mask).sum()*M_sphere.cub_coord.cell_size, Hb_slit)
    print('line: {0:12s} I/Ib Total: {1:6.4f} I/Ib Slit: {2:6.4f} Delta: {3:4.1f}%'.format(label, I_tot, I_slit, 
                                                                                           (I_slit-I_tot)/I_tot*100))
