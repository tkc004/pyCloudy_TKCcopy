#!/usr/bin/env python
# coding: utf-8



import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
home_dir = os.environ['HOME'] + '/'
import pyCloudy as pc
print(pc.__version__)




script_dir = Path(__file__).resolve().parent
cloudy_exe = None
for base_dir in (script_dir, *script_dir.parents):
    candidate = base_dir / 'Cloudy_exe' / 'Cloudy' / 'c22.02' / 'source' / 'cloudy.exe'
    if candidate.exists():
        cloudy_exe = candidate
        break
if cloudy_exe is None:
    raise FileNotFoundError('Could not find Cloudy_exe/Cloudy/c22.02/source/cloudy.exe')
pc.config.cloudy_exe = str(cloudy_exe)




temp_model_dir = script_dir / 'temp_models'
temp_model_dir.mkdir(exist_ok=True)
dir_ = str(temp_model_dir) + '/'
pc.print_make_file(dir_)




def set_models(dir_, model_name):
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
    a = 2.
    b = 1.0
    thetas = np.linspace(0., 90., 6)
    thetas_rad = np.pi / 180. * thetas
    fact_elli = a * b / np.sqrt((b * np.sin(thetas_rad))**2 + (a * np.cos(thetas_rad))**2)
    rs_in = 16.5 + np.log10(fact_elli)
    densities = 4 - np.log10(fact_elli) * 2
    
    model = pc.CloudyInput()
    model.set_BB(80000., 'q(H)', 47.3)
    model.set_grains()
    model.set_emis_tab(emis_tab)
    
    for theta, r_in, density in zip(thetas, rs_in, densities):
        model.model_name = '{0}/{1}_{2:.0f}'.format(dir_, model_name,theta)
        model.set_cste_density(density)
        model.set_radius(r_in)
        model.set_theta_phi(theta)
        model.print_input(to_file = True, verbose = False)




def def_profiles(m3d):
    """
    This uses the default velocity law (polynome) and default profile (gaussian)
    """
    m3d.set_velocity(params = [20.,60.])
    m3d.config_profile(size_spectrum = 51, vel_max = 50, v_turb = 0.01)    




def def_profiles_user(m3d):
    """
    Use this to define your own expansion velocity
    """
    def velo_polynome(params):
        """
        USer defined expansion velocity
        """
        # params is a 2 elements table, the first element is a table of parameters, the second one the cob_coord
        # which is needed to know r, x, y and z to define the velocity.
        coeffs = params[0]
        cub_coord = params[1]
        tmp = 0.
        for i, coeff in enumerate(coeffs): 
            # for each parameter we add the corresponding coeff * R**power
            tmp = tmp + coeff * cub_coord.r**i
        tmp = tmp / cub_coord.r
        # to avoid the singularity:
        tt = (cub_coord.r == 0.)
        tmp[tt] = 0
        # Projecting on each one of the 3 axes to obtain the velocity components
        vel_x = tmp * cub_coord.x / np.max(cub_coord.x)
        vel_y = tmp * cub_coord.y / np.max(cub_coord.y)
        vel_z = tmp * cub_coord.z / np.max(cub_coord.z)
        return vel_x, vel_y, vel_z
    
    def Hb_prof(x, zeta_0):
        """
        The Hbeta profile is sum of 2 blocks of lines (actually 3 + 4 lines)
        """
        res1 = .41 /zeta_0 / np.sqrt(np.pi) * np.exp(-(((x-2.7)/zeta_0)**2))
        res2 = .59 /zeta_0 / np.sqrt(np.pi) * np.exp(-(((x+2.0)/zeta_0)**2))
        return res1 + res2

    m3d.set_velocity(velocity_law='user', params = [[20.,60.], m3d.cub_coord], user_function = velo_polynome)
    m3d.config_profile(size_spectrum = 41, vel_max = 25, profile_function = Hb_prof, v_turb = 0.01)




def plot_profiles(m3d, x_pos, y_pos):
    plt.plot(m3d.vel_tab,m3d.get_profile('H__1_486132A', axis='x')[:,x_pos,y_pos] * 5, label = r'H$\beta$')
    plt.plot(m3d.vel_tab,m3d.get_profile('N__2_658345A', axis='x')[:,x_pos,y_pos] * 5, label = r'[NII]$\lambda$6584')
    plt.plot(m3d.vel_tab,m3d.get_profile('O__3_500684A', axis='x')[:,x_pos,y_pos], label = r'[OIII]$\lambda$5007')
    plt.legend()




def other_plots(m3d, proj_axis):
    plt.subplot(331)
    plt.imshow(m3d.get_emis('H__1_486132A').sum(axis = proj_axis)*m3d.cub_coord.cell_size)
    plt.title('Hb')
    plt.colorbar()
    
    plt.subplot(332)
    plt.imshow(m3d.get_emis('N__2_658345A').sum(axis = proj_axis)*m3d.cub_coord.cell_size)
    plt.title('[NII]')
    plt.colorbar()
    
    plt.subplot(333)
    plt.imshow(m3d.get_emis('O__3_500684A').sum(axis = proj_axis)*m3d.cub_coord.cell_size)
    plt.title('[OIII]')
    plt.colorbar()
    
    plt.subplot(334)
    plt.imshow(m3d.get_emis('N__2_658345A').sum(axis = proj_axis)/m3d.get_emis('H__1_486132A').sum(axis = proj_axis))
    plt.title('[NII]/Hb')
    plt.colorbar()
    
    plt.subplot(335)
    plt.imshow(m3d.get_emis('O__3_500684A').sum(axis = proj_axis)/m3d.get_emis('H__1_486132A').sum(axis = proj_axis))
    plt.title('[OIII]/Hb')
    plt.colorbar()
    
    plt.subplot(336)
    plt.imshow(m3d.get_ionic('O',1)[n_cut,:,:])
    plt.title('O+ cut')
    plt.colorbar()
    
    plt.subplot(337)
    plt.scatter(m3d.get_ionic('O',1).ravel(),m3d.get_ionic('N',1).ravel()/m3d.get_ionic('O',1).ravel(),
                c=np.abs(m3d.cub_coord.theta.ravel()), edgecolors = 'none')
    plt.title('Colored by |Theta|')
    plt.xlabel('O+ / O')
    plt.ylabel('N+/O+ / N/O')
    plt.colorbar()
    
    plt.subplot(338)
    plt.scatter(m3d.get_ionic('O',1).ravel(),m3d.get_ionic('N',1).ravel()/m3d.get_ionic('O',1).ravel(),
                c=m3d.relative_depth.ravel(),vmin = 0, vmax = 1, edgecolors = 'none')
    plt.title('Colored by position in the nebula')
    plt.xlabel('O+ / O')
    plt.ylabel('N+/O+ / N/O')
    plt.colorbar()
    
    plt.subplot(339)
    C1 = (m3d.get_ionic('N',1)/m3d.get_ionic('O',1)*m3d.get_ionic('N',2))
    C2 = (m3d.get_ionic('N',2))
    tt = (m3d.get_ionic('O',1) == 0)
    C1[tt] = 0
    C2[tt] = 0
    V = C1.sum(axis = proj_axis) / C2.sum(axis = proj_axis)
    plt.imshow(V)
    plt.colorbar()
    plt.title('N+/O+ / N/O weighted by NII')
    plt.contour(V,levels=[1.0])




model_name = "M3D_1"
pc.log_.calling = 'Model3D : ' + model_name
pc.log_.level = 3




dim = 101
n_cut = int((dim-1) /2)
proj_axis = 0




set_models(dir_, model_name)




pc.print_make_file(dir_ = dir_)
pc.run_cloudy(dir_ = dir_, n_proc = 6, model_name = model_name, use_make = True)




liste_of_models = pc.load_models('{0}/{1}'.format(dir_, model_name), list_elem=['H', 'He', 'C', 'N', 'O', 'Ar', 'Ne'],  
                                           read_cont = False, read_grains = False)




M=liste_of_models[0]
M.emis_labels




m3d = pc.C3D(liste_of_models, dims = [dim, dim, dim], angles = [45,45,0], plan_sym = True)




def_profiles(m3d)




plt.figure(figsize=(10,10))
plot_profiles(m3d, 55, 55)




plt.figure(figsize=(10,10))
plot_profiles(m3d, 55, 55)
def_profiles_user(m3d)
plt.plot(m3d.vel_tab,m3d.get_profile('H__1_486132A', axis='x')[:,55,55] * 5, ':b', label = r'H$\beta$')




plt.figure(figsize=(15,15))
other_plots(m3d, proj_axis)




im = m3d.get_RGB(list_emis = ['N__2_658345A', 'O__3_500684A', 'H__1_486132A'])
plt.figure(1, figsize=(10,10))
plt.imshow(im)




im = m3d.get_RGB(list_emis = ['N__2_658345A', 'O__3_500684A', 'H__1_486132A'])
plt.figure(1, figsize=(15,15))
plt.imshow(im)
m3d.plot_profiles(ref = 3, i_fig = 1, Nx=20, Ny=20)




f, ax = plt.subplots()
N2map = m3d.get_emis('N__2_658345A').sum(axis = proj_axis)
Hbmap = m3d.get_emis('H__1_486132A').sum(axis = proj_axis)
O3map = m3d.get_emis('O__3_500684A').sum(axis = proj_axis)
masks = []
for mapl in (Hbmap, O3map, N2map):
    masks.append(mapl > 0.01 * mapl.max())
mask = np.logical_and.reduce(masks)
ax.scatter(np.log10((N2map/Hbmap)[mask]), np.log10((O3map/Hbmap)[mask]))
ax.scatter(np.log10(N2map[mask].sum()/Hbmap[mask].sum()), np.log10(O3map[mask].sum()/Hbmap[mask].sum()), marker='*', s=200, color='red')
ax.set_xlabel('log10([NII]/Hb)')
ax.set_ylabel('log10([OIII/Hb)');
