import numpy as np
from .useful import *

def mur_to_vrad(mu_r,parallax):
    return mu_r*astro_unit/ parallax



#accumulated shift over period delta_t compared to standard epoch
def secular_shift(sss,v_r,delta_t):
    
    theta=np.arctan(sss[3]/sss[4])

    mu_total=(sss[3]**2+sss[4]**2)**0.5

    delta_pos= -1*v_r*sss[2]*mu_total*delta_t*abs(delta_t) /astro_unit_z
    delta_asc_star=delta_pos*np.cos(theta)
    delta_dec=delta_pos*np.sin(theta)

    
    return np.array([delta_asc_star,delta_dec]) # in mas 



#shift of parallax and proper motion per year
def secular_acceleration(sss,v_r):
    
    theta=np.arctan(sss[3]/sss[4])

    mu_total=(sss[3]**2+sss[4]**2)**0.5
    # Scalar formula for secular acceleration
    dot_mu = -2*v_r*sss[2]*mu_total/astro_unit_z 

    delta_mu_asc=dot_mu *np.cos(theta)
    delta_mu_dec=dot_mu *np.sin(theta)

    delta_par=-1*sss[2]**2 *v_r/astro_unit_z


    
    return np.array([delta_par,delta_mu_asc,delta_mu_dec]) # in mas/yr , mas/yr^2 , mas/yr^2


