import numpy as np
import scipy.constants
from .useful import *
from .vectorastrometry import spherical_to_cartesian, normal_triad

def mu_to_v(parallax,mu):
    return mu*astro_unit/ parallax



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





def ltd_approx(v,format="years"): #lighttimedifference

    
    if format=="years":
        delta_t=(v*365.25*24*60*60)/lightyear  # in 1/YEARS 
        return delta_t
    if format=="days":

        delta_t_day=(v*24*60*60)/lightyear # in 1/DAYS 
        return delta_t_day



def ltd_accurate(t,sss,v_rad,format="years"):

    pc_to_km = 3.085677581e13
    d0_pc=1000/sss[2] #distance in pc using parallax
    d0=d0_pc*pc_to_km #distance in km

    

    r0_vec=spherical_to_cartesian(d0,sss[0],sss[1]) #cartesian vector to source
    r0=np.linalg.norm(r0_vec)

    

    v_asc,v_dec=mu_to_v(sss[2],np.array([sss[3],sss[4]])) #mu to v using parallax
    triad=np.array(normal_triad(sss[0],sss[1]))

    if format=="years":
        v0_vec=np.array([v_asc,v_dec,v_rad])*365.25*24*60**2 #spherical velocities in km/year
    if format=="days":
        v0_vec=np.array([v_asc,v_dec,v_rad])*24*60**2
    
    v0_cartesian=triad @ v0_vec #cartesian velocity



   


    tau=1000*(np.linalg.norm(r0_vec+v0_cartesian*t)-r0)/(scipy.constants.c) #c in kms, v in kms r0 in pc


    return tau


    




