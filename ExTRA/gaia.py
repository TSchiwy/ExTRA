import numpy as np
import matplotlib.pyplot as plt
from .astrometry import *
from .hipparcos import *
from .useful import *

#### Psi=pi/2 - theta #psi=HIP2 theta=GAIADR4

def angle_trafo(theta):
    psi=-theta+np.pi/2

    psi=((psi + 180) % 360) - 180
    return psi

def gaia_JD(gaia_ad,Sepoch=None):
    if Sepoch==None:
        Sepoch=J2017()
    A3,A4,A5,A6,A7,A8,A9=gaia_ad
    frac=A7/A4
    epoch=frac+2017.5
    JD=J2017()+(epoch-2017.5)*365.25 #JD for standard epoch J2017.5
    return JD


def res_to_orbit_gaia(residuals,gaia_ad,orbitfit):
    
     
    #Hipparchos
    

    t=gaia_JD(gaia_ad)
    res=gaia_ad.copy()
    res[-2]=residuals

    

    res_2D=hip_2d(res)

    hip_x=res_2D[0]
    hip_x_err=res_2D[1]
    hip_y=res_2D[2]
    hip_y_err=res_2D[3]

    #print(res_2D)
    orb_x,orb_y=orbit(*orbitfit,t)

    res_orb_x=hip_x+orb_x
    res_orb_y=hip_y+orb_y

    return res_orb_x,hip_x_err,res_orb_y,hip_y_err