import numpy as np
import matplotlib.pyplot as plt
from .astrometry import *
from .hipparcos import *

#### Psi=pi/2 - theta #psi=HIP2 theta=GAIADR4

def angle_trafo(theta):
    psi=-theta+np.pi/2

    psi=((psi + 180) % 360) - 180
    return psi