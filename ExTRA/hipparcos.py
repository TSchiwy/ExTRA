import numpy as np
import matplotlib.pyplot as plt
from .astrometry import *

#This function calculates the new abscissa residual if you change a parameter compared to the hipparchos solution.
#For example if parallax_hipparchos=40mas , parallax_model=10mas, you will get a parameter_residual=-30mas.
#This is done for all the standard model parameters.
#The value of hipp_derivation, is given in the IAD of hipparchos.Namely IA3-IA7
#They are different in the first and the second reduction
#we will be working with the second reduction ONLY

#the first old_residual is the hipparchos_residual given in mas
def abs_res(old_res,parameter_model,parameter_hipp,hipp_derivation):
    parameter_residual=parameter_model-parameter_hipp
    
    new_res=old_res
    for i in range(5):
        new_res=new_res-hipp_derivation[i]*parameter_residual[i]
    
    return new_res


def JD_hip(A4,A7): #finds JD of measurement
    frac=A7/A4
    epoch=frac+1991.25
    JD=2451545.0+(epoch-2000.0)*365.25 #JD for standard epoch J2000
    return JD


def scanangle_hip(A3,A4): #finds the scanangle between NORTH and RGC(reference great circle) in radians
    angle=np.arctan2(A4,A3)
    angle=angle
    return angle

#the scanangle opens in the direction of RA




##############################################




def rotation_counterclockwise(x,y,theta):
    x_new=+np.cos(theta)*x-np.sin(theta)*y
    y_new=+np.sin(theta)*x+np.cos(theta)*y
    #y_new=-1*y_new
    return x_new,y_new
    
    


def hip_2d(A3,A4,A8,A9): #rotates hip measurement into Dec RA* system
    """rotates a hip measurement into Ra* ,Dec, returns x,x_err,y,y_err
    
    
    
    
    
    
    
    
    
    
    
    
    
    """
    scanangle=scanangle_hip(A3,A4)
    x=rotation_counterclockwise(A8,0,scanangle)[0]
    y=rotation_counterclockwise(A8,0,scanangle)[1]
    x_err=rotation_counterclockwise(A9,0,scanangle)[0]
    y_err=rotation_counterclockwise(A9,0,scanangle)[1]
    #print(np.degrees(scanangle))
    return np.array([x,x_err,y,y_err])


#This plots the dimension hipparchos DIDNT measure in, its just nice for plots but thats it
#inputs: x0,y0 positon of measurement,x_err,y_err error in direction,get this from hip measurement,
#size= length of line
def plot_hip(x0,x_err,y0,y_err,size=1,co="r",s=10,linecolor="lightgrey"):
#this part just adds the rotation of the error so they can be plot in the right direction

    norm=1/((x_err**2 +y_err**2)**0.5)
    x1=x0-y_err*norm*size
    y1=y0+x_err*norm*size
    x2=x0+y_err*norm*size
    y2=y0-x_err*norm*size
#then plots a line for each measurement in size and angle of the error    
    plt.scatter(x0,y0,color=co,s=s)
    plt.plot([x1,x2],[y1,y2],color=linecolor)
    return
    





#This plots the dimension hipparchos measured
#inputs: x0,y0 positon of measurement,x_err,y_err error in direction,get this from hip measurement,
def plot_hip_err(x0,x_err,y0,y_err,co="r",s=1,linecolor="lightgrey"):
#this part just adds the rotation of the error so they can be plot in the right direction

    x1=x0-x_err
    y1=y0-y_err
    x2=x0+x_err
    y2=y0+y_err
#then plots a line for each measurement in size and angle of the error    
    plt.scatter(x0,y0,color=co,s=s)
    plt.plot([x1,x2],[y1,y2],color=linecolor)
    return






#hipp measurements with error, t_hip is an array of the hipparchos measurement times
#ALL INPUTS MUST BE FROM THE HIPPARCHOS CATALOGUE!
def hip_measurement(asc,dec,parallax,mu_a_star,mu_d,t_HIP,earth,A3,A4,A8,A9,tangential=1):
    """Computes RAW 2D positions of hipparcos data given the hipparcos standard solution"""
    
    #asc_star=asc*np.cos(dec)#hipp catalogue gives us asc, we need asc_star for the model
    
    


    
    #this sections rotates abscissa measurements into the Dec,RA* frame.
    d=hip_2d(A3,A4,A8,A9)
    d=np.array(d)
    x=d[0]
    x_err=d[1]
    y=d[2]
    y_err=d[3]
    
    #to get the actual hipp measurement, we need to subtract the abscissa residual
    #from the hipparchos standard model solution
    #J1991.25, or 2448349.0625JD
    if tangential==0:
        
        x_mod,y_mod=standard_model(asc,dec,parallax,mu_a_star,mu_d,t_HIP,earth,Sepoch=2448349.0625,tangential=0)
    
        x0=x_mod-x
        y0=y_mod-y
    if tangential==1:
        
        x0=x
        y0=y
    return x0,x_err,y0,y_err    


def hip_residuals(HIP_data,hip_stand,stand_fit,orbit_fit):
         
        #Hipparchos
        #First, residual due to corrections
        A3,A4,A5,A6,A7,A8,A9=HIP_data
        
        hip_ad=np.array([A3,A4,A5,A6,A7])
        
        
        
        t_HIP=JD_hip(A4,A7)
        hip_stand=np.array(hip_stand)
        
        #Since gaia has the values for J2016, we need to recompute asc_star and dec for year
        #J1991 to calculate the new residuals
        
        
        t_1991=2448349.0625
        t_2016=2457389.0
        
        asc_91,dec_91=pos_recalc(stand_fit[0],stand_fit[1],stand_fit[2],stand_fit[3],t_2016,t_1991)
    
        

                                        
        
        #The derivations in A3 compute the change for the abscissa with given asc_star! 
        #we need to multiply the catalogue values with the cos of their respective declination.
        asc_91_star=asc_91*np.cos(np.radians(dec_91))
        corr1991=(asc_91_star,dec_91,stand_fit[2],stand_fit[3],stand_fit[4])

        #print(corr1991)
        
        hip_stand[0]=hip_stand[0]*np.cos(np.radians(hip_stand[1]))
        
        hip_stand[1]=hip_stand[1]
        A8=np.array(A8)
        
        
                                         
                                         
                                         
        #the new residuals for gaia catalogue standard parameters:
        
        c_res_hip=abs_res(A8,corr1991,hip_stand,hip_ad)#abs_residual=a8

        #print(c_res_hip)
        
 
        #Now we have the remaining residuals, where the orbit will be fit too.
        #To do that, we need to correct the hip residuals again, this time for an orbit:
        x_O,y_O=orbit(orbit_fit[0],orbit_fit[1],orbit_fit[2],orbit_fit[3],orbit_fit[4],orbit_fit[5],orbit_fit[6],t_HIP)
        
        residuals_hip=c_res_hip-(A3*x_O+A4*y_O)
   
        return residuals_hip #returns res