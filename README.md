# ExTRA
The **Ex**oplanetary **T**ool for **R**adial Velocities and **A**strometry is a simple collection of functions to work with Hipparcos and Radial Velocity data.
It provides Loglikelihood functions for both solo and combined RV/HIP fits. It can also be used to visualise Hipparcos data 2 dimensionally as well as for creating mock RV and HIP Data.



## Data Structure
### HIP
First of all it is important to note that this code was made for the **SECOND** Hipparcos Reduction. The data of the first one wont work with this code.

The Hipparcos IAD is generally given in coulmns as:
A1 : Orbit number of Hipparcos (not important)
A2 : Epoch (Year-1991.25) 
A3 : cos()
A4 : sin()
A8 : Abscissa residual [mas]
A9 : Formal error on abscissa residual [mas]

By multiplying the Epoch (A2) onto A3 and A4 we get :

A6 : cos()*t
A7 : sin()*t

which will be useful when recomputing proper motions later.

One can also compute the timestamps of the measurements in JD:
t_HIP: 2451545.0+(A7/A4+1991.25-2000.0)*365.25 [JD]

If you want to read out Hipparcos data of the second reduction you can do so using:

ExTRA.hip_read("path")

which return the data in the format of:

data=A_3,A_4,A_5,A_6,A_7,A_8,A_9,t_HIP





