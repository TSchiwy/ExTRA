# ExTRA
**Ex**oplanetary **T**ool for **R**adial Velocities and **A**strometry is a simple collection of functions to work with **Hipparcos** and **Radial Velocity** data.

It provides log-likelihood functions for both solo and combined RV/HIP fits. It can also be used to visualize Hipparcos data in 2 dimensions, as well as for creating mock RV and HIP data.

---

## Data Structure

### HIP

First of all, it is important to note that this code was made for the **SECOND** Hipparcos Reduction. The data from the first reduction won't work with this code.

The Hipparcos IAD is generally given in columns as:

- **A1** : Orbit number of Hipparcos (not important)
- **A2** : Epoch (Year-1991.25)
- **A3** : cos(θ)
- **A4** : sin(θ)
- **A8** : Abscissa residual [mas]
- **A9** : Formal error on abscissa residual [mas]

By multiplying the Epoch (A2) with A3 and A4 we get:

- **A6** : cos(θ) * t
- **A7** : sin(θ) * t

These will be useful when recomputing proper motions later.

### Reading Hipparcos Data

If you want to read Hipparcos data from the second reduction, you can do so using the function:

```python
ExTRA.hip_read("path")

which return the data in the format of:

data[0]=A_3,A_4,A_5,A_6,A_7,A_8,A_9  
data[1]=t_HIP





