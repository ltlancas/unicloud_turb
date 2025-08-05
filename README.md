# unicloud_turb
Generate a gaussian random velocity field for initial conditions for a turbulent molecular cloud.

To generate the 100 realizations, run `python generate.py`, which will generate a set of `.h5` HDF5 outputs storing the Cartesian cell-centered coordinates and corresponding velocity components. For examples of how to read these outputs see `generate.ipynb`. 

## How to use this to generate a standard GMC initial condition
To realize the standard GMC setup agreed upon at Numerical Recipes in Star Formation:
 1. Initialize a uniform-density sphere of radius `R` in your simulation setup in a box of side-length `4R`, filled with a static, ambient medium of density 1/1000 the density of the sphere.
 2. Re-scale and re-center the coordinates provided here so that your sphere is inscribed in the box the random velocities are defined on. That is, the velocity cube provided should be used for the "inner" cube of side length `2R` centered on the center of the domain, see schematic below.
 3. Interpolate velocities from the box to your uniform density sphere (within a distance `R` from the center of the domain), making sure to leave the velocities outside this sphere as zero.
<img width="914" height="600" alt="Screenshot 2025-08-05 at 11 08 12 AM" src="https://github.com/user-attachments/assets/c272af65-bb81-4c5d-90a4-f2ef0b421418" />
