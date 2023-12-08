# Pressure driven flow between two cylinders.
## What it is
This program allows the user to visualize pressure driven flow between two cylinders where the outer cylinder is moving.
## How to use
Download this repository and run main.jl using the Julia programing language. After running, .gif animations will be in the Animations folder.
### Modifying simulation parameters
Physical parameters for the simulation are found on lines 3-11 of main.jl Here are the default parameters and their meaning.
###
Inner cylinder radius:
`R1 = 1 `
###
Outer cylinder radius:
`R2 = 2`
###
Fluid kinematic viscosity:
`nu = .1`
###
Fluid density:
`rho = 1.0`
###
Pressure gradient:
`dPdz = -1.0`
###
Function describing the outer cylinder's angular velocity as a function of time:
`
function omega(t)
    return 0.3*sin(2*pi*t)
end
`
###
Timespan of simulation. Here from 0-5 seconds:
`tspan = (0,5)`
### Modifying visualization parameters
These parameters affect how long the program takes to run as well as how good the animations will look. Please note that running the program with the default parameters will produce high quality visualizations, but computation will likely take several (15+) minutes. Visualization parameters are found on lines 15-25 of main.jl Here are the default parameters and their meaning.
###
Total number of frames in each .gif. Large performance impact:
`nframes = 300`
###
Number of eigenvalues calculated for the tangential velocity profile. Medium performance impact:
`number_eigenvalues = 45`
###
Number of divisions for numeric integration. Massive performance and quality impact:
`integralresolution = 50000`
###
Number of points between R1 and R2 to visualize. Large performance impact:
`Rresolution = 350`
###
Framerate of animations. No performance impact:
`FPS = 30`
###
File names for each of the five animations. No performance impact:
`Vθ_3Dfilename = "Vθ_3Danimation.gif"`
`Vθ_2Dfilename = "Vθ_2Danimation.gif"`
`Vz_3Dfilename = "Vz_3Danimation.gif"`
`Vz_2Dfilename = "Vz_2Danimation.gif"`
`V_3Dfilename = "V_3Danimation.gif"`
###
Number of velocity profiles to show on the V_3Danimation.gif. Medium to low performance impact:
`nprofiles = 8`
