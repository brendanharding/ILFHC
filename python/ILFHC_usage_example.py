"""This script demonstrates some very basic usage of the 
Inertial Lift Force Helper Class (or ILFHC for short).
The class methods are themselves documented via 
```help(InertialLiftForceHelper)```.

Note: use of the helper class requires the following files 
to be located relative to the working directory:
 - 'curved_duct_lift_data/square_lift_data.npz'  
 - 'curved_duct_lift_data/rect_2x1_lift_data.npz' 
 - 'curved_duct_lift_data/rect_4x1_lift_data.npz' 

Please ensure you cite our JFM paper (https://doi.org/10.1017/jfm.2019.323) 
if you use this code/data. This code is provided under an MIT license 
(see https://opensource.org/licenses/MIT). However, I would appreciate it
if you contact me and to let me know if you use this code/data.
Please also don't hesitate to contact me if you have any questions/queries.

Brendan Harding, 2019."""

# Import the helper class
from ILFHC import InertialLiftForceHelper

# Initialise an instance for a neutrally buoyant
# particle with radius a=0.10 within a duct having 
# bend radius R=80.0 and a square cross-section.
ILFH = InertialLiftForceHelper(0.10,80.0,'square')

# Print an estimate of the migration force at the 
# centre of the cross-section (includes both the 
# inertial lift and secondary flow contributions)
print(ILFH.migration_force(0.0,0.0))

# Generate a plot of the migration force
ILFH.plot_migration_force()

# Print the documentation for the class
# (including a description of the methods/functions available)
# (type q to quit in the terminal)
help(InertialLiftForceHelper) # or help(ILFH)
