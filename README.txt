Shaped trajectory solver for finding approximately optimal trajectory for low thrust 
Earth-Mars interplanetary mission

User guide:
 - Run code from traj.m to get individual trajectories 
     - Set arrival and departure time t.arrive and t.depart 
     - Set target planets for Earth -> Mars or Mars -> Earth (P1 = E, P2 = M) or (P1 = M, P2 = E)
     - Run and read off relevant values
     - If shaped trajectory cannot be found code returns error
 - Run code from porkchop.m to get delta V and required thrust maps for a range of dates
     - Set values over which to calculate delta V and thrust map 
     - Set target planets for Earth -> Mars or Mars -> Earth (P1 = E, P2 = M) or (P1 = M, P2 = E)
      - Run and read off relevant values
      - Trajectories which cannot be resolved are returned as NaN (error = true)
 - Run code from TOF_porkchop.m to get delta V and required thrust maps for a range of departure dates and TOF
     - Set departure time range and step
     - Set TOF range and step
     - Run and read off relevant values
      - Trajectories which cannot be resolved should return appropriate error message in error array

Changelog:
v1.0: First working instance of code 
v1.1: Improved user interface and allowed conversion from canonical units to SI units
v1.2: Added porkchop function and fixed bug with theta_f
v1.2.1: Instanced to keep track of changes
v1.3: Added TOF porkchop function and realistic thrust/power requirement plots for individual trajectories
v1.3.1: Instanced to keep track of changes
v2.0: Up-revisioned to prevent conflict implementing 3D solver

Planned updates:

Bug fixes:
 - Angle theta_f previously was limited to between 0<theta_f<pi meaning trajectories where 2*pi>theta_f>pi would not return a result
     - Fixed by adding new method taking into account direction of cross product of planet vectors, this is likely to breakdown if difference in orbital inclination is large
     - New issue with integrator not solving correctly
 - Issue with N_rev not being added to theta_f for multiple rev orbits
 - Fixed bug with spacecraft mass lost from rocket equation
 - Major bug fix with theta_f angle being incorrectly altered