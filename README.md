# interplanetary_trajectory

MATLAB tool which can be used to define interplanetary missions between any two planets of the solar system, using the method of patched conics. The software allows the user to specify the desired planets and dates for departure and arrival and will output the results in the console together with a 3D graphical representation of the heliocentric trajectory.
The method of patched conics divides the mission into three phases: the hyperbolic departure phase, the heliocentric trajectory phase and the arrival phase. It is assumed that the spacecraft departs from a circular parking orbit around the departure planet and embarks on a heliocentric trajectory up to the arrival planet where it will go on a circular capture orbit. Inside the spheres of influence of the planets the only force acting on the spacecraft is assumed to be the gravitational pull of the planet, and outside of them the gravitational pull of the sun, any other forces or perturbations being neglected.


[![View Interplanetary Mission Design on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://ch.mathworks.com/matlabcentral/fileexchange/66192-interplanetary-mission-design)
