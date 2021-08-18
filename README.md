Run my simulation using the Main.py file and run the functions in that file to run the code.

Main.py: File used to run the project. File where you can change the initial conditions of the simulation, what tests are run, what method is used and various other things of how the simulation can be ran. Prints results of the simulation ran in a readable format.
Particle.py: Class to setup a body instance with an associated name, velocity, acceleration, mass and position. Contains some basic methods that can update/calculate these attributes of the class.
Simulation.py: File containing a class to setup a simulation instance with an associated list od bodies. Contains the main bulk of the code where a lot of the calculations take place. Also calculates changes in linear and angular momentum as well as can plot the simulation.
Planets.py: Fetches the masses, positions and velocities of all the planets in the solar system at a arbitrary time. Returns this as a list where the elements in the list are instances of the Particle class.
Analysis.py: Returns the ran simulation class instance and calculates the percentage differences between my simulation and where the planets actually were according to the JPL.
