from Particle import Particle
from Simulation import Simulation
from Planets import getPlanets
from Analysis import detAccuracy
from astropy.time import Time, TimeDelta
import numpy as np
import copy

def analyticalTest():
    #defining the initial conditions
    Earth=Particle(
    position=np.zeros(3),
    velocity=np.zeros(3),
    name="Earth",
    mass=5.972e24
    )
    Sat=Particle(
    position=np.array([6371e3+36000e3,0,0]),
    velocity=np.array([0,(Earth.G*Earth.mass/(6371e3+36000e3))**0.5,0]),
    acceleration=np.array([-0.22201795, 0, 0]),
    name="Satellite",
    mass=3500
    )
    #instance of the simulation with the initial conditions defined above, the time step, time interval and method can easily be changed for testing purposes
    analyticalSim=Simulation(
    name="Analytical_Test",
    Bodies=[Earth, Sat],
    timeIntervals=10850,
    deltaT=8.,
    method=4,
    testLM=1,
    testAM=1
    )
    #solutions to the analytical simulation are just its initial conditions as the total time is equal to the orbital time period
    sol=[copy.deepcopy(Earth), copy.deepcopy(Sat)]
    #runs the simulation, calculates its inaccuracies and prints all the relevant information in a readable format
    analyticalSim.run()
    analyticalAcc=np.array([100*abs(1-analyticalSim.Bodies[x].__getattribute__(variable)[i]/sol[x].__getattribute__(variable)[i])
            for x in range(len(analyticalSim.Bodies))
            for variable in ["position", "velocity", "acceleration"]
            for i in range(3)
            if sol[x].__getattribute__(variable)[i] != 0
        ])
    print(analyticalSim, f"\nThe Satellite's Position, Velocity and Acceleration after {test.timeIntervals*test.deltaT} seconds are:")
    print(Sat)
    print(f"Max Inaccuracy: {np.amax(analyticalAcc)}%\nMean Inaccuracy: {np.mean(analyticalAcc)}%\nChange in Linear Momentum: {analyticalSim.LMC}\nChange in Angular Momentum: {analyticalSim.AMC}")

#initial conditions for running a simulation and comparing the results to the JPL
t1=Time("2019-04-11 11:00:00", scale="tdb")
test=Simulation(
    name="Test",
    Bodies=getPlanets(t1),
    timeIntervals=5000, 
    deltaT=8., 
    method=1,
    testLM=1,
    testAM=1
)
#runs the simulation instance test and calculates the accuracies of my simulation when comparing it to the JPL then prints the relevant information in a readable format
def jplTest(t1, test):
    #inputs: starting time, timeIntervals, deltaT, method, testLM, testAM
    t2=t1+TimeDelta(test.timeIntervals*test.deltaT, format="sec")
    test, accuracies=detAccuracy(t2, test)
    print(test, f"\nThe Bodies Positions, Velocities and Accelerations after {test.timeIntervals*test.deltaT} seconds are:")
    for body in test.Bodies:
        print(body)
    print(f"Max Inaccuracy: {np.amax(accuracies)}%\nMean Inaccuracy: {np.mean(accuracies)}%\nChange in Linear Momentum: {test.LMC}\nChange in Angular Momentum: {test.AMC}")

#just some examples of how the code can be ran
if __name__ == "__main__":
    #analyticalTest()
    jplTest(t1, test)
