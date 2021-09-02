from Planets import getPlanets
from Simulation import Simulation
from astropy.time import Time
import numpy as np

def detAccuracy(t=Time("2019-04-21 11:00:00", scale="tdb"), test=Simulation()):
    """
    Returns the ran instance of Simulation given and percentage differences between every planet, 
    after the test has been ran, and what the planets actual locations were according to the JPL.

    Parameters
    ----------
    t: class
        Instance of the class "astropy.time.core.Time".
        The exact time which the simulation will be at once it has finished running in the required format for the astropy.time module.
    test: class
        Instance of the class "Simulation" that has not been ran before.
    
    Returns
    -------
    test: instance of the class "Simulation"
        The instance has been ran over a given number of time intervals for a given time step using the method in the instance of the class.
    ndarray
        Contains percentage differences between every planet, after the test has been ran, and what the planets actual locations were according to the JPL.
    """
    test.run()
    Bodies=getPlanets(t)
    return test, np.array([100*abs(1-test.Bodies[x].__getattribute__(variable)[i]/Bodies[x].__getattribute__(variable)[i])
        for x in range(len(Bodies))
        for variable in ["position", "velocity"]
        for i in range(3)
    ])
