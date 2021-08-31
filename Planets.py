from astropy.time import Time
from astropy.coordinates import get_body_barycentric_posvel, solar_system_ephemeris
from astropy.constants import G
from spiceypy import sxform, mxvg
from poliastro import constants
from Particle import Particle
import numpy as np

def getPlanets(t=Time("2019-04-11 11:00:00", scale="tdb")):
    """
    Fetches the masses, positions and velocities of all the planets in the solar system at a arbitrary time, t.
    Returns this as a list where the elements in the list are instances of the Particle class.

    Parameters
    ----------
    t: class
        Instance of the class "astropy.time.core.Time".
        An arbitrary time in the required format for the astropy.time module.
    
    Returns
    -------
    Bodies: list
        A list of all the planets, in the form of Particle instances, in the solar system at the specified time (t) from the JPL.
    """
    Bodies=[]
    for planet in [x for x in solar_system_ephemeris.bodies if x != "earth-moon-barycenter"]:
        pos, vel=get_body_barycentric_posvel(planet, t, ephemeris="jpl")
        def metric(x):
            if x==pos:
                return "m"
            elif x==vel:
                return "m/s"
        statevec=[x.xyz[i].to(metric(x)).value for x in [pos, vel] for i in range(3)]
        statevececl=mxvg(sxform("J2000", "ECLIPJ2000", t.jd), statevec, 6, 6)
        planet=Particle(
            position=np.array([statevececl[0], statevececl[1], statevececl[2]]),
            velocity=np.array([statevececl[3], statevececl[4], statevececl[5]]),
            name=planet.capitalize(),
            mass=(getattr(constants, "GM_{}".format(planet)) / G).value
        )
        Bodies.append(planet)
    return Bodies