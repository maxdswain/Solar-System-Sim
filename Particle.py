import numpy as np

class Particle:
    """
    Class to setup a body instance with an associated name, velocity, acceleration, mass and position.

    Parameters
    ----------
    name: str
        Name of the body.
    mass: float
        Mass of the body in kg.
    position: ndarray
        The position of the body in cartesian coordinates.
    velocity: ndarray
        The velocity of the body in cartesian coordinates.
    acceleration: ndarray
        The acceleration of the body in cartesian coordinates.
    """
    def __init__(
        self,
        position=np.zeros(3),
        velocity=np.zeros(3),
        acceleration=np.zeros(3),
        name='Ball',
        mass=1.0,
    ):
        self.name=name
        self.mass=mass
        self.position=position.copy().astype(np.float)
        self.velocity=velocity.copy().astype(np.float)
        self.acceleration=acceleration.copy().astype(np.float)
        self.G=6.6743e-11

    def __str__(self):
        return "Particle: {0}\n     Mass: {1:.3e}\n     Position: {2}\n     Velocity: {3}\n     Acceleration: {4}".format(
            self.name, self.mass,self.position, self.velocity, self.acceleration
        )

    #calculates the position and velocity of the body after one time step using the Euler method
    def updateE(self, deltaT):
        self.position+=self.velocity*deltaT
        self.velocity+=self.acceleration*deltaT

    #calculates the position and velocity of the body after one time step using the Euler-Cromer method
    def updateEC(self, deltaT):
        self.velocity+=self.acceleration*deltaT
        self.position+=self.velocity*deltaT

    #calculates the position and velocity of the body after one time step using the Euler-Richardson method
    def updateER(self, deltaT, amid):
        vmid=self.velocity+0.5*deltaT*self.acceleration
        self.velocity+=amid*deltaT
        self.position+=vmid*deltaT

    #calculates the position and velocity of the body after one time step using the Verlet method
    def updateV(self, deltaT, endAcceleration):
        self.position+=self.velocity*deltaT+0.5*self.acceleration*deltaT**2
        self.velocity+=0.5*deltaT*(endAcceleration+self.acceleration)

    #calculates the gravitational acceleration between this body and another arbitrary body
    def updateGravitationalAcceleration(self, body):
        diff=np.subtract(self.position, body.position)
        r=np.linalg.norm(diff)
        return ((-self.G*body.mass)/(r**2))*(diff/(r))

    def linearMomentum(self):
        return np.multiply(self.mass, self.velocity)

    def angularMomentum(self):
        return np.cross(self.position, np.multiply(self.mass, self.velocity))