import numpy as np
import copy
import matplotlib.pyplot as plt
from Planets import getPlanets

class Simulation:
    """
    Class to setup a simulation instance with an associated list of bodies.

    Parameters
    ----------
    name: str
        Name of the simulation.
    Bodies: list
        An arbitrary list of bodies that can be simulated. 
        Defaults to getting a list of all the planets in the solar system from getPlanets()
    method: int
        An integer defining the numerical approximation method to be used: 1 is the Euler method,
        2 is the Euler-Cromer method, 3 is the Euler-Richardson method, 4 is the Verlet method.
        The default is to use the Verlet method.
    testLM: int
        An integer defining whether the linear momentum before, 
        after and the change in linear momentum will be calculated using the calcLinearMomentum method.
    testAM: int
        An integer defining whether the angular momentum before, 
        after and the change in angular momentum will be calculated using the calcAngularMomentum method.
    timeIntervals: int
        The number of time steps that will be ran if the simulation is ran.
    deltaT: float
        The value of the time step that will be used if the simulation is ran.
    """
    def __init__(self, name="Verlet", Bodies=getPlanets(), timeIntervals=14400, deltaT=60., method=4, testLM=0, testAM=0):
        self.name=name
        self.Bodies=Bodies
        self.timeIntervals=timeIntervals
        self.deltaT=deltaT
        self.setMethod(method)
        self.testLM=testLM
        self.testAM=testAM
        self.LMB=self.calcLinearMomentum()
        self.AMB=self.calcAngularMomentum()
    
    def __str__(self):
        #Sets appropriate values for method being simulated and what conservation tests are being performed
        if self.testLM==0:
            LM="No"
        else:
            LM="Yes"
        if self.testAM==0:
            AM="No"
        else:
            AM="Yes"
        #Uses list comprehension as well as the input data to give all the important details about the simulation being run when its printed
        return "Name: {0}\nMethod: {1}\nBodies: {2}\nLinear Momentum Tested?: {3}\nAngular Momentum Tested?: {4}\nTime: {5}s over {6} intervals with a difference {7}s betweeb them.".format(
            self.name, self.mName, ", ".join([body.name for body in self.Bodies]), LM, AM,
            self.deltaT*self.timeIntervals, self.timeIntervals, self.deltaT
        )

    def setMethod(self, method):
        if method==1:
            self.run=self.methodE
            self.mName="Euler"
        elif method==2:
            self.run=self.methodEC
            self.mName="Euler-Cromer"
        elif method==3:
            self.run=self.methodER
            self.mName="Euler-Richardson"
        elif method==4:
            self.run=self.methodV
            self.mName="Verlet"
        else:
            raise ValueError("Invalid method inputted, method must be 1-4.")
    
    def methodE(self):
        Data=[]
        for x in range(1,self.timeIntervals+1):
            #first calculates the gravitational acceleration for all bodies, then calculates their position and velocity using the Euler method
            for body in self.Bodies:
                body.acceleration=np.sum([body.updateGravitationalAcceleration(self.Bodies[i]) for i in range(len(self.Bodies)) if self.Bodies[i] !=body], axis=0)
            for body in self.Bodies:
                body.updateE(self.deltaT)
            #stores each body being simulated at the respective time to a list every 100 time intervals
            if (x-1)%100 ==0:
                tempList=[copy.deepcopy(member) for member in self.Bodies]
                tempList.insert(0, x*self.deltaT)
                Data.append(tempList)
        #calculates change in linear and angular momentum and saves the data to a file
        self.LMC=abs((np.linalg.norm(self.calcLinearMomentum())-np.linalg.norm(self.LMB)))
        self.AMC=abs((np.linalg.norm(self.calcAngularMomentum())-np.linalg.norm(self.AMB)))
        np.save(self.name, Data, allow_pickle=True)

    def methodEC(self):
        Data=[]
        for x in range(1,self.timeIntervals+1):
            #first calculates the gravitational acceleration for all bodies, then calculates their position and velocity using the Euler-Cromer method
            for body in self.Bodies:
                body.acceleration=np.sum([body.updateGravitationalAcceleration(self.Bodies[i]) for i in range(len(self.Bodies)) if self.Bodies[i] !=body], axis=0)
            for body in self.Bodies:
                body.updateEC(self.deltaT)
            #stores each body being simulated at the respective time to a list every 100 time intervals
            if (x-1)%100 ==0:
                tempList=[copy.deepcopy(member) for member in self.Bodies]
                tempList.insert(0, x*self.deltaT)
                Data.append(tempList)
        #calculates change in linear and angular momentum and saves the data to a file 
        self.LMC=abs((np.linalg.norm(self.calcLinearMomentum())-np.linalg.norm(self.LMB)))
        self.AMC=abs((np.linalg.norm(self.calcAngularMomentum())-np.linalg.norm(self.AMB)))
        np.save(self.name, Data, allow_pickle=True)
    
    def methodER(self):
        Data=[]
        for x in range(1,self.timeIntervals+1):
            #first calculates the gravitational acceleration for all bodies, then the position midpoint between the current position and the position after
            #one time step for all bodies.
            for body in self.Bodies:
                body.acceleration=np.sum([body.updateGravitationalAcceleration(self.Bodies[i]) for i in range(len(self.Bodies)) if self.Bodies[i] !=body], axis=0)
            for body in self.Bodies:
                body.pmid=body.position+0.5*self.deltaT*body.velocity
            #Calculates the acceleration mid point for one body,
            #then uses that to calculate the velocity and position of that body using the Euler-Richardson method,
            #then repeats this process for all bodies
            for body in self.Bodies:
                amid=np.zeros(3)
                for i in range(len(self.Bodies)):
                    if self.Bodies[i] !=body:
                        diff=body.pmid-self.Bodies[i].pmid
                        r=np.linalg.norm(diff)
                        amid+=((-body.G*self.Bodies[i].mass)/(r**2))*((diff)/(r))
                body.updateER(self.deltaT, amid)
            #stores each body being simulated at the respective time to a list every 100 time intervals
            if (x-1)%100 ==0:
                tempList=[copy.deepcopy(member) for member in self.Bodies]
                tempList.insert(0, x*self.deltaT)
                Data.append(tempList)
        #calculates change in linear and angular momentum and saves the data to a file
        self.LMC=abs((np.linalg.norm(self.calcLinearMomentum())-np.linalg.norm(self.LMB)))
        self.AMC=abs((np.linalg.norm(self.calcAngularMomentum())-np.linalg.norm(self.AMB)))
        np.save(self.name, Data, allow_pickle=True)

    def methodV(self):
        Data=[]
        for x in range(1,self.timeIntervals+1):
            #first calculates the gravitational acceleration for all bodies, then the position after one time step using the Euler-Cromer method, for every body
            for body in self.Bodies:
                body.acceleration=np.sum([body.updateGravitationalAcceleration(self.Bodies[i]) for i in range(len(self.Bodies)) if self.Bodies[i] !=body], axis=0)
            for body in self.Bodies:
                body.pV=body.position+(body.velocity+body.acceleration*self.deltaT)*self.deltaT
            #Calculates the acceleration after one time step for one body, 
            #then uses that to calculate the velocity and position of that body using the Verlet method, then repeats this process for all bodies
            for body in self.Bodies:
                endAcceleration=np.zeros(3)
                for i in range(len(self.Bodies)):
                    if self.Bodies[i] !=body:
                        diff=body.pV-self.Bodies[i].pV
                        r=np.linalg.norm(diff)
                        endAcceleration+=((-body.G*self.Bodies[i].mass)/(r**2))*((diff)/(r))
                body.updateV(self.deltaT, endAcceleration)
            #stores each body being simulated at the respective time to a list every 100 time intervals
            if (x-1)%100 ==0:
                tempList=[copy.deepcopy(member) for member in self.Bodies]
                tempList.insert(0, x*self.deltaT)
                Data.append(tempList)
        #calculates change in linear and angular momentum and saves the data to a file
        self.LMC=abs((np.linalg.norm(self.calcLinearMomentum())-np.linalg.norm(self.LMB)))
        self.AMC=abs((np.linalg.norm(self.calcAngularMomentum())-np.linalg.norm(self.AMB)))
        np.save(self.name, Data, allow_pickle=True)

    def calcLinearMomentum(self):
        if self.testLM==1:
            return np.sum([body.linearMomentum() for body in self.Bodies], axis=0)
        else:
            return np.zeros(3)

    def calcAngularMomentum(self):
        if self.testAM==1:
            return np.sum([body.angularMomentum() for body in self.Bodies], axis=0)
        else:
            return np.zeros(3)

    #plots all the bodies in the simulation at their current position
    def plot(self):
        DataIn=np.load("{}.npy".format(self.name), allow_pickle=True)
        ax, length=plt.subplot(), range(len(DataIn))
        for x in range(1,len(DataIn[0])):
            ax.plot([DataIn[i][x].position[0] for i in length], 
                [DataIn[i][x].position[1] for i in length],
                label=DataIn[0][x].name
            )
        ax.set_xlabel(r"$x$ position [m]")
        ax.set_ylabel(r"$y$ position [m]")
        ax.legend()
        plt.savefig("{}.png".format(self.name), dpi=250)
        plt.show()