#Project by Nachiket Karve and Kaustav Kashyap Das, 2019

import numpy as np
import matplotlib.pyplot as plt

G = 6.67*10**(-11)  # Gravitational Constant
mEarth = 5.972*10**(24) # Mass of Earth
rEarth = 6378100  # Radius of Earth
deltaT = 0.01*(rEarth**3/(G*mEarth))**0.5  # Timestep

# Stores the position and velocity information of the body
class body:
    
    # Initializes the object
    def __init__(self, x, y, velX, velY, mass, time):
        self.x = x
        self.y = y
        self.velX = velX
        self.velY = velY
        self.mass = mass
        self.time = time
        self.r = (x**2 + y**2)**(0.5)
        self.theta = np.angle(complex(x, y))
    
    # Updates the position and velocity of the body after time 'delta'
    # 'fx' and 'fy' are the forces in the x and y directions
    # 'boostX' and 'boostY' are arrays which specifies the values of the boosts to be given
    # 'start' and 'end' arrays specify the start and end times for each boost
    def updateBody(self, fx, fy, delta, boostX, boostY, start, end):
        prevBody = self
        self.x = self.x + self.velX*delta + 0.5*fx(self, boostX, start, end)*(delta)**2/(self.mass)
        self.y = self.y + self.velY*delta + 0.5*fy(self, boostY, start, end)*(delta)**2/(self.mass)
        self.velX = self.velX + 0.5*(fx(self, boostX, start, end) + fx(prevBody, boostX, start, end))*delta/(self.mass)
        self.velY = self.velY + 0.5*(fy(self, boostY, start, end) + fy(prevBody, boostY, start, end))*delta/(self.mass)
        self.time = self.time + delta
        self.r = (self.x**2 + self.y**2)**(0.5)
        self.theta = np.angle(complex(self.x, self.y))

# Returns the force in the x direction on the body 'sat' by all the other bodies
def fXSat(sat, boost, start, end):
    giveBoost = 0
    for i in range(len(boost)):
        if sat.time <= end[i] and sat.time >= start[i]:
            giveBoost = boost[i]
            break
    force = 0
    for otherBody in bodies:
        if otherBody != sat:
            force = force - G*otherBody.mass*sat.mass*(sat.x - otherBody.x)/((sat.x - otherBody.x)**2 + (sat.y - otherBody.y)**2)**(3.0/2)
    return force + giveBoost

# Returns the force in the y direction on the body 'sat' by all the other bodies
def fYSat(sat, boost, start, end):
    giveBoost = 0
    for i in range(len(boost)):
        if sat.time <= end[i] and sat.time >= start[i]:
            giveBoost = boost[i]
            break
    force = 0
    for otherBody in bodies:
        if otherBody != sat:
            force = force - G*otherBody.mass*sat.mass*(sat.y - otherBody.y)/((sat.x - otherBody.x)**2 + (sat.y - otherBody.y)**2)**(3.0/2)
    return force + giveBoost

# Returns the cosine of the angle the velocity of 'sat' makes with the x-axis
def velXAng(sat):
    return sat.velX/(sat.velX**2 + sat.velY**2)**0.5

# Returns the sine of the angle the velocity of 'sat' makes with the x-axis
def velYAng(sat):
    return sat.velY/(sat.velX**2 + sat.velY**2)**0.5

# Defines the initial orbit of the satellite
satPerigee = 6714700.0
satApogee = 74715000.0
satSemiMajor = (satPerigee + satApogee)/2
satPerigeeVel = (G*mEarth*((2/satPerigee) - (1/satSemiMajor)))**0.5

# Defines the initial position of the Moon
MoonPos = 385000000.0
MoonVel = 1000.0
mMoon = 7.342*10**(22)
theta0 = -0.803589979409645

# Defines the initial velocity of the Earth
EarthVel = mMoon*MoonVel/mEarth

# Creates the bodies 'sat', 'Moon', and 'Earth'
sat = body(satPerigee, 0.0, 0.0, satPerigeeVel + EarthVel, 1000, 0.0)
Moon = body(MoonPos*np.cos(theta0), MoonPos*np.sin(theta0), -MoonVel*np.sin(theta0), MoonVel*np.cos(theta0), mMoon, 0.0)
Earth = body(0.0, 0.0, EarthVel*np.sin(theta0), -EarthVel*np.cos(theta0), 5.972*10**(24), 0.0)

# Creates an array which stores all the bodies
bodies = [sat, Earth, Moon]

# Creates arrays which store the state of the satellite, Moon, and Earth
length = 300000
x1 = np.zeros(length)
y1 = np.zeros(length)
t1 = np.zeros(length)
x2 = np.zeros(length)
y2 = np.zeros(length)
t2 = np.zeros(length)
x3 = np.zeros(length)
y3 = np.zeros(length)
t3 = np.zeros(length)

# Defines the values of the boosts in Newtons
boostVal1 = 50
boostVal2 = 50
boostVal3 = 50
boostVal4 = 500
boostVal5 = -50
boostVal6 = -50

# Defines the times at which the boosts are provided in seconds
boostStartTime1 = 80957.7444938
boostStartTime2 = 385483.47625
boostStartTime3 = 1050594.40474
boostStartTime4 = 1412382.1439953118
boostStartTime5 = 1540707.1502211716
boostStartTime6 = 1595588.284959275

# Defines the duration of the boosts in seconds
boostTime1 = 1200*deltaT
boostTime2 = 330*deltaT
boostTime3 = 115*deltaT
boostTime4 = 350*deltaT
boostTime5 = 100*deltaT
boostTime6 = 200*deltaT

# Iterate over time
for i in range(length):
    x1[i] = sat.x
    y1[i] = sat.y
    t1[i] = sat.time
    
    x2[i] = Moon.x
    y2[i] = Moon.y
    t2[i] = Moon.time
    
    x3[i] = Earth.x
    y3[i] = Earth.y
    t3[i] = Earth.time
    
    sat.updateBody(fXSat, fYSat, deltaT, [boostVal1*velXAng(sat), boostVal2*velXAng(sat), boostVal3*velXAng(sat), 0, boostVal5*velXAng(sat), boostVal6*velXAng(sat)], [boostVal1*velYAng(sat), boostVal2*velYAng(sat), boostVal3*velYAng(sat), -boostVal4, boostVal5*velYAng(sat), boostVal6*velYAng(sat)], [boostStartTime1 - boostTime1/2.0, boostStartTime2 - boostTime2/2.0, boostStartTime3 - boostTime3/2.0, boostStartTime4, boostStartTime5 - boostTime5/2.0, boostStartTime6 - boostTime6/2.0], [boostStartTime1 + boostTime1/2.0, boostStartTime2 + boostTime2/2.0, boostStartTime3 + boostTime3/2.0, boostStartTime4 + boostTime4, boostStartTime5 + boostTime5/2.0, boostStartTime6 + boostTime6/2.0])
    Moon.updateBody(fXSat, fYSat, deltaT, [0], [0], [0], [0])
    Earth.updateBody(fXSat, fYSat, deltaT, [0], [0], [0], [0])

# Plot the trajectories of all the bodies
# For making the video, we exported the arrays to MATLAB.
plt.plot(x1, y1, x2, y2, x3, y3)
