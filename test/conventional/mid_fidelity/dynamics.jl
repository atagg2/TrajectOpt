using TrajectOpt
using Plots
using FLOWMath
using CCBlade

# general parameters
m = 1.36
I = .0111
cog = [.5, 0.0, 0.0]
L = 4.0
rho = 1.225
mu = 1.81e-5
g = 9.81

# wing
xleWing = [0.0, 0.0]
yleWing = [0.0, 5.0]
zleWing = [0.0, 0.0]
cWing = [.8, .8]
twistWing = [0.0, 0.0]
phi = [0.0, 0.0]
camberWing = fill((xc) -> 0, 2) # camberline function for each section

# horizontal stabilizer
xleTail = [L, L]
yleTail = [0.0, 1.25]
zleTail = [0.0, 0.0]
cTail = [0.56, 0.56]*.6
twistTail = [0.0, 0.0]
camberTail = fill((xc) -> 0, 2) # camberline function for each section

# elevator
xleElevator = [L, L]
yleElevator = [0.0, 1.25]
zleElevator = [0.0, 0.0]
cElevator = [0.56, 0.56]*.3
twistElevator = [0.0, 0.0]
camberElevator = fill((xc) -> 0, 2)

# rotor
pos = [0.0, 0.0, 0.0]
Rtip = 10/2.0 * 0.0254 
Rhub = 0.10*Rtip
n = 2  
propgeom = [
0.15   0.130   32.76
0.20   0.149   37.19
0.25   0.173   33.54
0.30   0.189   29.25
0.35   0.197   25.64
0.40   0.201   22.54
0.45   0.200   20.27
0.50   0.194   18.46
0.55   0.186   17.05
0.60   0.174   15.97
0.65   0.160   14.87
0.70   0.145   14.09
0.75   0.128   13.39
0.80   0.112   12.84
0.85   0.096   12.25
0.90   0.081   11.37
0.95   0.061   10.19
1.00   0.041   8.99
]
r = propgeom[:, 1] * Rtip
cRotor = propgeom[:, 2] * Rtip
twistRotor = propgeom[:, 3] * pi/180
af = AlphaAF("files/naca4412.dat")

#create model
environment = Environment(rho, mu, g)
inertia = Inertia(m, I, cog)
wing = VLMSurface(xleWing, yleWing, zleWing, cWing, twistWing, camberWing)
tail = VLMSurface(xleTail, yleTail, zleTail, cTail, twistTail, camberTail)
elevator = VLMSurface(xleElevator, yleElevator, zleElevator, cElevator, twistElevator, camberElevator)
surfaces = [wing, tail, elevator]
rotors = [CCBladeRotor(pos, Rtip, Rhub, n, r, cRotor, twistRotor, af)]
parameters = Parameters(environment, inertia, surfaces, rotors)
forces = forces_conventional_mid_fidel(parameters)
plane = Model(parameters, forces)

x0 = [25.0000000000548
   6.96392087655856e-15
  -5.999998907734639e-9
   0.007078442538802229
 227.29577828161197
 320.43033695848163]

u = [ 0.000250088702418036
-0.001405209044108578]

tSpan = [0 1]
t = range(0, stop = tSpan[2], length = 100)
deflections = u[2]*ones(length(t))
# deflections[40:48] .= -.5
uSpline = [Akima(t,u[1]*ones(length(t))), Akima(t, deflections)]

path = simulate(x0, uSpline, plane, tSpan)
plot_simulation(path, uSpline)
