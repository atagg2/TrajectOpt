using TrajectOpt
using Plots
using FLOWMath
using DelimitedFiles
using CCBlade

#general system
m = 1.36
I = .0111
X = .5
L = 4.0
rho = 1.225
mu = 1.81e-5
g = 9.81
cog = [-.5, 0.0, 0.0]

# wing
xleWing = [0.0, 0.2]
yleWing = [0.0, 5.0]
zleWing = [0.0, 0.0]
cWing = [1.0, 0.6]
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
xleElevator = [L + cTail[1], L + cTail[1]]
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

#state and input
# x = [
# 10       # Vinf
# 0       # flight path angle
# 0       # body angle derivative
# .1       # body angle
# 0       # x position
# 0       # y position
# ]   
# u = [
# 0       # thrust
# 0     # elevator deflection
# ]

x = [216.74617936267427
-4.892643981926575e-10
-1.9617151912584933e-20
 0.0005885655318472827
 0.4043820707410358
 1.2600833289370368]

u = [2.073892654039904e-5
-0.00011684186012880983]

alphas = -5:.5:15
xForceTrajectOpt = zeros(length(alphas))
yForceTrajectOpt = zeros(length(alphas))
momentTrajectOpt = zeros(length(alphas))

for i in 1:length(alphas)
    x[4] = alphas[i]
    F, M = plane.forces(x,u)
    xForceTrajectOpt[i] = F[1]
    yForceTrajectOpt[i] = F[2]
    momentTrajectOpt[i] = M
end

A = readdlm("files/xflr_full_results.txt")

Cd = A[:,6]
Cl = A[:,3]
Cm = A[:,9]

Sref = 1
cref = 1
D = Cd .*.5*rho*x[1]^2*Sref
L = Cl .*.5*rho*x[1]^2*Sref
momentXflr = Cm .*.5*rho*x[1]^2*Sref*cref

xForceXflr = u[1]*cosd.(alphas) .- D
yForceXflr = u[1]*sind.(alphas) .+ L .- m*g

plot(alphas, xForceTrajectOpt, xlabel = "Angle of Attack (degrees)", ylabel = "X Force (N)", label = "TrajectOpt")
# scatter!(alphas, xForceXflr, label = "Xflr")

# plot(alphas, yForceTrajectOpt, xlabel = "Angle of Attack (degrees)", ylabel = "Y Force (N)", label = "TrajectOpt")
# scatter!(alphas, yForceXflr, label = "Xflr")

# plot(alphas, momentTrajectOpt, xlabel = "Angle of Attack (degrees)", ylabel = "Moment (Nm)", label = "TrajectOpt")
# scatter!(alphas, momentXflr, label = "Xflr")