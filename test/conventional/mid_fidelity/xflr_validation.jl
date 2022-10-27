using TrajectOpt
using Plots
using FLOWMath
using DelimitedFiles

#general system
m = 1.36
I = .0111
X = .5
L = 4.0
rho = 1.225
mu = 1.81e-5
g = 9.81

# wing
xleWing = [0.0, 0.2]
yleWing = [0.0, 5.0]
zleWing = [0.0, 0.0]
cWing = [1.0, 0.6]
twistWing = [0.0, 0.0]
phi = [0.0, 0.0]
camberWing = fill((xc) -> 0, 2) # camberline function for each section

# horizontal stabilizer
xleTail = [0.0, 0.14]
yleTail = [0.0, 1.25]
zleTail = [0.0, 0.0]
cTail = [0.7, 0.42]
twistTail = [0.0, 0.0]
camberTail = fill((xc) -> 0, 2) # camberline function for each section
cElevatorFraction = .25
bElevatorFraction = .9

#create model
parameters = ConventionalMidFidel(m, I, X, L, xleWing, yleWing, zleWing, cWing, twistWing,
  camberWing, xleTail, yleTail, zleTail, cTail, twistTail, camberTail, cElevatorFraction, 
  bElevatorFraction, rho, mu, g)
forces = conventional_forces_constructor(parameters)
plane = HighFidel(parameters, forces)

#state and input
x = [
10       # Vinf
0       # flight path angle
0       # body angle derivative
.1       # body angle
0       # x position
0       # y position
]   
u = [
5       # thrust
0     # elevator deflection
]

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

A = readdlm("T1-10 m_s-VLM2.txt")

Cd = A[:,6]
Cl = A[:,3]
Cm = A[:,9]

Sref = 8
cref = .817
D = Cd .*.5*rho*x[1]^2*Sref
L = Cl .*.5*rho*x[1]^2*Sref
momentXflr = Cm .*.5*rho*x[1]^2*Sref*cref

xForceXflr = u[1]*cosd.(alphas) .- D
yForceXflr = u[1]*sind.(alphas) .+ L .- m*g

plot(alphas, xForceTrajectOpt, xlabel = "Angle of Attack (degrees)", ylabel = "X Force (N)", label = "TrajectOpt")
scatter!(alphas, xForceXflr, label = "Xflr")

plot(alphas, yForceTrajectOpt, xlabel = "Angle of Attack (degrees)", ylabel = "Y Force (N)", label = "TrajectOpt")
scatter!(alphas, yForceXflr, label = "Xflr")

# plot(alphas, momentTrajectOpt, xlabel = "Angle of Attack (degrees)", ylabel = "Moment (Nm)", label = "TrajectOpt")
# scatter!(alphas, momentXflr, label = "Xflr")