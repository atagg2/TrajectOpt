using TrajectOpt
using Plots
using FLOWMath
using DelimitedFiles
#========Lower Fidelity System==========#
# #physical parameters
# m = 1.36
# I = .0111
# X = .5
# L = 4.0     #tail moment arm
# SWing = 8.0
# bWing = 10.0

# STail = 1.4
# bTail = 2.5
# rho = 1.225
# mu = 1.81e-5
# g = 9.81


#physical system
m = 1.36
I = .0111
cog = [.5, 0.0, 0.0]
SWing = 8.0
bWing = 10.0
rWing = [0.0, 0.0, 0.0]
STail = 1.4
bTail = 2.5
rTail = [4.0, 0.0, 0.0]
rRotor = [0.0, 0.0, 0.0]
rho = 1.225
mu = 1.81e-5
g = 9.81



# planeParameters = ConventionalLowFidel(m, I, X, L, SWing, bWing, STail, bTail, rho, mu, g)

#polar
alphas = -5:.5:15
Res = [342000, 1026000]
Cds = [
    0.005236	0.005236
    0.004246	0.004246
    0.003359	0.003359
    0.002574	0.002574
    0.001893	0.001893
    0.001315	0.001315
    0.000842	0.000842
    0.000474	0.000474
    0.000211	0.000211
    5.3e-5	5.3e-5
    0.0	0.0
    5.3e-5	5.3e-5
    0.000211	0.000211
    0.000474	0.000474
    0.000842	0.000842
    0.001315	0.001315
    0.001893	0.001893
    0.002574	0.002574
    0.003359	0.003359
    0.004246	0.004246
    0.005236	0.005236
    0.006327	0.006327
    0.007518	0.007518
    0.00881	0.00881
    0.0102	0.0102
    0.011687	0.011687
    0.013271	0.013271
    0.01495	0.01495
    0.016723	0.016723
    0.018589	0.018589
    0.020546	0.020546
    0.022593	0.022593
    0.024728	0.024728
    0.026949	0.026949
    0.029255	0.029255
    0.031644	0.031644
    0.034115	0.034115
    0.036664	0.036664
    0.039291	0.039291
    0.041993	0.041993
    0.044769	0.044769
    
]

Cls = [
    -0.459485	-0.459485
    -0.413715	-0.413715
    -0.367888	-0.367888
    -0.322012	-0.322012
    -0.276091	-0.276091
    -0.230133	-0.230133
    -0.184144	-0.184144
    -0.13813	-0.13813
    -0.092097	-0.092097
    -0.046052	-0.046052
    0.0	0.0
    0.046052	0.046052
    0.092097	0.092097
    0.13813	0.13813
    0.184144	0.184144
    0.230133	0.230133
    0.276091	0.276091
    0.322012	0.322012
    0.367888	0.367888
    0.413715	0.413715
    0.459485	0.459485
    0.505194	0.505194
    0.550834	0.550834
    0.596399	0.596399
    0.641884	0.641884
    0.687282	0.687282
    0.732588	0.732588
    0.777794	0.777794
    0.822897	0.822897
    0.867888	0.867888
    0.912764	0.912764
    0.957517	0.957517
    1.002142	1.002142
    1.046634	1.046634
    1.090986	1.090986
    1.135193	1.135193
    1.179249	1.179249
    1.22315	1.22315
    1.266888	1.266888
    1.31046	1.31046
    1.35386	1.35386
]

Cms = [
    -0.114238	-0.114238
    -0.102914	-0.102914
    -0.091558	-0.091558
    -0.080174	-0.080174
    -0.068766	-0.068766
    -0.057337	-0.057337
    -0.045891	-0.045891
    -0.03443	-0.03443
    -0.022959	-0.022959
    -0.011481	-0.011481
    0.0	0.0
    0.011481	0.011481
    0.022959	0.022959
    0.03443	0.03443
    0.045891	0.045891
    0.057337	0.057337
    0.068766	0.068766
    0.080174	0.080174
    0.091558	0.091558
    0.102914	0.102914
    0.114238	0.114238
    0.125528	0.125528
    0.136779	0.136779
    0.147989	0.147989
    0.159154	0.159154
    0.17027	0.17027
    0.181334	0.181334
    0.192343	0.192343
    0.203294	0.203294
    0.214182	0.214182
    0.225005	0.225005
    0.23576	0.23576
    0.246443	0.246443
    0.257051	0.257051
    0.267581	0.267581
    0.278029	0.278029
    0.288392	0.288392
    0.298668	0.298668
    0.308852	0.308852
    0.318943	0.318943
    0.328936	0.328936
]

wing_polar = polar_constructor(Cds, Cls, Cms, alphas, Res)


Cds = [
    0.008234	0.008234
    0.006677	0.006677
    0.005282	0.005282
    0.004048	0.004048
    0.002976	0.002976
    0.002068	0.002068
    0.001324	0.001324
    0.000745	0.000745
    0.000331	0.000331
    8.3e-5	    8.3e-5
    0.0	        0.0
    8.3e-5	    8.3e-5
    0.000331	0.000331
    0.000745	0.000745
    0.001324	0.001324
    0.002068	0.002068
    0.002976	0.002976
    0.004048	0.004048
    0.005282	0.005282
    0.006677	0.006677
    0.008234	0.008234
    0.00995	    0.00995
    0.011824	0.011824
    0.013854	0.013854
    0.01604	    0.01604
    0.018379	0.018379
    0.02087	    0.02087
    0.023511	0.023511
    0.0263	    0.0263
    0.029234	0.029234
    0.032312	0.032312
    0.03553	    0.03553
    0.038888	0.038888
    0.042381	0.042381
    0.046008	0.046008
    0.049765	0.049765
    0.05365	    0.05365
    0.057659	0.057659
    0.06179	    0.06179
    0.06604	    0.06604
    0.070405	0.070405
]

Cls = [
    -0.353425	-0.353425
-0.318282	-0.318282
-0.283076	-0.283076
-0.247815	-0.247815
-0.212504	-0.212504
-0.177151	-0.177151
-0.141763	-0.141763
-0.106347	-0.106347
-0.07091	-0.07091
-0.035458	-0.035458
0.0	        0.0
0.035458	0.035458
0.07091	    0.07091
0.106347	0.106347
0.141763	0.141763
0.177151	0.177151
0.212504	0.212504
0.247815	0.247815
0.283076	0.283076
0.318282	0.318282
0.353425	0.353425
0.388498	0.388498
0.423494	0.423494
0.458407	0.458407
0.49323	    0.49323
0.527955	0.527955
0.562578	0.562578
0.597089	0.597089
0.631484	0.631484
0.665756	0.665756
0.699898	0.699898
0.733904	0.733904
0.767767	0.767767
0.801482	0.801482
0.835041	0.835041
0.86844	    0.86844
0.901672	0.901672
0.93473	    0.93473
0.96761	    0.96761
1.000305	1.000305
1.032811	1.032811
]

Cms = [
    -0.064541	-0.064541
    -0.058143	-0.058143
    -0.051728	-0.051728
    -0.045296	-0.045296
    -0.038851	-0.038851
    -0.032394	-0.032394
    -0.025927	-0.025927
    -0.019452	-0.019452
    -0.012971	-0.012971
    -0.006487	-0.006487
    0.0	         0.0
    0.006487	0.006487
    0.012971	0.012971
    0.019452	0.019452
    0.025927	0.025927
    0.032394	0.032394
    0.038851	0.038851
    0.045296	0.045296
    0.051728	0.051728
    0.058143	0.058143
    0.064541	0.064541
    0.07092	    0.07092
    0.077276	0.077276
    0.08361	    0.08361
    0.089917	0.089917
    0.096198	0.096198
    0.102449	0.102449
    0.108668	0.108668
    0.114855	0.114855
    0.121007	0.121007
    0.127122	0.127122
    0.133198	0.133198
    0.139233	0.139233
    0.145227	0.145227
    0.151175	0.151175
    0.157078	0.157078
    0.162933	0.162933
    0.168739	0.168739
    0.174493	0.174493
    0.180193	0.180193
    0.185839	0.185839    
]

tail_polar = polar_constructor(Cds, Cls, Cms, alphas, Res)


environment = Environment(rho, mu, g)
inertia = Inertia(m, I, cog)
wing_polar = polar_constructor(Cds, Cls, Cms, alphas, Res)
wing = SimpleSurface(SWing, bWing, rWing, wing_polar)
tail_polar = polar_constructor(Cds, Cls, Cms, alphas, Res)
tail = SimpleSurface(STail, bTail, rTail, tail_polar)
surfaces = [wing, tail]
rotors = [SimpleRotor(rRotor)]
parameters = Parameters(environment, inertia, surfaces, rotors)
forces = forces_conventional_low_fidel(parameters)
planeLowFidel = Model(parameters, forces)

# tail_polar_function = polar_constructor(Cds, Cls, Cms, alphas, Res)
# planeForces = conventional_forces_constructor(wing_polar_function, tail_polar_function, planeParameters)
# planeLowFidel = LowFidel(planeParameters, planeForces)

#========Higher Fidelity System==========#
#general system
m = 1.36
I = .0111
X = .5
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
rRotor = [0.0, 0.0, 0.0]

# #create model
# parameters = ConventionalMidFidel(m, I, X, L, xleWing, yleWing, zleWing, cWing, twistWing,
#   camberWing, xleTail, yleTail, zleTail, cTail, twistTail, camberTail, cElevatorFraction, 
#   bElevatorFraction, rho, mu, g)
# forces = conventional_forces_constructor(parameters)
# planeMidFidel = HighFidel(parameters, forces)

#create model
environment = Environment(rho, mu, g)
inertia = Inertia(m, I, cog)
wing = VLMSurface(xleWing, yleWing, zleWing, cWing, twistWing, camberWing)
tail = VLMSurface(xleTail, yleTail, zleTail, cTail, twistTail, camberTail)
elevator = VLMSurface(xleElevator, yleElevator, zleElevator, cElevator, twistElevator, camberElevator)
surfaces = [wing, tail, elevator]
rotors = [SimpleRotor(rRotor)]
parameters = Parameters(environment, inertia, surfaces, rotors)
forces = forces_conventional_mid_fidel(parameters)
planeMidFidel = Model(parameters, forces)

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

xForceLowFidel = zeros(length(alphas))
yForceLowFidel = zeros(length(alphas))
momentLowFidel = zeros(length(alphas))

xForceHighFidel = zeros(length(alphas))
yForceHighFidel = zeros(length(alphas))
momentHighFidel = zeros(length(alphas))

for i in 1:length(alphas)
    x[4] = alphas[i]
    
    F, M = planeLowFidel.forces(x,u)
    xForceLowFidel[i] = F[1]
    yForceLowFidel[i] = F[2]
    momentLowFidel[i] = M

    F, M = planeMidFidel.forces(x,u)
    xForceHighFidel[i] = F[1]
    yForceHighFidel[i] = F[2]
    momentHighFidel[i] = M
end

plot(alphas, xForceLowFidel, xlabel = "Angle of Attack (degrees)", ylabel = "X Force (N)", label = "Low Fidelity", markershape = :circle)
plot!(alphas, xForceHighFidel, label = "Mid Fidelity", markershape = :circle)

plot(alphas, yForceLowFidel, xlabel = "Angle of Attack (degrees)", ylabel = "Y Force (N)", label = "Low Fidelity", markershape = :circle)
plot!(alphas, yForceHighFidel, label = "Mid Fidelity", markershape = :circle, legend = :topleft)

plot(alphas, momentLowFidel, xlabel = "Angle of Attack (degrees)", ylabel = "Moment (Nm)", label = "Low Fidelity", markershape = :circle)
plot!(alphas, momentHighFidel, label = "Mid Fidelity", markershape = :circle)