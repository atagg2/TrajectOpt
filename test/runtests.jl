using TrajectOpt
using Test
using Plots

@testset "dynamics.jl" begin
    x = [10,90,0,90,0,0]
    m = 1
    I = .0111
    X = .0086
    s = .25
    S = .05
    b = .508
    c = S/b
    rho = 1.225
    mu = 1.81e-5

    Cds = rand(10)
    Cls = rand(10)
    Cms = rand(10)
    alphas = rand(10)
    Res = rand(10)

    #create model
    polar_function = polar_constructor(Cds, Cls, Cms, alphas, Res)
    CRC3Parameters = BiWingTailSitter(m, I, X, s, S, b, rho, mu)
    CRC3Forces = biwing_tailsitter_forces_constructor(polar_function, CRC3Parameters)
    CRC3 = LowFidel(CRC3Parameters, CRC3Forces)

    F = [1,1]
    M = 1
    t = [0,1]
    uSpline = [FM.Akima(t, u[1]), FM.Akima(t, u[1])]
    dx = dynamics_2D!(zeros(6), x, (uSpline, CRC3), t[1])
    @test dx[1] == 1                            #Linear acceleration test (Vinf_dot)
    @test dx[2] == .1*180/pi                    #Flight angular velocity test (gamma_dot)
    @test dx[3] == 90.09009009009009*180/pi     #Body angular acceleration test (theta_dot_dot)
    

    F = [0,-2]
    M = 0
    x = [10,0,0,0,0,0]
    dt = .1
    t = 0:dt:31.4
    for i in 1:length(t)
        dx = dynamics(x,F,M,CRC3)
        x += dx*dt
        #println(x[6])
    end

    @test x[2] ≈ -360 atol=1           #Circular path test (gamma)

    F = [2,1]
    Fnew = deepcopy(F)
    M = 0
    x = [1,0,0,0,0,0]
    x0 = deepcopy(x)
    dt = .01
    t = 0:dt:30
    for i in 1:length(t)
        dx = dynamics(x,Fnew,M,CRC3)
        x += dx*dt
        vinf, gamma, thetadot, theta, posx, posy = x
        Fnew = [F[1]*cosd(gamma)+F[2]*sind(gamma), -F[1]*sind(gamma)+F[2]*cosd(gamma)]
        Fcheck = [Fnew[1]*cosd(gamma) - Fnew[2]sind(gamma), Fnew[1]*sind(gamma) + Fnew[2]*cosd(gamma)]
        # @show x
    end
    @show x[5] x[6]                                                         
    @test x[5] ≈ x0[1]*cosd(x0[2])*t[end] + .5*F[1]/m*t[end]^2 atol=1       #X position test (posx)
    @test x[6] ≈ x0[1]*sind(x0[2])*t[end] + .5*F[2]/m*t[end]^2 atol=1       #Y position test (posy)


    #forces test
    alphas = -16:10
    Res = [66000, 133000, 200000, 266000, 333000]
    Cds = [
    0.304701  0.301227  0.304131  0.302987  0.296036
    0.281761  0.275781  0.278783  0.277891  0.269864
    0.253494  0.248554  0.253074  0.250785  0.242528
    0.211665  0.222885  0.231003  0.225033  0.216555
    0.176881  0.193029  0.206483  0.200704  0.189419
    0.112822  0.155752  0.173539  0.171199  0.15488
    0.058415  0.087785  0.10945   0.122136  0.093951
    0.041396  0.042747  0.049159  0.031724  0.02201
    0.036257  0.028788  0.022803  0.019944  0.017844
    0.032727  0.024693  0.019051  0.016838  0.015141
    0.030291  0.021906  0.016717  0.01486   0.01347
    0.029063  0.019919  0.015366  0.013882  0.01273
    0.029331  0.018575  0.014366  0.013189  0.012308
    0.031229  0.019491  0.015074  0.013911  0.013054
    0.03375   0.021395  0.016957  0.015725  0.014813
    0.036992  0.023949  0.019622  0.018385  0.017461
    0.041376  0.027202  0.022948  0.021741  0.020835
    0.046949  0.031303  0.026969  0.025766  0.02487
    0.053645  0.036241  0.031837  0.030562  0.029616
    0.066135  0.041887  0.037457  0.036145  0.035157
    0.088572  0.04825   0.043756  0.042413  0.041393
    0.1054    0.055273  0.050705  0.049326  0.048276
    0.121486  0.062934  0.058286  0.056868  0.055784
    0.138769  0.071247  0.066484  0.065026  0.063902
    0.156057  0.080188  0.075288  0.073771  0.072602
    0.174349  0.089748  0.084694  0.083111  0.081885
    0.192409  0.099923  0.094685  0.093033  0.091745
    ]

    Cls = [
    -0.710933  -0.710933  -0.710933  -0.710933  -0.710933
    -0.641963  -0.641963  -0.641963  -0.641963  -0.641963
    -0.57245   -0.57245   -0.57245   -0.57245   -0.57245
    -0.502444  -0.502444  -0.502444  -0.502444  -0.502444
    -0.431996  -0.431996  -0.431996  -0.431996  -0.431996
    -0.361157  -0.361157  -0.361157  -0.361157  -0.361157
    -0.28998   -0.28998   -0.28998   -0.28998   -0.28998
    -0.218517  -0.218517  -0.218517  -0.218517  -0.218517
    -0.146823  -0.146823  -0.146823  -0.146823  -0.146823
    -0.074951  -0.074951  -0.074951  -0.074951  -0.074951
    -0.002957  -0.002957  -0.002957  -0.002957  -0.002957
    0.069105   0.069105   0.069105   0.069105   0.069105
    0.141179   0.141179   0.141179   0.141179   0.141179
    0.21321    0.21321    0.21321    0.21321    0.21321
    0.285143   0.285143   0.285143   0.285143   0.285143
    0.356923   0.356923   0.356923   0.356923   0.356923
    0.428495   0.428495   0.428495   0.428495   0.428495
    0.499803   0.499803   0.499803   0.499803   0.499803
    0.570795   0.570795   0.570795   0.570795   0.570795
    0.641415   0.641415   0.641415   0.641415   0.641415
    0.711613   0.711613   0.711613   0.711613   0.711613
    0.781334   0.781334   0.781334   0.781334   0.781334
    0.850528   0.850528   0.850528   0.850528   0.850528
    0.919145   0.919145   0.919145   0.919145   0.919145
    0.987134   0.987134   0.987134   0.987134   0.987134
    1.05445    1.05445    1.05445    1.05445    1.05445
    1.12104    1.12104    1.12104    1.12104    1.12104
    ]

    Cms = [
    -0.294631  -0.294352  -0.294585  -0.294493  -0.293935
    -0.275525  -0.27507   -0.275299  -0.275231  -0.27462
    -0.255848  -0.255493  -0.255818  -0.255653  -0.25506
    -0.235137  -0.235896  -0.236445  -0.236042  -0.235468
    -0.214901  -0.215925  -0.216777  -0.216411  -0.215696
    -0.192933  -0.19547   -0.196521  -0.196383  -0.195418
    -0.171772  -0.173382  -0.174569  -0.175265  -0.17372
    -0.152728  -0.152797  -0.15312   -0.15224   -0.15175
    -0.134169  -0.133824  -0.133548  -0.133416  -0.133319
    -0.115541  -0.115205  -0.114969  -0.114876  -0.114805
    -0.096834  -0.09652   -0.096325  -0.096256  -0.096204
    -0.078068  -0.077765  -0.077615  -0.077566  -0.077528
    -0.059267  -0.058958  -0.058838  -0.058804  -0.058778
    -0.040446  -0.04016   -0.040053  -0.040025  -0.040004
    -0.021594  -0.021347  -0.021259  -0.021234  -0.021216
    -0.002735  -0.002533  -0.002466  -0.002446  -0.002432
    0.016103   0.01626    0.016308   0.016321   0.016331
    0.034904   0.035009   0.035038   0.035046   0.035052
    0.053654   0.053694   0.053704   0.053707   0.053709
    0.072346   0.072295   0.072285   0.072282   0.07228
    0.091053   0.09079    0.090761   0.090752   0.090745
    0.109708   0.109159   0.109109   0.109094   0.109083
    0.128279   0.127381   0.12731    0.127288   0.127271
    0.146768   0.145434   0.14534    0.145312   0.145289
    0.165131   0.1633     0.163181   0.163145   0.163117
    0.18337    0.180956   0.180812   0.180767   0.180732
    0.20143    0.198386   0.198213   0.198159   0.198117
    ]

u = [1,1]
alpha = 0
Re = 13300
vinf = Re*mu/c/rho
x = [vinf,10,0,10,0,0]
D = 2*Cds[17,2]*S*.5*rho*x[1]^2 
L = 2*Cls[17,2]*S*.5*rho*x[1]^2
M1 = 2*Cms[17,2]*S*.5*rho*x[1]^2*c + L*X
F1 = [2*u[1]*cosd(alpha) + 2*u[2]*cosd(alpha) - D - m*9.81sind(x[2]), L - m*9.81cosd(x[2])] 
@show F1 M1 

polar_funct = polar_constructor(Cds, Cls, Cms, alphas, Res)
CRC3 = LowFidel(polar_funct, m, I, X, s, S, b, rho, mu)
F2,M2 = forces(x,u,CRC3)
@show F2 M2

@test F1 ≈ F2 atol= .01         #Forces test
@test M1 ≈ M2 atol= .00001      #Moment test
end