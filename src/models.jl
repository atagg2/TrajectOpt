abstract type Model end
abstract type Parameters end

struct ConventionalLowFidel{TF} <: Parameters
    m::TF
    I::TF
    X::TF
    L::TF  #tail moment arm
    SWing::TF
    bWing::TF
    STail::TF
    bTail::TF
    rho::TF
    mu::TF
    g::TF
end

struct ConventionalMidFidel{TF1, TF2, TF3, TF4} <: Parameters
    #mass and shape
    m::TF1
    I::TF1
    X::TF1
    L::TF1 #tail moment arm
    #wing 
    xleWing::TF2
    yleWing::TF2
    zleWing::TF2
    cWing::TF2
    twistWing::TF2
    camberWing::TF3
    #tail
    xleTail::TF2
    yleTail::TF2
    zleTail::TF2
    cTail::TF2
    twistTail::TF2
    camberTail::TF4
    #elevator
    cElevatorFraction::TF1
    bElevatorFraction::TF1
    #environment
    rho::TF1
    mu::TF1
    g::TF1
end

struct BiWingTailSitterLowFidel{TF} <: Parameters
    m::TF
    I::TF
    X::TF
    s::TF  #distance between wings
    S::TF
    b::TF
    rho::TF
    mu::TF
    g::TF
end

struct LowFidel{TF1, TF2} <: Model
    parameters::TF1
    forces::TF2  
end

struct HighFidel{TF1, TF2} <: Model
    parameters::TF1
    forces::TF2  
end

function polar_constructor(Cds,Cls,Cms,alphas,Res)
    function polar_function(alpha, Re)         #interpolation
        Cd = interp2d(linear, alphas, Res, Cds, alpha, Re)
        Cl = interp2d(linear, alphas, Res, Cls, alpha, Re)
        Cm = interp2d(linear, alphas, Res, Cms, alpha, Re)
        return [Cd;Cl;Cm]
    end
    return polar_function
end

function conventional_forces_constructor(wing_polar_function, tail_polar_function, parameters::ConventionalLowFidel)
    function forces_conventional_low_fidel(x, u)
        #current state and inputs
        vinf, gamma, thetadot, theta, posx, posy = x
        thrust, deflection = u
        #get physical and environmental parameters
        SWing = parameters.SWing
        bWing = parameters.bWing
        STail = parameters.STail
        bTail = parameters.bTail
        cWing = SWing/bWing
        cTail = STail/bWing
        L = parameters.L
        X = parameters.X
        m = parameters.m
        rho = parameters.rho
        mu = parameters.mu
        #current wing aerodynamics
        alphaWing = theta - gamma
        ReWing = rho*vinf*cWing/mu
        CdWing, ClWing, CmWing = wing_polar_function(alphaWing, ReWing)
        #current tail aerodynamics
        alphaTail = alphaWing + deflection
        ReTail = ReWing*cTail/cWing
        CdTail, ClTail, CmTail = tail_polar_function(alphaTail, ReTail)
        #forces on each surface
        f = [[CdWing, CdTail] [ClWing, ClTail] [CmWing, CmTail]]'
        f[1:2,1] *= 1/2*rho*vinf^2*SWing
        f[1:2,2] *= 1/2*rho*vinf^2*STail
        f[3,1] *= 1/2*rho*vinf^2*SWing*cWing
        f[3,2] *= 1/2*rho*vinf^2*STail*cTail
        #add thrust
        # f[1,:] .-= thrust*cosd(alpha)
        # f[2,:] .+= thrust*sind(alpha)
        #transform forces into total force and moment
        # f[1,:] .*= -1
        F = [thrust*cosd(alphaWing) - sum(f[1,:]) - m*9.81*sind(gamma), thrust*sind(alphaWing) + sum(f[2,:]) - m*9.81*cosd(gamma)]
        M = sum(f[3,:]) - X*(f[1,1]*sind(alphaWing) + f[2,1]*cosd(alphaWing)) - (X+L)*(f[1,2]*sind(alphaWing) + f[2,2]*cosd(alphaWing))
        return F, M
    end
    return forces_conventional_low_fidel
end

function biwing_tailsitter_forces_constructor(polar_function, parameters::BiWingTailSitterLowFidel)
    function forces_biwing_tailsitter(x, u)
        vinf, gamma, thetadot, theta, posx, posy = x
        #get physical and environmental parameters
        S = parameters.S  
        b = parameters.b 
        c = S/b
        s = parameters.s
        X = parameters.X
        m = parameters.m
        rho = parameters.rho
        mu = parameters.mu
        #current top state
        alpha = theta - gamma
        v = vinf - s/2*thetadot*cosd(alpha)
        Re = rho*v*c/mu
        #calculate aerodynamic forces
        #top aerodynamics at current time step 
        CdTop, ClTop, CmTop = polar_function(alpha,Re)
        #current bottom state
        v = vinf + s/2*thetadot*cosd(alpha)
        Re = rho*v*c/mu
        #bottom aerodynamics current time step 
        CdBottom, ClBottom, CmBottom = polar_function(alpha,Re)
        f = [[CdTop, CdBottom] [ClTop, ClBottom] [CmTop, CmBottom]]'
        #denormalization
        f[1:2,:] *= 1/2*rho*v^2*S               # change v here!!
        f[3,:] *= 1/2*rho*v^2*S*c
        #add thrust 
        # f[1,:] .-= 2*u*cosd(alpha)
        # f[2,:] .+= 2*u*sind(alpha)
        #forces in the x and y directions with respect to the inertial frame
        # f[1,:] .*= -1
        F = [sum(u)*cosd(alpha) - sum(f[1,:]) - m*9.81*sind(gamma), sum(u)*sind(alpha) + sum(f[2,:]) - m*9.81*cosd(gamma)]
        #total moment due to forces and torques
        M = sum(f[3,:]) + s/2*(u[2] - u[1] + (f[1,1] - f[1,2])*cosd(alpha) + (f[2,2] - f[2,1])*sind(alpha)) - 
            X*((f[2,1] + f[2,2])*cosd(alpha) + (f[1,1] + f[1,2])*sind(alpha))
        return F, M
    end
    return forces_biwing_tailsitter
end

function conventional_forces_constructor(parameters::ConventionalMidFidel)
    # general parameters
    L = parameters.L
    X = parameters.X
    m = parameters.m
    rho = parameters.rho
    b = parameters.bElevatorFraction
    c = parameters.cElevatorFraction

    # wing discretization parameters
    ns = 12
    nc = 6
    spacing_s = Uniform()
    spacing_c = Uniform()
    # wing geometry (right half of the wing)
    xle = parameters.xleWing
    yle = parameters.yleWing
    zle = parameters.zleWing
    chord = parameters.cWing
    twist = parameters.twistWing
    phi = [0.0, 0.0]    # horizontal wing
    fc = parameters.camberWing #fill((xc) -> 0, 2) # camberline function for each section - see what to do about this
    # reference parameters
    Vinfref = 10
    Sref = 8#(chord[1] + chord[2])*(yle[2] - yle[1])/2
    cref = .817#chord[1]
    bref = 10#2*yle[2]
    rref = [-.5, 0.0, 0.0]
    ref = Reference(Sref, cref, bref, rref, Vinfref)
    # construct surface
    gridWing, wing = wing_to_surface_panels(xle, yle, zle, chord, twist, phi, ns, nc;
        fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)

    # tail discretization parameters
    ns = 6
    nc = 3
    spacing_s = Uniform()
    spacing_c = Uniform()
    # tail geometry (right half of the wing)
    xle = [parameters.xleTail[1], parameters.xleTail[2]*b, parameters.xleTail[2]*b, parameters.xleTail[2]]
    yle = [parameters.yleTail[1], parameters.yleTail[2]*b, parameters.yleTail[2]*b, parameters.yleTail[2]]
    zle = [parameters.zleTail[1], parameters.zleTail[2]*b, parameters.zleTail[2]*b, parameters.zleTail[2]]
    chord = [parameters.cTail[1]*(1-c), parameters.cTail[1]*(1-c) - xle[2], parameters.cTail[1] - (parameters.cTail[1] - parameters.cTail[2])*b, parameters.cTail[2]]        
    twist = [parameters.twistTail[1], parameters.twistTail[2]*b, parameters.twistTail[2]*b, parameters.twistTail[2]]
    phi = [0.0, 0.0, 0.0, 0.0]
    fc = [parameters.camberTail[1],parameters.camberTail[1],parameters.camberTail[2], parameters.camberTail[2]] 
    # construct surface
    gridTail, tail = wing_to_surface_panels(xle, yle, zle, chord, twist, phi, ns, nc;
        fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)
    VL.translate!(gridTail, [L, 0.0, 0.0])
    VL.translate!(tail, [L, 0.0, 0.0])

    # elevator geometry (right half of the wing)
    xle = [0.0, 0.0]
    yle = parameters.yleTail*b
    zle = parameters.zleTail
    chord = [parameters.cTail[1]*c, chord[3] - chord[2]]
    twist = parameters.twistTail
    phi = [0.0, 0.0]
    fc = parameters.camberTail #
    # construct surface
    gridElevator, elevator = wing_to_surface_panels(xle, yle, zle, chord, twist, phi, ns, nc;
        fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)
    LElevator = L + parameters.cTail[1] - chord[1]
    VL.translate!(gridElevator, [LElevator, 0.0, 0.0])
    VL.translate!(elevator, [LElevator, 0.0, 0.0])

    # create vector containing all surfaces
    surfaces = [wing]#, tail, elevator]
    # we can use symmetry since the geometry and flow conditions are symmetric about the X-Z axis
    symmetric = fill(true, length(surfaces))
    # freestream parameters
    alpha = 0
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinfref, alpha, beta, Omega)
    # perform steady state analysis
    aircraft = steady_analysis(surfaces, ref, fs; symmetric=symmetric)
        
    function forces_conventional_high_fidel(x, u)
        # get current state and input
        Vinf, gamma, thetadot, theta, posx, posy = x
        alpha = theta - gamma
        thrust, deflection = u

        # reset the freestream
        fs = Freestream(Vinf, alpha*pi/180, beta, Omega)

        # reset the elevator geometry
        R = [
            cosd(deflection) 0 sind(deflection)
            0 1 0
            -sind(deflection) 0 cosd(deflection)
        ]
        VL.translate!(gridElevator, [-LElevator, 0.0, 0.0])
        for i in 1:length(gridElevator[1,:,1])
            for j in 1:length(gridElevator[1,1,:])
                gridElevator[:,i,j] = R*gridElevator[:,i,j]
            end
        end
        VL.translate!(gridElevator, [LElevator, 0.0, 0.0])
        VL.update_surface_panels!(elevator, gridElevator)

        # update aircraft
        surfaces = [wing]#, tail, elevator]
        steady_analysis!(aircraft, surfaces, ref, fs; symmetric=symmetric)

        # retrieve near-field forces
        CF, CM = body_forces(aircraft; frame=Wind())

        # extract 2D force and moment coefficients
        Cd, _, Cl = CF
        _, Cm, _ = CM
        @show Cm

        # denormalize
        F = [Cd,Cl]*.5*rho*Vinfref^2*Sref
        M = Cm*.5*rho*Vinfref^2*Sref*cref

        # add thrust and gravity
        F[1] *= -1
        F[1] += thrust*cosd(alpha) - m*9.81*sind(gamma)
        F[2] += thrust*sind(alpha) - m*9.81*cosd(gamma)

        # reset elevator to 0
        VL.translate!(gridElevator, [-LElevator, 0.0, 0.0])
        for i in 1:length(gridElevator[1,:,1])
            for j in 1:length(gridElevator[1,1,:])
                gridElevator[:,i,j] = R'*gridElevator[:,i,j]
            end
        end
        VL.translate!(gridElevator, [LElevator, 0.0, 0.0])
        VL.update_surface_panels!(elevator, gridElevator)

        # f = [[CdAircraft, CdElevator] [ClAircraft, ClElevator] [CmAircraft, CmElevator]]'
        # f[1:2,:] *= 1/2*rho*Vinf^2*Sref
        # f[3,:] *= 1/2*rho*Vinf^2*Sref*cref
        
        # F = [thrust*cosd(alphaAircraft) - sum(f[1,:]) - m*9.81*sind(gamma), thrust*sind(alphaAircraft) + sum(f[2,:]) - m*9.81*cosd(gamma)]
        # M = sum(f[3,:]) - X*(f[1,1]*sind(alphaAircraft) + f[2,1]*cosd(alphaAircraft)) - (X+LElevator)*(f[1,2]*sind(alphaAircraft) + f[2,2]*cosd(alphaAircraft))
        #double check moment

        return F, M
    end
    return forces_conventional_high_fidel
end