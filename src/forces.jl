struct Inputs{T1, T2}
    rotorCommands::Vector{T1}
    controlSurfaceCommands::Vector{T2}
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

function forces_low_fidel(parameters::Parameters)
    function forces(x, u::Inputs)
        vinf, gamma, thetadot, theta, posx, posy = x
        alpha = theta - gamma
        thrusts = u.rotorCommands
        deflections = u.controlSurfaceCommands
        inertia = parameters.inertia
        environment = parameters.environment
        X = inertia.cog[1]
        m = inertia.m
        rho = environment.rho
        mu = environment.mu
        g = environment.g

        F = zeros(0)
        M = 0
        for i in 1:length(parameters.surfaces)
            s = parameters.surfaces[i]
            alphaNew = alpha + deflections[i] 
            c = s.S/s.b
            Re = rho*vinf*c/mu
            Cd, Cl, Cm = s.polar(alphaNew, Re)
            D = Cd*.5*rho*vinf^2*s.S
            L = Cl*.5*rho*vinf^2*s.S
            F[1] -= D
            F[2] += L
            
            M += Cm*.5*rho*vinf^2*s.S*c - r[1]*(D*sind(alpha) + 
                L*cosd(alpha)) + r[3]*(D*cosd(alpha) - L*sind(alpha)) 
        end
        for i in 1:length(parameters.rotors)
           F[1] += thrusts[i]*cosd(alpha)
           F[2] += thrusts[i]*sind(alpha)
           M += thrusts[i]*parameters.rotors[i].r[3]
        end

        return F, M
    end
    return forces
end

function forces_conventional_low_fidel(parameters::Parameters)
    function forces(x, u)
        #current state and inputs
        vinf, gamma, thetadot, theta, posx, posy = x
        thrust, deflection = u
        #get physical and environmental parameters
        inertia = parameters.inertia
        environment = parameters.environment
        X = inertia.cog[1]
        m = inertia.m
        rho = environment.rho
        mu = environment.mu
        g = environment.g
        #get wing
        wing = parameters.surfaces[1]
        wing_polar = wing.polar
        SWing = wing.S
        bWing = wing.b
        cWing = SWing/bWing
        #get tail
        tail = parameters.surfaces[2]
        L = tail.r[1]
        tail_polar = tail.polar
        STail = tail.S
        bTail = tail.b
        cTail = STail/bTail
        #current wing aerodynamics
        alphaWing = theta - gamma
        ReWing = rho*vinf*cWing/mu
        CdWing, ClWing, CmWing = wing_polar(alphaWing, ReWing)
        #current tail aerodynamics
        alphaTail = alphaWing + deflection
        ReTail = ReWing*cTail/cWing
        CdTail, ClTail, CmTail = tail_polar(alphaTail, ReTail)
        #forces on each surface
        f = [[CdWing, CdTail] [ClWing, ClTail] [CmWing, CmTail]]'
        f[1:2,1] *= 1/2*rho*vinf^2*SWing
        f[1:2,2] *= 1/2*rho*vinf^2*STail
        f[3,1] *= 1/2*rho*vinf^2*SWing*cWing
        f[3,2] *= 1/2*rho*vinf^2*STail*cTail
        #add thrust
        F = [thrust*cosd(alphaWing) - sum(f[1,:]) - m*g*sind(gamma), 
             thrust*sind(alphaWing) + sum(f[2,:]) - m*g*cosd(gamma)]
        M = sum(f[3,:]) - X*(f[1,1]*sind(alphaWing) + 
            f[2,1]*cosd(alphaWing)) - 
            (X+L)*(f[1,2]*sind(alphaWing) + 
            f[2,2]*cosd(alphaWing))
        return F, M
    end
    return forces
end

function forces_conventional_mid_fidel(parameters::Parameters)
    # general parameters
    L = parameters.surfaces[2].xle[1]
    cog = parameters.inertia.cog
    X = cog[1]
    m = parameters.inertia.m
    rho = parameters.environment.rho
    # wing discretization parameters
    ns = 12
    nc = 6
    spacing_s = Uniform()
    spacing_c = Uniform()
    # surfaces
    surfaces = Array{Array{SurfacePanel{Float64}}}(undef, 3,1)
    grids = Array{Array{Float64}}(undef, 3,1)
    for i in 1:length(parameters.surfaces)
        s = parameters.surfaces[i]
        phi = [0.0, 0.0]    # horizontal wing
        fc = s.camber
        grids[i], surfaces[i] = wing_to_surface_panels(s.xle, s.yle, s.zle, s.chord, s.twist, phi, ns, nc;
            fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)
    end
    surfaces = [surfaces[1],surfaces[2],surfaces[3]]      #there is probably a better way to do this
    # we can use symmetry since the geometry and flow conditions are symmetric about the X-Z axis
    symmetric = fill(true, length(surfaces))
    # reference parameters
    Vinfref = 10
    Sref = 1#(chord[1] + chord[2])*(yle[2] - yle[1])/2
    cref = 1#chord[1]
    bref = 1#2*yle[2]
    rref = cog
    ref = Reference(Sref, cref, bref, rref, Vinfref)
    # freestream parameters
    alpha = 0
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinfref, alpha, beta, Omega)
    # perform steady state analysis
    aircraft = steady_analysis(surfaces, ref, fs; symmetric=symmetric)
        
    function forces(x, u)
        # get current state and input
        Vinf, gamma, thetadot, theta, posx, posy = x
        alpha = theta - gamma
        omega, deflection = u
        thrust, _ = omega_2_thrust(omega, Vinf, parameters.environment, parameters.rotors[1])
        # reset the freestream
        fs = Freestream(Vinf, alpha*pi/180, beta, Omega)
        # reset the elevator geometry
        R = [
            cosd(deflection) 0 sind(deflection)
            0 1 0
            -sind(deflection) 0 cosd(deflection)
        ]
        gridElevator = grids[3]
        VL.translate!(grids[3], [-L, 0.0, 0.0])
        @show R[1][1]
        for i in 1:length(gridElevator[1,:,1])
            for j in 1:length(gridElevator[1,1,:])
                gridElevator[:,i,j] = R*grids[3][:,i,j]
            end
        end
        VL.translate!(gridElevator, [L, 0.0, 0.0])
        VL.update_surface_panels!(surfaces[3], gridElevator)
        # update aircraft
        steady_analysis!(aircraft, surfaces, ref, fs; symmetric=symmetric)
        # retrieve near-field forces
        CF, CM = body_forces(aircraft; frame=Wind())
        # extract 2D force and moment coefficients  
        Cd, _, Cl = CF
        _, Cm, _ = CM
        # denormalize
        F = [Cd,Cl]*.5*rho*Vinf^2*Sref
        M = Cm*.5*rho*Vinf^2*Sref*cref
        # add thrust and gravity
        F[1] *= -1
        # @show F thrust alpha gamma
        F[1] += thrust*cosd(alpha) - m*9.81*sind(gamma)
        # @show F
        F[2] += thrust*sind(alpha) - m*9.81*cosd(gamma)
        # reset elevator to 0
        VL.translate!(gridElevator, [-L, 0.0, 0.0])
        for i in 1:length(gridElevator[1,:,1])
            for j in 1:length(gridElevator[1,1,:])
                gridElevator[:,i,j] = R'*gridElevator[:,i,j]
            end
        end
        VL.translate!(gridElevator, [L, 0.0, 0.0])
        VL.update_surface_panels!(surfaces[3], gridElevator)

        return F, M
    end
    return forces
end

function vlm_surface_setup(parameters::Parameters)
    return aircraft
end

function omega_2_thrust(omega, Vinf, environment::Environment, rotor::CCBladeRotor)
    CCRotor = CC.Rotor(rotor.Rhub, rotor.Rtip, rotor.n)     #this is confusing because there are two types of rotors
    sections = Section.(rotor.r, rotor.chord, rotor.theta, Ref(rotor.af))
    op = simple_op.(Vinf, omega, rotor.r, environment.rho)
    BEMSolution = CC.solve.(Ref(CCRotor), sections, op)
    thrust, torque = thrusttorque(CCRotor, sections, BEMSolution)
    return thrust, torque
end