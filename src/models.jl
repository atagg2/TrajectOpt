abstract type Surface end
abstract type Rotor end

struct Environment{T}
    rho::T
    mu::T
    g::T
end

struct Inertia{T}
    m::T
    I::T
    cog::Vector{T}
end

struct SimpleSurface{T1, T2} <: Surface
    S::T1
    b::T1
    r::Vector{T1}
    polar::T2
end

struct SimpleRotor{T} <: Rotor
    r::Vector{T}
end

struct VLMSurface{T1, T2} <: Surface
    xle::Vector{T1}
    yle::Vector{T1}
    zle::Vector{T1}
    chord::Vector{T1}
    twist::Vector{T1}
    camber::T2
end

struct CCBladeRotor{T1, T2, T3} <: Rotor
    pos::Vector{T1}
    Rtip::T1
    Rhub::T1
    n::T2
    r::Vector{T1} 
    chord::Vector{T1}
    theta::Vector{T1}
    af::T3
end

struct Parameters
    environment::Environment
    inertia::Inertia
    surfaces::Vector{Surface}
    rotors::Vector{Rotor}
end

struct Model{T}
    parameters::Parameters
    forces::T
end
