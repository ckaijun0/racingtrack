using LinearAlgebra
using Plots
using Interpolations
using Random
include("track.jl")

function optimize(X,U)
    return X,U,T
end

# Parameters for user input
vmax = 50 # m/s (111.8mph)
amax = 3 # m/s^2
θmax = 45 # degrees 
ωmax = 15 # ω=θ/t (assuming full lock in 3seconds)
m = 1500 # kg
friction_coeff = 0.7 # 0.7=dry road, 0.4=wet road
N = 50 # number of points on track
ϵ = 1e-5 # permitted error on position due to floating point errors

# Derived parameters
Fc = friction_coefficient*m*9.81 # friction = f(normal_force)
Ac = Fc/m # centripetal acceleration

# State variables vector
# X1 = x-coordinate
X = [ones(5) for _ in 1:N] 
U = [ones(2) for _ in 1:N] # Input variables vector

mutable struct create_track 
    bounds
    start_point
end

track = create_track(create_track_width(centerline),start_point)

# Optimize
function optimize(func, cons, car_limits, x0, n_iter, track)
    shistory, uhistory, fhistory, chistory = quad_penalty_particle_swarm(func, cons, car_limits, x0, n_iter, track)

end

# Method #1 - Quadratic penalty + Particle swarm optimization
function quad_penalty_particle_swarm(func, cons, car_limits, x0, n_iter, track)
    return shistory, uhistory, fhistory, chistory
end

# Objective
function objective(X)
    return sum()
end

# Constraints
function constraints(s, u, track_bound)
    coord = [s[1]; s[2]]
    bound1, bound2 = track_bound[1], track_bound[2]
    v, θ = s[3], s[4]
    ω, a = u[1], u[2]
    g = [1/ϵ*(dist(coord,bound1))-1;
        v/vmax-1;
        a/amax-1;
        θ/θmax-1;
        ω/ωmax-1;
        v*ω/Fc/m-1]
    return g
end

