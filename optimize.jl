using LinearAlgebra
using Plots
using Interpolations
using Random
include("track.jl")

function optimize(X,U)
    return X,U,T
end

# Specify the centerline track coordinates here (don't duplicate start and end points)
x_points = [1; 2; 3; 3; 2; 1; 0; 0]
y_points = [0; 0; 1; 2; 3; 3; 2; 1]
center_line = [x_points, y_points]
width = 0.1

# Parameters for user input
vmax = 50 # m/s (111.8mph)
amax = 3 # m/s^2
θmax = 45 # degrees 
ωmax = 15 # ω=θ/t (assuming full lock in 3seconds)
m = 1500 # kg
friction_coefficient = 0.7 # 0.7=dry road, 0.4=wet road
N = lastindex(x_points)+1 # number of points on track + 1 (repeats start point)
ϵ = 1e-5 # permitted error on position due to floating point errors

# Derived parameters
Fc = friction_coefficient*m*9.81 # friction = f(normal_force)
Ac = Fc/m # centripetal acceleration

# State variables vector
S = ones(5, N)
U = ones(2, N)

track_bound = create_track_width(center_line, width)
# S[:,1] = track_bound[3][1]

# Optimize
function optimize(func, cons, car_limits, x0, n_iter, track)
    shistory, uhistory, fhistory, chistory = quad_penalty_particle_swarm(func, cons, car_limits, x0, n_iter, track)
    return shistory, uhistory, fhistory, chistory
end

# Method #1 - Quadratic penalty + Particle swarm optimization
function quad_penalty_particle_swarm(func, cons, car_limits, x0, n_iter, track)
    return shistory, uhistory, fhistory, chistory
end

# Compute State vector, S, given S[1] and Input vector, U
function compute_state(initial_point, S, U)
    S[:,1] = [initial_point; 0; 0; 0] # Give initial point, 
    N = lastindex(S[1,:])
    for n = 1:N
        S = [semi_holonomic_model(S[n],) for n in 1:(N-1)]
    end
    return S
end

# Objective: Total time
function total_time(S)
    return sum(S[5,:])
end

# Constraints
function constraints(S, U, track_bound; ϵ=1e-5)
    coord = S[1:2,:]
    bound1, bound2 = reduce(hcat,track_bound[1])', reduce(hcat,track_bound[2])' # Convert vector of vector to array
    v, θ = S[3,:], S[4,:]
    ω, a = U[1,:], U[2,:]
    # Computes distance from each point to checkpoint lines
    N = lastindex(S[1, :])
    all_distances = [abs(norm(coord[:,n]-bound1[:,n])+norm(coord[:,n]-bound2[:,n])-norm(bound1[:,n]-bound2[:,n])) for n in 1:N]
    # Creates a vector with raw constraint penalties: size (N+1)*5 rows x 1 column
    constrain_vector = [all_distances.-ϵ;
                        v./vmax.-1;
                        a./amax.-1;
                        θ./θmax.-1;
                        ω./ωmax.-1;
                        v.*ω./(Fc*m).-1]
    return constrain_vector
end

c = constraints(S, U, track_bound; ϵ=1e-5)
# print(c)

