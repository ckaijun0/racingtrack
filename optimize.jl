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
        S = [semi_holonomic_model(S[n],U[n]) for n in 1:(N-1)]
    end
    return S
end

# Create the semi holonomic model, to calculate S+1 from S and U
function semi_holonomic_model(S,U)
    S_next = copy(S)
    S_next[1] = S[1] + (S[3]*S[5]+0.5*U[2]*S[5].^2).*cos(U[1])
    S_next[2] = S[2] + (S[3]*S[5]+0.5*U[2]*S[5].^2).*sin(U[1])
    S_next[3] = S[3] + U[2].*S[5]
    S_next[4] = S[4] + U(1)
    d = sqrt((S_next[1]-S[1]).^2 + (S_next[2]-S[2]).^2)
    S_next[5] = max((-S(3)+sqrt(S[3].^2-2*U[2]*d))/U[2],(-S(3)-sqrt(S[3].^2-2*U[2]*d))/U[2])
    return S_next
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

# This is Reader's PSO for project2

## Calling:
# f_p = quadratic_penalty_function2(f,c)
# N = 12
# v_range = (-3,-1)
# population = initialize_population(x0, N, v_range)
# xhistory = particle_swarm_optimization(f_p, population, n; w=0.7, c1=1.2, c2=1.2)

function quadratic_penalty_function(f, c; ρ=9999999999)
    function h(x)
        penalty = sum(((c(x).>0) .* c(x)).^2)
        return f(x) + (ρ.^8) * penalty
    end
    return h
end

mutable struct Particle 
    x
    v
    x_best
end

function initialize_population(X0, N, v_range)
    population = Particle[]
    for i in 1:N
        x = copy(X0)
        v = rand(length(X0)) .* (v_range[2] - v_range[1]) .+ v_range[1]
        x_best = copy(X0)
        push!(population, Particle(x, v, x_best))
    end
    return population
end



function particle_swarm_optimization(f, population, k_max; w=1, c1=1, c2=1)
    n = length(population[1].x)
   
    x_best, y_best = copy(population[1].x_best), Inf
    xhistory = copy(x_best) 
    for P in population
        y = f(P.x)
        if y < y_best; x_best[:], y_best = P.x, y; end 
    end
    for k in 1 : k_max/((6*length(population)+2))
        for P in population
            r1, r2 = rand(n), rand(n)
            P.x += P.v
            P.v = w*P.v + c1*r1.*(P.x_best - P.x) + c2*r2.*(x_best - P.x)
            y = f(P.x)
            if y < y_best; x_best[:], y_best = P.x, y; end 
            if y < f(P.x_best); P.x_best[:] = P.x; end
        end
        append!(xhistory, x_best)
    end 
    return xhistory 
end