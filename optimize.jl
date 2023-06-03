# @counted function simple2(x::Vector)
#     return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
# end


# @counted function simple2_constraints(x::Vector)
#     return [(x[1]-1)^3 - x[2] + 1,
#             x[1] + x[2] - 2]
# end

# function simple2_init()
#     return rand(2) .* 2.0 .- 1.0
# end

# For this file we need to output f, c, x0, and maybe n

using LinearAlgebra
using Plots
using Interpolations
using Random
include("track.jl")
include("problem.jl")

# Initial design point (x0)
track_bound = build_track()
U = create_U0(track_bound)
S = compute_state(U, track_bound)

# Initial design point (x0)
x0 = U

# Initialize n
n = 5000

# Constraint function (c)
function con(U)
    S = compute_state(U, track_bound)
    constraint_vector = compute_track_car_constraints(S, U, track_bound; ϵ=1e-5)
    return constraint_vector
end

# Objective function (f)
function fun(U)
    design_point = compute_state(U, track_bound)
    total_time = compute_total_time(design_point)
    # display(total_time)
    return total_time
end

# # Optimize
# function optimize(func, cons, car_limits, x0, n_iter, track)
#     shistory, uhistory, fhistory, chistory = quad_penalty_particle_swarm(func, cons, car_limits, x0, n_iter, track)
#     return shistory, uhistory, fhistory, chistory
# end

# # Method #1 - Quadratic penalty + Particle swarm optimization
# function quad_penalty_particle_swarm(func, cons, car_limits, x0, n_iter, track)
#     return shistory, uhistory, fhistory, chistory
# end



# This is based on Reader's PSO for project2

## Calling:
function optimize(f, c, x0, n)
        f_p = quadratic_penalty_function(f,c)
        N = 200
        v_range = (0,100)
        population = initialize_population(x0, N, v_range)
        xhistory = particle_swarm_optimization(f_p, population, n; w=0.7, c1=1.2, c2=1.2)
        # display(xhistory)
        x_best = xhistory[(end-length(x0)+1) : end]
    return x_best
end

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

# function initialize_population(X0, N, v_range)
#     population = Particle[]
#     for i in 1:N
#         x = copy(X0)
#         # v = rand(length(X0)) .* (v_range[2] - v_range[1]) .+ v_range[1]
#         v = rand(2,lastindex(X0[1,:])) .* (v_range[2] - v_range[1]) .+ v_range[1]
#         x_best = copy(X0)
#         push!(population, Particle(x, v, x_best))
#     end
#     return population
# end

function initialize_population(X0, N, v_range)
    population = Particle[]
    for i in 1:N
        x = copy(X0)
        v = rand(size(U,1), size(U,2)) .* (v_range[2] - v_range[1]) .+ v_range[1]
        # display(v)
        x_best = copy(X0)
        push!(population, Particle(x, v, x_best))
    end
    return population
end

function particle_swarm_optimization(f, population, k_max; w=1, c1=1, c2=1)
    # n = length(population[1].x)
    n1 = size((population[1].x),1)
    n2 = size(population[1].x,2)
   
    x_best, y_best = copy(population[1].x_best), Inf
    xhistory = [copy(x_best) ]
    for P in population
        y = f(P.x)
        if y < y_best; x_best[:], y_best = P.x, y; end 
    end
    for k in 1 : k_max/((6*length(population)+2))
        for P in population
            r1, r2 = rand(n1, n2), rand(n1, n2)
            # display(P.x)
            # display(P.v)
            # display(P.x_best)
            # display(n1)
            P.x .+= P.v
            P.v = w.*P.v + c1.*r1.*(P.x_best - P.x) + c2.*r2.*(x_best - P.x)
            y = f(P.x)
            if y < y_best; x_best[:], y_best = P.x, y; end
            if y < f(P.x_best); P.x_best[:] = P.x; end
        end
        push!(xhistory, x_best)
        # append!(xhistory, x_best)
    end 
    return xhistory 
end


optimize(fun, con, x0, n)