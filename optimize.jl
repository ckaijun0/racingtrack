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
U_initial = create_U0(track_bound)
# U[1,:] = [0.02 for _ in 1:lastindex(U[1,:])] # Artificial deviation from the optimal line
# U_initial = copy(U)
S_initial = compute_state(U_initial, track_bound)
total_time_initial = compute_total_time(S_initial)

# Initial design point (x0)
x0 = U_initial

# Initialize n
n = 50

# # Constraint function (c)
# function con(U1)
#     S = compute_state(U1, track_bound)
#     constraint_vector = compute_track_car_constraints(S, U1, track_bound; ϵ=1e-5)
#     return constraint_vector
# end

# # Objective function (f)
# function fun(U2)
#     design_point = compute_state(U2, track_bound)
#     total_time = compute_total_time(design_point)
#     #display(design_point)
#     #display(total_time)
#     return abs(total_time)
# end

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
        f_p = interior_point_method(f,c)
        N = 200
        x_range = (-10,10)
        v_range = (-4,4)
        population = initialize_population(x0, N, x_range, v_range)
        xhistory, fhistory = particle_swarm_optimization(f_p, population, n; w=0.1, c1=0.7, c2=0.7)
        x_best = xhistory[end]
        
        # f_p = mixed_penalty_function(f,c)
        # xhistory, fhistory, chistory = hooke_jeeves(f_p,x0,c,n,0.1,0.5)

        # α=0.1
        # β=0.8
        # xhistory  = nesterov_momentum!(α, β, v, f_p, c, x0)
        # push!(xhistory, xnext)
        # push!(fhistory, f(xnext))

    return xhistory, fhistory #, chistory
end

# function nesterov_momentum!(α, β, v, f, g, x)
#     xhistory = [copy(x)]
#     fhistory = [f(x)]

#     for i = 1:n/2-1
#         g = g(x_next+β*v)
#         gnorm = g/sqrt(sum(dot(g,g)))
#         v[:] = β*v - α*gnorm
#         x_next = x + v
#     end
#     return xhistory, fhistory
# end

function quadratic_penalty_function(f, c; ρ=999999999999999)
    function h(x)
        penalty = sum(((c(x).>0) .* c(x)).^2)
        return f(x) + (ρ.^9999999999999999) * penalty
    end
    return h
end

function mixed_penalty_function(f, c; ρ=999999999)
    function h(x)
        penalty = sum(((c(x).>0) .* c(x)).^2)
        return f(x) + (ρ.^99999999) * penalty + 999999 * sum(c(x).>0 .* c(x))
    end
    return h
end

function interior_point_method(f, c; ρ=99999999999)
    function h(x)
        # penalty = -sum((c(x).>-1) .* log.(-c(x)))
        penalty = -sum(1/c(x))
        return f(x) + 1/ρ*penalty
    end
    return h
end
mutable struct Particle 
    x
    v
    x_best
end

function initialize_population(X0, N, x_range, v_range)
    population = Particle[]
    for i in 1:N
        x = copy(X0)
        x = x + rand(size(X0,1), size(X0,2)) .* (x_range[2] - x_range[1]) .+ x_range[1]
        v = rand(size(X0,1), size(X0,2)) .* (v_range[2] - v_range[1]) .+ v_range[1]
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
    xhistory = [copy(x_best)]
    fhistory = [fun(x_best)]
    for P in population
        y = fun(P.x)
        if y < y_best; 
            x_best[:], y_best = P.x, y; 
        end 
    end
    # for k in 1 : k_max/((6*length(population)+2))
    for k in 1 : k_max
        for P in population
            r1, r2 = rand(n1, n2), rand(n1, n2)
            # display(P.x)
            # display(P.v)
            # display(P.x_best)
            # display(n1)
            P.x += P.v
            P.v = w.*P.v + c1.*r1.*(P.x_best .- P.x) + c2.*r2.*(x_best .- P.x)
            y = fun(P.x)
            # display(fun(P.x))
            if y < y_best; 
                x_best[:], y_best = P.x, y; 
            end
            # display(y)
            if y < fun(P.x_best);
                # display(x_best) 
                # P.x_best[:] .= P.x;
                P.x_best = copy(P.x) 
            end
        end
        push!(xhistory, copy(x_best))
        push!(fhistory, copy(fun(x_best)))
        # append!(xhistory, x_best)
    end 
    return xhistory, fhistory
end

function hooke_jeeves(f, x, c, counts, α, γ)
    y = f(x)
    xhistory = [copy(x)]
    fhistory = [copy(fun(x))]
    c_vec = c(x)
    chistory = [copy(sum((c_vec.>0).*c_vec))]
    for k in 1:counts
        improved = false
        x_best, y_best = x, y

        for i in 1:length(x)
            for sgn in (-1, 1)
                x′ = copy(x)
                x′[i] += sgn * α
                y′ = f(x′)

                if y′ < y_best
                    x_best, y_best, improved = copy(x′), y′, true
                end
            end
        end

        x, y = copy(x_best), y_best

        if !improved
            α *= γ
        end
        # display(x_best)
        push!(xhistory, copy(x_best))
        push!(fhistory, copy(fun(x_best)))
        c_vec = copy(c(x_best))
        push!(chistory, sum((c_vec.>0).*c_vec))
    end

    return xhistory, fhistory, chistory
    
end

xhistory, fhistory = optimize(fun, con, x0, n)
U_optimal = xhistory[end]
total_time_optimal = fhistory[end]
# U_optimal = optimize(fun, con, x0, n)
S_optimal = compute_state(U_optimal, track_bound)
# time_taken_optimal = compute_total_time(S_optimal)
println("~~Time~")
display(total_time_optimal)
println("~~State~~")
display(S_optimal)
println("~~Inputs~~")
display(U_optimal)
plot_track(S_optimal, U_optimal, track_bound)
