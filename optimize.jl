# using LinearAlgebra
# using Plots
# using Interpolations
# using Random
# include("track.jl")

# This is Reader's PSO for project2

## Calling:
function optimize(f, g, c, x0, n, prob)
f_p = quadratic_penalty_function2(f,c)
N = 12
v_range = (-3,-1)
population = initialize_population(x0, N, v_range)
xhistory = particle_swarm_optimization(f_p, population, n; w=0.7, c1=1.2, c2=1.2)
x_best = xhistory[end-length(x0)+1 : end]
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