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

# function optimize(X,U)
#     return X,U,T
# end

#######################################################################################################################################

# Specify the centerline track coordinates here (don't duplicate start and end points)
x_points = [1; 2; 3; 3; 2; 1; 0; 0]
y_points = [0; 0; 1; 2; 3; 3; 2; 1]
center_line = [x_points, y_points]
width = 0.1

# Parameters for user input
vmax = 50 # m/s (111.8mph)
amax = 3 # m/s^2
θmax = 90 # degrees 
ωmax = 15 # ω=θ/t (assuming full lock in 3seconds)
m = 1500 # kg
friction_coefficient = 0.7 # 0.7=dry road, 0.4=wet road
N = lastindex(x_points)+1 # number of points on track + 1 (repeats start point)
ϵ = 1e-5 # permitted error on position due to floating point errors

# Derived parameters
Fc = friction_coefficient*m*9.81 # friction = f(normal_force)
Ac = Fc/m # centripetal acceleration

# State variables vector
# S = ones(5, N)
U = zeros(2, N)
U[1,:] = [pi/4 for _ = 1:9]

track_bound = create_track_width(center_line, width)

# Create the semi holonomic model, to calculate S+1 from S and U
function semi_holonomic_model(S, U, local_track_bound)
    S_next = copy(S)
    point_R, point_L = local_track_bound[1], local_track_bound[2]
    d = get_distance(S[1:2], point_R, point_L, S[4])
    # display(d)
    S_next[3] = sqrt(S[3]^2 + 2*U[2]*d)
    if U[2] == 0
        S_next[5] = d/S[3]
    else
        S_next[5] = (S_next[3]-S[3])/U[2]
    end
    S_next[1] = S[1] + d*cos(S[4])
    S_next[2] = S[2] + d*sin(S[4])
    S_next[4] = S[4] + U[1]
    return S_next
end 

# Compute State vector, S, given S[1] and Input vector, U
function compute_state(U, track_bound)
    N = lastindex(U[1,:])
    S = ones(5, N)
    # Assign initial state
    S[1,1] = track_bound[3][1][1]
    S[2,1] = track_bound[3][2][1]
    S[3,1] = 1
    S[4,1] = atan((track_bound[3][2][2]-track_bound[3][2][1])/(track_bound[3][1][2]-track_bound[3][1][1]+1e-6))
    S[5,1] = 0
    for n = 1:(N-1)
        local_track_bound = [[track_bound[1][1][n+1];track_bound[1][2][n+1]],[track_bound[2][1][n+1];track_bound[2][2][n+1]]]
        S[:,n+1] = semi_holonomic_model(S[:,n],U[:,n], local_track_bound)
    end
    return S
end

# Get distance from point
function get_distance(point,point_R,point_L,θ)
    x0, y0 = point[1], point[2]
    xR, yR = point_R[1], point_R[2]
    xL, yL = point_L[1], point_L[2]
    m = (yR-yL)/(xR-xL)
    c = -m*xR+yR
    d = (m*x0-y0+c)/(sin(θ)-m*cos(θ))
    return d
end

# Compute total time
function compute_total_time(S)
    return sum(S[5,:])
end

# Constraints
function compute_track_car_constraints(S, U, track_bound; ϵ=1e-5)
    coord = S[1:2,:]
    bound1, bound2 = reduce(hcat,track_bound[1])', reduce(hcat,track_bound[2])' # Convert vector of vector to array
    v, _ = S[3,:], S[4,:]
    θ, a = U[1,:], U[2,:]
    ω = θ./S[5,:]
    # Computes distance from each point to checkpoint lines
    N = lastindex(S[1, :])
    all_distances = [abs(norm(coord[:,n]-bound1[:,n])+norm(coord[:,n]-bound2[:,n])-norm(bound1[:,n]-bound2[:,n])) for n in 1:N]
    # Creates a vector with raw constraint penalties: size (N+1)*5 rows x 1 column
    track_car_constraints = [all_distances.-ϵ;
                            v./vmax.-1;
                            a./amax.-1;
                            abs.(θ)./θmax.-1;
                            abs.(ω)./ωmax.-1;
                            v.*ω./(Fc*m).-1]
    return track_car_constraints
end

# Initial design point (x0)
x0 = U

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

# Initialize n

n = 100

#######################################################################################################################################


using LinearAlgebra
using Plots
using Interpolations
using Random
include("track.jl")


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
        N = 20
        v_range = (0,100)
        population = initialize_population(x0, N, v_range)
        xhistory = particle_swarm_optimization(f_p, population, n; w=0.7, c1=1.2, c2=1.2)
        # display(xhistory)
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
    xhistory = copy(x_best) 
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
        append!(xhistory, x_best)
    end 
    return xhistory 
end


optimize(fun, con, x0, n)