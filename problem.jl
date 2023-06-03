using LinearAlgebra
using Plots
using Interpolations
using Random
include("track.jl")
include("plotting.jl")

function build_track()
    # Specify the centerline track coordinates here (don't duplicate start and end points)
    # x_points = [1; 2; 3; 3; 2; 1; 0; 0]
    # y_points = [0; 0; 1; 2; 3; 3; 2; 1]
    x_points = [1; 1.5; 3; 3.1; 3.2; 3.3; 3.4; 3.5; 4; 4.1; 3.9; 0; 0.5; 1]
    y_points = [0; 0.5; 0; -0.2; -0.3; -0.35; -0.3; -0.2; 0; 2; 2.5; 2.5; 0.5; 0.25]
    width = 0.1

    center_line = [x_points, y_points]
    track_bound = create_track_width(center_line, width)
    return track_bound
end

function get_angle(x0,y0,x1,y1)
    return atan((y1-y0),(x1-x0))
end

# Create input variables vector
function create_U0(track_bound)
    center_line = track_bound[3]
    N = lastindex(track_bound[1][1])
    U = zeros(2, N)
    heading_angle = [get_angle(center_line[1][n],center_line[2][n],center_line[1][n+1],center_line[2][n+1]) for n = 1:(N-1)]
    U[1,1:(N-2)] = [heading_angle[n+1] - heading_angle[n] for n in 1:(N-2)]
    U[1,(N-1)] = heading_angle[1] - heading_angle[N-1]
    U[2,:] = [0.05*n for n=1:N]
    return U
end

# Create the semi holonomic model, to calculate S+1 from S and U
function semi_holonomic_model(S, U, local_track_bound)
    S_next = copy(S)
    point_R, point_L = local_track_bound[1], local_track_bound[2]
    d = get_distance(S[1:2], point_R, point_L, S[4])
    if (S[3]^2 + 2*U[2]*abs(d)) < 0
        S_next[3] = -sqrt(-(S[3]^2 + 2*U[2]*abs(d)))
    else
        S_next[3] = sqrt(S[3]^2 + 2*U[2]*abs(d))
    end
    if U[2] == 0
        S_next[5] = abs(d)/S[3]
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
    S[4,1] = atan((track_bound[3][2][2]-track_bound[3][2][1]),(track_bound[3][1][2]-track_bound[3][1][1]))
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
    if xR-xL == 0
        m = 1e12 # some abitrarily large number
    else
        m = (yR-yL)/(xR-xL)
    end
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
    # Parameters for user input
    vmax = 50 # m/s (111.8mph)
    amax = 3 # m/s^2
    θmax = 90 # degrees 
    ωmax = 15 # ω=θ/t (assuming full lock in 3seconds)
    m = 1500 # kg
    friction_coefficient = 0.7 # 0.7=dry road, 0.4=wet road
    ϵ = 1e-5 # permitted error on position due to floating point errors

    # Derived parameters
    N = lastindex(x_points)+1 # number of points on track + 1 (repeats start point)
    Fc = friction_coefficient*m*9.81 # friction = f(normal_force)

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
    #display(total_time)
    return total_time
end

# track_bound = build_track()
# U = create_U0(track_bound)
# S = compute_state(U, track_bound)
# print(S)
# Plot track (function is in plotting.jl)
# plot_track(S, U, track_bound)


