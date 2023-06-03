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
    # display(U[2])
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