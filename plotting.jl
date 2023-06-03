using LinearAlgebra
using Plots

# # Generate example data
# x = [1; 2; 3; 3; 2; 1; 0; 0; 1]
# y = [0; 0; 1; 2; 3; 3; 2; 1; 0]
# a = [1; 2; 0.5; -1; 1; 0.5; -0.5; 1]  # acceleration
# v = [sum(a[1:n]) for n = 1:8]         # velocity (scalar)
# θ = [pi/4*(n-1) for n = 1:8]          # heading angle

# # Plotting
# arrow_scale = 0.2
# plt = quiver(x, y, quiver=(arrow_scale.*v.*cos.(θ), arrow_scale.*v.*sin.(θ)), 
#              c = :black, linewidth = 2, label = "velocity [m/s]", legend = false)
# scatter!(x[1:(end-1)], y[1:(end-1)], markersize = 8, marker_z = a, colorbar = true, 
#          c = cgrad(:RdYlGn), label = "acceleration [m/s2]")
# plot!(track_bound[1][1],track_bound[1][2], color = :darkblue, label="track bound")
# plot!(track_bound[2][1],track_bound[2][2], color = :darkblue, label="track bound")
# display(plt)

function plot_track(S, U, track_bound)
    x = S[1,1:(end-1)]
    y = S[2,1:(end-1)]
    v = S[3,1:(end-1)]
    θ = S[4,1:(end-1)]
    a = U[2,1:(end-1)]

    # Plotting
    arrow_scale = 0.2
    plt = scatter(x, y, markersize = 8, marker_z = a, colorbar = true, 
            c = cgrad(:RdYlGn), label = "acceleration [m/s2]")
    quiver!(x, y, quiver=(arrow_scale.*v.*cos.(θ), arrow_scale.*v.*sin.(θ)), 
            c = :black, linewidth = 2, label = "velocity [m/s]", legend = false)
    plot!(track_bound[1][1],track_bound[1][2], color = :darkblue, label="track bound")
    plot!(track_bound[2][1],track_bound[2][2], color = :darkblue, label="track bound")
    display(plt)
end