using LinearAlgebra
using Plots
using Interpolations
using Random

# Specify the centerline track coordinates here (don't duplicate start and end points)
x = [1; 2; 3; 3; 2; 1; 0; 0]
y = [0; 0; 1; 2; 3; 3; 2; 1]
middle_line = [x,y]

# Creates the width given an input set of points 
function create_track_width(middle_line, width)
    xm, ym = copy(middle_line[1]), copy(middle_line[2])
    xmr, ymr = copy(middle_line[1]), copy(middle_line[2])
    vector_length = lastindex(xm)

    # Copy the start and end points to the end and start of the vector respectively
    reverse!(push!(reverse!(push!(xm,xm[1])),xm[2]))
    reverse!(push!(reverse!(push!(ym,ym[1])),ym[2]))

    # Initialize vector for bound coordinates
    xo, yo, xi, yi = zeros(vector_length), zeros(vector_length), zeros(vector_length), zeros(vector_length)

    # Create the bounds: prev-from previous point to current point, ahead-from current point to next point
    for i = 2:(vector_length+1)
        grad_prev = [(xm[i+1] - xm[i]);(ym[i+1] - ym[i])]
        grad_ahead = [(xm[i] - xm[i-1]);(ym[i] - ym[i-1])]
        
        # Get orthogonal vector to centerline
        ortho_prev = [-grad_prev[2];grad_prev[1]]
        ortho_ahead = [-grad_ahead[2];grad_ahead[1]]
        ortho = ortho_prev/norm(ortho_prev) + ortho_ahead/norm(ortho_ahead)
        xo[i-1] = xm[i]+(ortho/norm(ortho)*width)[1]
        xi[i-1] = xm[i]-(ortho/norm(ortho)*width)[1]
        yo[i-1] = ym[i]+(ortho/norm(ortho)*width)[2]
        yi[i-1] = ym[i]-(ortho/norm(ortho)*width)[2]        
    end
    push!(xo,xo[1])
    push!(yo,yo[1])
    push!(xi,xi[1])
    push!(yi,yi[1])
    push!(xmr,xmr[1])
    push!(ymr,ymr[1])
    return [xo,yo], [xi,yi], [xmr,ymr]
end

function dist(point1,point2)
    return norm(point1-point2)
end

outer_line, inner_line, center_line = create_track_width(middle_line, 0.1)

plot(center_line[1],center_line[2])
plot!(outer_line[1],outer_line[2])
plot!(inner_line[1],inner_line[2])
savefig("track.png")