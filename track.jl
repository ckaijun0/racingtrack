using LinearAlgebra
using Plots
using Interpolations
using Random

x = [1; 2; 3; 3; 2; 1; 0; 0]
y = [0; 0; 1; 2; 3; 3; 2; 1]
middle_line = [x,y]

# Creates the width given an input set of points 
function create_track_width(middle_line, width)
    xm, ym = middle_line[1], middle_line[2]
    vector_length = lastindex(xm)
    push!(xm,xm[1])
    push!(ym,ym[1])
    xo, yo, xi, yi = zeros(vector_length), zeros(vector_length), zeros(vector_length), zeros(vector_length)
    for i = 1:vector_length
        grad = [(xm[i+1] - xm[i]);(ym[i+1] - ym[i])]
        ortho = [-grad[2];grad[1]]
        xo[i] = xm[i]+(ortho/norm(ortho)*width)[1]
        xi[i] = xm[i]-(ortho/norm(ortho)*width)[1]
        yo[i] = ym[i]+(ortho/norm(ortho)*width)[2]
        yi[i] = ym[i]-(ortho/norm(ortho)*width)[2]        
    end
    push!(xo,xo[1])
    push!(yo,yo[1])
    push!(xi,xi[1])
    push!(yi,yi[1])
    return [xo,yo], [xi,yi]
end

outer_line, inner_line = create_track_width(middle_line, 0.1)

plot(middle_line[1],middle_line[2])
plot!(outer_line[1],outer_line[2])
plot!(inner_line[1],inner_line[2])
savefig("track.png")