using LinearAlgebra
using Plots
using Interpolations
using Random

# Specify the centerline track coordinates here (don't duplicate start and end points)
# x = [1.0; 2; 3; 3; 2; 1; 0; 0]
# y = [0.0; 0; 1.0; 2; 3; 3; 2; 1]
# middle_line = [x,y]

# Artificially extend the start and ends of the track to permit initial and final change in heading angles
function artificial_zones(xa, ya)
    artificial_zone_percent = 0.01 
    end_zone_x = xa[end] + (xa[end] - xa[end-1])*artificial_zone_percent
    end_zone_y = ya[end] + (ya[end] - ya[end-1])*artificial_zone_percent
    start_zone_x = xa[1] + (xa[1] - xa[2])*artificial_zone_percent
    start_zone_y = ya[1] + (ya[2] - ya[1])*artificial_zone_percent
    return end_zone_x, end_zone_y, start_zone_x, start_zone_y
end

# Creates the width given an input set of points 
function create_track_width(middle_line, width, complete_loop)
    xm, ym = copy(middle_line[1]), copy(middle_line[2])
    xmr, ymr = copy(middle_line[1]), copy(middle_line[2])
    vector_length = lastindex(xm)

    # Copy the start and end points to the end and start of the vector respectively (if complete_loop = true)
    # Otherwise create artificial short zones for the end and start (if complete_loop = false)
    if complete_loop == true
        reverse!(push!(reverse!(push!(xm,xm[1])),xm[2]))
        reverse!(push!(reverse!(push!(ym,ym[1])),ym[2]))
    elseif complete_loop == false
        end_zone_xm, end_zone_ym, start_zone_xm, start_zone_ym = artificial_zones(xm, ym)
        reverse!(push!(xm,end_zone_xm))
        reverse!(push!(ym,end_zone_ym))
        reverse!(push!(xm,start_zone_xm))
        reverse!(push!(ym,start_zone_ym))       
    end

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
    if complete_loop == true
        push!(xo,xo[1])
        push!(yo,yo[1])
        push!(xi,xi[1])
        push!(yi,yi[1])
        push!(xmr,xmr[1])
        push!(ymr,ymr[1])
    elseif complete_loop == false
        end_zone_xo, end_zone_yo, start_zone_xo, start_zone_yo = artificial_zones(xo,yo)
        push!(xo,end_zone_xo)
        push!(yo,end_zone_yo)
        end_zone_xi, end_zone_yi, start_zone_xi, start_zone_yi = artificial_zones(xi,yi)
        push!(xi,end_zone_xi)
        push!(yi,end_zone_yi)
        end_zone_xmr, end_zone_ymr, start_zone_xmr, start_zone_ymr = artificial_zones(xmr,ymr)
        push!(xmr,end_zone_xmr)
        push!(ymr,end_zone_ymr)        
    end
    return [xo,yo], [xi,yi], [xmr,ymr]
end

# Compute euclidean distance between 2 points
function dist(point1,point2)
    return norm(point1-point2)
end

# outer_line, inner_line, center_line = create_track_width(middle_line, 0.1)

# plt = plot(center_line[1],center_line[2])
# plot!(outer_line[1],outer_line[2])
# plot!(inner_line[1],inner_line[2])
# savefig("track.png")
# display(plt)