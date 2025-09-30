using CrossSectionGeometry

R = 0.125 + 0.100
t = 0.100


A = 135.0

H = 4.0
D1 = 1.02 
D2 = 1.18 
D = D1 + D2 + 0.50
L = 0.39
dh_H = 0.710
dh_D1 = 0.431
dh_D2 = 0.531
de_H  = 0.8
de_D1 = D1/2
de_D2 = D2/2
hole_pitch_H = 2.0
hole_pitch_D1 = 2.0
hole_pitch_D2 = 2.0
hole_length_H = 1.086
hole_length_D1 = 0.531
hole_length_D2 = 0.531




A = deg2rad(A)  #convert 

    yc = t * cos(A/2)


    # #rib geometry, bottom of section  
    # I = π / 2   #assume this is how long the rib arc is, hard coded, always 

    # Tr = ((H_rib - t) - (R_rib_peak - t) * (1 - cos(I/2))) / sin(I/2)
    # T = (R_rib_peak - t) * tan(I/2)
    # T_total = Tr + T
    # ΔTx = T_total * cos(I/2)

    #calculate X 
    D3 = D - D2 - D1 
    X = (D3 + yc) * tan(π - A) + t

    segments = [L-t, D2-t-yc, (X-t)/sin(π-A), D1, H, D1, (X-t)/sin(π-A), D2-t-yc, L-t]

    θ = [-π/2, π, A, π, -π/2, 0.0, π-A, 0.0, -π/2]
    r = [R-t, R-t, R, R, R, R, R-t, R-t]
    n = [4, 3, 4, 3, 6, 3, 4, 3, 4]
    n_r = [3, 3, 3, 3, 3, 3, 3, 3] .* 3

    section_geometry = CrossSectionGeometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D-t, H - (X - L)))

    x = [section_geometry.centerline_node_XY[i][1] for i in eachindex(section_geometry.centerline_node_XY)]
    y = [section_geometry.centerline_node_XY[i][2] for i in eachindex(section_geometry.centerline_node_XY)]


    scatterlines(x, y)

    nodes = [x y zeros(Float64, length(x))]

    if (de_D1 !=0) & (dh_D1 !=0)

        #D1 flange hole minus
        hole_node_index_start = sum(n[1:5]) + 1 + sum(n_r[1:5]) + 1
        hole_node_index = [hole_node_index_start, hole_node_index_start + 1]
        hole_element_index_D1_minus = hole_node_index[2] - 1
        x[hole_node_index[1]] = de_D1 - dh_D1/2
        x[hole_node_index[2]] = de_D1 + dh_D1/2


        #D1 flange hole plus
        hole_node_index_start = sum(n[1:3]) + 1 + sum(n_r[1:3]) + 1
        hole_node_index = [hole_node_index_start, hole_node_index_start + 1]
        hole_element_index_D1_plus = hole_node_index[2] - 1
        x[hole_node_index[2]] = de_D1 - dh_D1/2
        x[hole_node_index[1]] = de_D1 + dh_D1/2

    else 
        hole_element_index_D1_plus = 0
        hole_element_index_D1_minus = 0
    end

    D1_hole_element_index = [hole_element_index_D1_plus; hole_element_index_D1_minus]


    if (de_D2 !=0) & (dh_D2 !=0)

        #D2 flange hole minus
        hole_node_index_start = sum(n[1:end-2]) + 1 + sum(n_r[1:end-1]) + 1
        hole_node_index = [hole_node_index_start, hole_node_index_start + 1]
        hole_element_index_D2_minus = hole_node_index[2] - 1
        x[hole_node_index[1]] = D - de_D2 - dh_D2/2
        x[hole_node_index[2]] = D - de_D2 + dh_D2/2

    
        #D2 flange hole plus
        hole_node_index = [n[1] + n_r[1] + 1 + 1, n[1] + n_r[1] + 2 + 1]
        hole_element_index_D2_plus = hole_node_index[2] - 1
        x[hole_node_index[2]] = D - de_D2 - dh_D2/2
        x[hole_node_index[1]] = D - de_D2 + dh_D2/2

    else

        hole_element_index_D2_plus = 0
        hole_element_index_D2_minus = 0

    end

    D2_hole_element_index = [hole_element_index_D2_plus; hole_element_index_D2_minus]


    if (de_H !=0) & (dh_H !=0)

        # H_flat = H - 2 * R

        #web hole plus
        hole_node_index_start = sum(n[1:4]) + 1 + sum(n_r[1:4]) + 1
        hole_node_index = [hole_node_index_start, hole_node_index_start+1]
        hole_element_index_H_plus = hole_node_index[2] - 1
        y[hole_node_index[1]] = H - de_H + dh_H/2
        y[hole_node_index[2]] = H - de_H - dh_H/2

        # xloc = t/2
        # yloc = 3*H/4
        # zloc = 0.0
        # atol_x = 0.0
        # atol_y = H_flat/n[5] * 1.05
        # atol_z = 0.0 

        # hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
        # hole_element_index_H_plus = hole_node_index[2] - 1

        # y[hole_node_index[1]] = de_H + dh_H/2
        # y[hole_node_index[2]] = de_H - dh_H/2


        # #web hole minus
        hole_node_index_start = sum(n[1:5]) + 1 + sum(n_r[1:4]) - 2
        hole_node_index = [hole_node_index_start, hole_node_index_start+1]
        hole_element_index_H_minus = hole_node_index[2] - 1
        y[hole_node_index[1]] = de_H + dh_H/2
        y[hole_node_index[2]] = de_H - dh_H/2


        #web hole minus
        # xloc = t/2
        # yloc = H/4
        # zloc = 0.0
        # atol_x = 0.0
        # atol_y = H_flat/n[5] * 1.05
        # atol_z = 0.0 

        # hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
        # hole_element_index_H_minus = hole_node_index[2] - 1

        # y[hole_node_index[1]] = de_H + dh_H/2
        # y[hole_node_index[2]] = de_H - dh_H/2



    else

        hole_element_index_H_plus = 0
        hole_element_index_H_minus = 0

    end


    H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]


    geometry = (coordinates = section_geometry, x=x, y=y, D1_hole_element_index = D1_hole_element_index, D2_hole_element_index = D2_hole_element_index, H_hole_element_index=H_hole_element_index)
