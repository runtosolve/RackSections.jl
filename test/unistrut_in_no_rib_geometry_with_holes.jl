using RackSections, CrossSectionGeometry, CairoMakie, LinesCurvesNodes 

H = 3.0 
D = 3.0
L1 = 0.75
L2 = 0.75
R = 0.125 + 0.100
t = 0.100
dh_H = 0.710
dh_D = 0.531
de_H  = H/2 - 0.700
de_D = 0.875



geometry = RackSections.Columns.unistrut_in_no_rib_geometry(H, D, L1, L2, R, t, dh_H, dh_D, de_H, de_D)


# ######




#     segments = [L2, L1, D, H, D, L1, L2]

#     θ = [0.0, π/2, π, -π/2, 0.0, π/2, π]
#     r = [R, R, R, R, R, R]
#     n = [4, 4, 3, 6, 3, 4, 4]
#     n_r = [5, 5, 5, 5, 5, 5]

#     section_geometry = CrossSectionGeometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D-L2, H-L1))

#     x = [section_geometry.centerline_node_XY[i][1] for i in eachindex(section_geometry.centerline_node_XY)]
#     y = [section_geometry.centerline_node_XY[i][2] for i in eachindex(section_geometry.centerline_node_XY)]

#     scatterlines(x, y)

#     nodes = [x y zeros(Float64, length(x))]

#     # if (de_D != 0.0) & (dh_D != 0)


#         #flange hole minus
#         D_flat = D - 2*R 
#         xloc = D/2
#         yloc = t/2
#         zloc = 0.0
#         atol_x = D_flat/3/2 * 1.05
#         atol_y = 0.0
#         atol_z = 0.0 

#         hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
#         hole_element_index_D_minus = hole_node_index[2] - 1

#         x[hole_node_index[1]] = de_D - dh_D/2
#         x[hole_node_index[2]] = de_D + dh_D/2

#         scatterlines(x, y)

#         #flange hole plus
#         # D_flat_plus = D - 2*R 
#         xloc = D/2
#         yloc = H - t/2
#         zloc = 0.0
#         atol_x = D_flat/3/2 * 1.05
#         atol_y = 0.0
#         atol_z = 0.0 

#         hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
#         hole_element_index_D_plus = hole_node_index[2] - 1

#         x[hole_node_index[1]] = de_D + dh_D/2
#         x[hole_node_index[2]] = de_D - dh_D/2

#         scatterlines(x, y)

#     # else 

#     #     hole_element_index_D_plus = 0
#     #     hole_element_index_D_minus = 0

#     # end


#     D_hole_element_index = [hole_element_index_D_plus; hole_element_index_D_minus]

#     #web hole plus

#     # if (de_H != 0.0) & (dh_H != 0)

#         H_flat = H - 2*R

#         # index_start = n[1] + n_r[1] + n[2] + n_r[2] + n[3] + n_r[3] + 1 + 1
#         # index_end = n[1] + n_r[1] + n[2] + n_r[2] + n[3] + n_r[3] + 2 + 1

#         # centerline_rib_length = y[index_start] - y[index_end]

#         xloc = t/2
        
#         yloc = R + 3/4 * H_flat 
#         zloc = 0.0
#         atol_x = 0.0
#         atol_y = H_flat / 2 / 3 / 2 * 1.05
#         atol_z = 0.0 

#         hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
#         hole_element_index_H_plus = hole_node_index[2] - 1

#         y[hole_node_index[1]] = H - de_H + dh_H/2
#         y[hole_node_index[2]] = H - de_H - dh_H/2

#         scatterlines(x, y)

#         #web hole minus
#         xloc = t/2
#         yloc = R + 1/4 * H_flat
#         zloc = 0.0
#         atol_x = 0.0
#         atol_y = H_flat / 2 / 3/ 2 * 1.05
#         atol_z = 0.0 

#         hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
#         hole_element_index_H_minus = hole_node_index[2] - 1

#         y[hole_node_index[1]] = de_H + dh_H/2
#         y[hole_node_index[2]] = de_H - dh_H/2
#         scatterlines(x, y)

#     else 

#         hole_element_index_H_plus = 0
#         hole_element_index_H_minus = 0

#     end

#     H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]

#     geometry = (coordinates = section_geometry, x=x, y=y, D_hole_element_index = D_hole_element_index, H_hole_element_index=H_hole_element_index)
