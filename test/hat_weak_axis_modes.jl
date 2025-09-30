using RackSections, RackSectionsAPI, ShowRackSections, CrossSectionGeometry, CairoMakie




member_type = "column"
section_type = "hat"
section_info = "with_holes_in_web"

E = 29500.0
ν = 0.30


R = 0.125 + 0.100
t = 0.100


A = 135.0

H = 4.0
D1 = 1.02 
D2 = 1.18 
D = D1 + D2 + 0.50
dh_H = 0.710
dh_D1 = 0.0
dh_D2 = 0.531
de_H  = 0.8
de_D1 = D1/2
de_D2 = D2/2
hole_pitch_H = 2.0
hole_pitch_D1 = 2.0
hole_pitch_D2 = 2.0
hole_length_H = 1.086
hole_length_D1 = 0.531
hole_length_D2 = 0.0

#######


        # A = deg2rad(A)  #convert 

        # yc = t * cos(A/2)
        

        # #calculate X 
        # D3 = D - D2 - D1 
        # X = (D3 + yc) * tan(π - A) + t
        
        # #define straight line segments 
        # segments = [D2-yc, (X-t)/sin(π-A), D1, H, D1, (X-t)/sin(π-A), D2-yc]
        
        # θ = [π, A, π, -π/2, 0.0, π-A, 0.0]
        # r = [R-t, R, R, R, R, R-t]
        # n = [3, 4, 3, 6, 3, 4, 3]
        # n_r = [5, 5, 5, 5, 5, 5]
        
        # section_geometry = CrossSectionGeometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D, H - (X-t)))
        
        # x = [section_geometry.centerline_node_XY[i][1] for i in eachindex(section_geometry.centerline_node_XY)]
        # y = [section_geometry.centerline_node_XY[i][2] for i in eachindex(section_geometry.centerline_node_XY)]
        
  
        # if (de_D1 != 0.0) & (dh_D1 != 0.0)

        #     #D1 flange hole minus
        #     hole_node_index_start = sum(n[1:7]) + 1 + sum(n_r[1:7]) + 1
        #     hole_node_index = [hole_node_index_start, hole_node_index_start + 1]
        #     hole_element_index_D1_minus = hole_node_index[2] - 1
        #     x[hole_node_index[1]] = de_D1 - dh_D1/2
        #     x[hole_node_index[2]] = de_D1 + dh_D1/2
            
            
        #     #D1 flange hole plus
        #     hole_node_index_start = sum(n[1:2]) + 1 + sum(n_r[1:2]) + 1
        #     hole_node_index = [hole_node_index_start, hole_node_index_start + 1]
        #     hole_element_index_D1_plus = hole_node_index[2] - 1
        #     x[hole_node_index[2]] = de_D1 - dh_D1/2
        #     x[hole_node_index[1]] = de_D1 + dh_D1/2

        # else 

        #     hole_element_index_D1_plus = 0
        #     hole_element_index_D1_minus = 0

        # end

        
        # D1_hole_element_index = [hole_element_index_D1_plus; hole_element_index_D1_minus]
        
        
        # if (de_D2 != 0.0) & (dh_D2 != 0.0)

        #     #D2 flange hole minus
        #     hole_node_index_start = sum(n[1:end-1]) + 1 + sum(n_r[1:end]) + 1
        #     hole_node_index = [hole_node_index_start, hole_node_index_start + 1]
        #     hole_element_index_D2_minus = hole_node_index[2] - 1
        #     x[hole_node_index[1]] = D - de_D2 - dh_D2/2
        #     x[hole_node_index[2]] = D - de_D2 + dh_D2/2
            
        #     #D2 flange hole plus
        #     hole_node_index = [2, 3]
        #     hole_element_index_D2_plus = hole_node_index[2] - 1
        #     x[hole_node_index[2]] = D - de_D2 - dh_D2/2
        #     x[hole_node_index[1]] = D - de_D2 + dh_D2/2

        # else

        #     hole_element_index_D2_plus = 0
        #     hole_element_index_D2_minus = 0

        # end
        
        # D2_hole_element_index = [hole_element_index_D2_plus; hole_element_index_D2_minus]
        

        # if (de_H != 0.0) & (dh_H != 0.0)

        #     #web hole plus
        #     hole_node_index_start = sum(n[1:3]) + 1 + sum(n_r[1:3]) + 1
        #     hole_node_index = [hole_node_index_start, hole_node_index_start+1]
        #     hole_element_index_H_plus = hole_node_index[2] - 1
        #     y[hole_node_index[1]] = H - de_H + dh_H/2
        #     y[hole_node_index[2]] = H - de_H - dh_H/2
            
        #     #web hole minus
        #     hole_node_index_start = sum(n[1:4]) + sum(n_r[1:3]) - 1
        #     hole_node_index = [hole_node_index_start, hole_node_index_start+1]
        #     hole_element_index_H_minus = hole_node_index[2] - 1
        #     y[hole_node_index[1]] = de_H + dh_H/2
        #     y[hole_node_index[2]] = de_H - dh_H/2

        # else

        #     hole_element_index_H_plus = 0
        #     hole_element_index_H_minus = 0

        # end

        
        # H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]
        










######


section_details = RackSections.Columns.HatInput(H, D1, D2, D, A, R, t, E, ν, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2, hole_pitch_H, hole_pitch_D1, hole_pitch_D2, hole_length_H, hole_length_D1, hole_length_D2)

properties = RackSections.Columns.hat(section_details)

api_figure_options = (max_pixel_size = 2048, cross_section_linecolor =:grey, signature_curve_linecolor=:blue)


all_figures_IO, all_figures, figure_labels = ShowRackSections.hat_with_rib_column(properties, api_figure_options)

# api_inputs = RackSectionsAPI.Inputs(member_type, section_type, section_info, section_details, create_output_binary, CUFSM_figure_files_bucket_name, create_CUFSM_MAT_files, CUFSM_MAT_files_bucket_name)
# event_data = JSON3.write(api_inputs)
# section_outputs = RackSectionsAPI.handle_event(event_data, String[])

# write_input_output_jsons(JSON_file_path, member_type, section_type, section_info, api_inputs, section_outputs)


