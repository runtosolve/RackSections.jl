using RackSections, RackSectionsAPI, ShowRackSections


E = 29500.0
ν = 0.30


member_type = "column"
section_type = "cee_with_lips"
section_info = "with_holes_in_web_and_flanges"
H = 3.0 
D = 3.0
L = 0.75
R = 0.125 + 0.100
t = 0.100
dh_H = 0.710
dh_D = 0.531
de_H  = H/2 - 0.700
de_D = 0.875
hole_pitch_H = 2.0
hole_pitch_D = 2.0
hole_length_H = 1.086
hole_length_D = 0.531;

# create_output_binary = true 
# CUFSM_figure_files_bucket_name = joinpath("epiq-cufsm-files", member_type, section_type, section_info, "CUFSM_figures")

# create_CUFSM_MAT_files = true
# CUFSM_MAT_files_bucket_name = joinpath("epiq-cufsm-files", member_type, section_type, section_info, "CUFSM_MAT")


section_details = RackSections.Columns.CeeLipsInput(H, D, L, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D)

properties = RackSections.Columns.cee_with_lips(section_details)

api_figure_options = (max_pixel_size = 2048, cross_section_linecolor =:grey, signature_curve_linecolor=:blue)


all_figures_IO, all_figures, figure_labels = ShowRackSections.cee_with_lips_column(properties, api_figure_options)

# api_inputs = RackSectionsAPI.Inputs(member_type, section_type, section_info, section_details, create_output_binary, CUFSM_figure_files_bucket_name, create_CUFSM_MAT_files, CUFSM_MAT_files_bucket_name)
# event_data = JSON3.write(api_inputs)
# section_outputs = RackSectionsAPI.handle_event(event_data, String[])

# write_input_output_jsons(JSON_file_path, member_type, section_type, section_info, api_inputs, section_outputs)
