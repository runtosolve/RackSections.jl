module RackSections

using CrossSection, RMI, LinesCurvesNodes, Parameters, CUFSM, AISIS100





@with_kw struct CeeLipsInput

    H::Float64
    D::Float64
    L::Float64
    R::Float64
    t::Float64
    E::Float64
    ν::Float64
    dh_H::Float64
    dh_D::Float64
    de_H::Float64
    de_D::Float64
    hole_pitch_H::Float64
    hole_pitch_D::Float64
    hole_length_H::Float64
    hole_length_D::Float64

end


@with_kw struct CeeLips

    input::CeeLipsInput
    section_info::NamedTuple{(:geometry, :x, :y, :D_hole_element_index, :H_hole_element_index), Tuple{NamedTuple{(:center, :left, :right), Tuple{Vector{Tuple{Float64, Float64}}, Vector{Tuple{Float64, Float64}}, Vector{Tuple{Float64, Float64}}}}, Vector{Float64}, Vector{Float64}, Vector{Int64}, Vector{Int64}}}
    properties::CUFSM.SectionPropertiesObject
    net_properties::CUFSM.SectionPropertiesObject
    Lnp_H::Float64
    Lnp_D::Float64
    tg_H::Float64
    tg_D::Float64
    tg::Vector{Float64}
    td_H::Float64
    td_D::Float64
    td::Vector{Float64}
    local_buckling_P::CUFSM.Model
    Pcrℓ::Float64
    local_buckling_Mxx::CUFSM.Model
    Mcrℓ_xx::Float64
    local_buckling_Mzz_neg::CUFSM.Model
    Mcrℓ_zz_neg::Float64
    local_buckling_Mzz_pos::CUFSM.Model
    Mcrℓ_zz_pos::Float64 
    distortional_buckling_P::CUFSM.Model
    Pcrd::Float64
    distortional_buckling_Mxx::CUFSM.Model
    Mcrd::Float64

end


@with_kw struct RectangularTubeInput

    H::Float64
    D::Float64
    R::Float64
    t::Float64
    E::Float64
    ν::Float64
    dh_H::Float64
    dh_D::Float64
    de_H::Float64
    de_D::Float64
    hole_pitch_H::Float64
    hole_pitch_D::Float64
    hole_length_H::Float64
    hole_length_D::Float64

end

@with_kw struct RectangularTube

    input::RectangularTubeInput
    section_info::NamedTuple{(:geometry, :x, :y, :D_hole_element_index, :H_hole_element_index), Tuple{NamedTuple{(:center, :left, :right), Tuple{Vector{Tuple{Float64, Float64}}, Vector{Tuple{Float64, Float64}}, Vector{Tuple{Float64, Float64}}}}, Vector{Float64}, Vector{Float64}, Vector{Int64}, Vector{Int64}}}
    properties::CUFSM.SectionPropertiesObject
    net_properties::CUFSM.SectionPropertiesObject
    Lnp_H::Float64
    Lnp_D::Float64
    tg_H::Float64
    tg_D::Float64
    tg::Vector{Float64}
    local_buckling_P::CUFSM.Model
    Pcrℓ::Float64
    local_buckling_Mxx::CUFSM.Model
    Mcrℓ_xx::Float64
    local_buckling_Mzz_neg::CUFSM.Model
    Mcrℓ_zz_neg::Float64
    local_buckling_Mzz_pos::CUFSM.Model
    Mcrℓ_zz_pos::Float64 

end


@with_kw struct StepBeamInput

    H::Float64
    D::Float64
    W::Float64
    L::Float64
    R::Float64
    t::Float64
    E::Float64
    ν::Float64

end


@with_kw struct StepBeam

    input::StepBeamInput
    section_info::NamedTuple{(:geometry, :x, :y), Tuple{NamedTuple{(:center, :left, :right), Tuple{Vector{Tuple{Float64, Float64}}, Vector{Tuple{Float64, Float64}}, Vector{Tuple{Float64, Float64}}}}, Vector{Float64}, Vector{Float64}}}
    properties::CUFSM.SectionPropertiesObject
    local_buckling_P::CUFSM.Model
    Pcrℓ::Float64
    local_buckling_Mxx::CUFSM.Model
    Mcrℓ_xx::Float64
    local_buckling_Mzz_neg::CUFSM.Model
    Mcrℓ_zz_neg::Float64
    local_buckling_Mzz_pos::CUFSM.Model
    Mcrℓ_zz_pos::Float64 

end



function cee_with_lips(H, D, L, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D)

    input = CeeLipsInput(H, D, L, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D)

    section_info = RackSections.cee_with_lips_geometry(H, D, L, R, t, dh_H, dh_D, de_H, de_D)

    #gross section properties 
    gross_section_properties = CrossSection.Properties.open_thin_walled(section_info.geometry.center, fill(t, length(section_info.x)-1)) 

    #calculate reduced thickness at hole
    Lnp_H = hole_pitch_H - hole_length_H
    Lnp_f = hole_pitch_D - hole_length_D

    kg = 0.60
    kd = 0.80

    tg_H = RMI.v2021.eq8_2__1(kg, t, Lnp_H, hole_pitch_H)
    tg_f = RMI.v2021.eq8_2__1(kg, t, Lnp_f, hole_pitch_D)

    td_H = RMI.v2021.eq8_2__2(kd, t, Lnp_H, hole_pitch_H)
    td_f = RMI.v2021.eq8_2__2(kd, t, Lnp_f, hole_pitch_D)

    #define cross-section element thicknesses
    num_elem = length(section_info.x) - 1
    tg = fill(t, num_elem)
    tg[section_info.D_hole_element_index] .= tg_f
    tg[section_info.H_hole_element_index] .= tg_H

    td = fill(t, num_elem)
    td[section_info.D_hole_element_index] .= td_f
    td[section_info.H_hole_element_index] .= td_H

    #net section properties 
    net_section_properties = CrossSection.Properties.open_thin_walled(section_info.geometry.center, tg) 


    #elastic buckling properties 

    #local buckling, compression
    P = 1.0
    Mxx = 0.0
    Mzz = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 5)
    local_buckling_P = CUFSM.Tools.open_section_analysis(section_info.x, section_info.y, tg, lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs)
    Pcrℓ = minimum(CUFSM.Tools.get_load_factor(local_buckling_P))

    #local buckling, Mxx 
    P = 0.0
    Mxx = 1.0
    Mzz = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 5)
    local_buckling_Mxx = CUFSM.Tools.open_section_analysis(section_info.x, section_info.y, tg, lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs)
    Mcrℓ_xx  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Mxx))


    #local buckling, Mzz_neg 
    P = 0.0
    Mxx = 0.0
    Mzz = 1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 5)
    local_buckling_Mzz_neg = CUFSM.Tools.open_section_analysis(section_info.x, section_info.y, tg, lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs)
    Mcrℓ_zz_neg  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Mzz_neg))

    #local buckling, Mzz_pos
    P = 0.0
    Mxx = 0.0
    Mzz = -1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.25*minimum([H, D]), 1.25*minimum([H,D]), 5)
    local_buckling_Mzz_pos = CUFSM.Tools.open_section_analysis(section_info.x, section_info.y, tg, lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs)
    Mcrℓ_zz_pos  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Mzz_pos))

    #distortional buckling 

    #find approximate distortional buckling half-wavelength
    CorZ = 0
    b = D
    d = L
    θ = 90.0
    Af,Jf,Ixf,Iyf,Ixyf,Cwf,xof,hxf,hyf,yof = AISIS100.v16.table23131(CorZ,t,b,d,θ)

    ho = H
    μ = ν
    Lm = 999999999.0
    Lcrd, not_used = AISIS100.v16.app23334(ho, μ, t, Ixf, xof, hxf, Cwf, Ixyf, Iyf, Lm)


    #distortional buckling, compression 
    P = 1.0
    Mxx = 0.0
    Mzz = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*Lcrd, 1.5*Lcrd, 9)
    distortional_buckling_P = CUFSM.Tools.open_section_analysis(section_info.x, section_info.y, td, lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs)
    Pcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_P))


    #distortional buckling, Mxx
    P = 0.0
    Mxx = 1.0
    Mzz = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*Lcrd, 2.0*Lcrd, 9)
    distortional_buckling_Mxx = CUFSM.Tools.open_section_analysis(section_info.x, section_info.y, td, lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs)
    Mcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_Mxx))


    #gather up everything 
    properties = CeeLips(input, section_info, gross_section_properties, net_section_properties, Lnp_H, Lnp_D, tg_H, tg_D, tg, td_H, td_f, td, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Mzz_neg, Mcrℓ_zz_neg, local_buckling_Mzz_pos, Mcrℓ_zz_pos, distortional_buckling_P, Pcrd, distortional_buckling_Mxx, Mcrd)

    return properties 

end




function cee_with_lips_geometry(H, D, L, R, t, dh_H, dh_D, de_H, de_D)

    segments = [L, D, H, D, L]
    θ = [π/2, π, -π/2, 0.0, π/2]
    r = [R, R, R, R]
    n = [4, 4, 5, 4, 4]
    n_r = [3, 3, 3, 3]

    section_geometry = CrossSection.Geometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D, H-L))

    x = [section_geometry.center[i][1] for i in eachindex(section_geometry.center)]
    y = [section_geometry.center[i][2] for i in eachindex(section_geometry.center)]

    nodes = [x y zeros(Float64, length(x))]


    #flange hole minus
    D_flat_minus = segments[4] - 2*R 
    xloc = t/2 + (R-t/2) + D_flat_minus/n[4] * 1.5
    yloc = t/2
    zloc = 0.0
    atol_x = D_flat_minus/n[4] * 0.55
    atol_y = 0.0
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_D_minus = hole_node_index[2] - 1

    x[hole_node_index[1]] = de_D - dh_D/2
    x[hole_node_index[2]] = de_D + dh_D/2

    #flange hole plus
    D_flat_plus = segments[2] - 2*R 
    xloc = t/2 + (R-t/2) + D_flat_plus/n[2] * 1.5
    yloc = H - t/2
    zloc = 0.0
    atol_x = D_flat_plus/n[2] * 0.55
    atol_y = 0.0
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_D_plus = hole_node_index[2] - 1

    x[hole_node_index[1]] = de_D + dh_D/2
    x[hole_node_index[2]] = de_D - dh_D/2

    D_hole_element_index = [hole_element_index_D_plus; hole_element_index_D_minus]


    #web hole plus
    H_flat = segments[3] - 2*R 
    xloc = t/2
    yloc = H - t/2 - (R-t/2) - H_flat/n[3] * 1.5
    zloc = 0.0
    atol_x = 0.0
    atol_y = H_flat/n[3] * 0.55
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_H_plus = hole_node_index[2] - 1

    y[hole_node_index[1]] = H - de_H + dh_H/2
    y[hole_node_index[2]] = H - de_H - dh_H/2

    #web hole minus
    xloc = t/2
    yloc = H - t/2 - (R-t/2) - H_flat/n[3] * 3.5
    zloc = 0.0
    atol_x = 0.0
    atol_y = H_flat/n[3] * 0.55
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_H_minus = hole_node_index[2] - 1

    y[hole_node_index[1]] = de_H + dh_H/2
    y[hole_node_index[2]] = de_H - dh_H/2

    H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]

    section_info = (geometry = section_geometry, x=x, y=y, D_hole_element_index = D_hole_element_index, H_hole_element_index=H_hole_element_index)

    return section_info

end




function rectangular_tube_geometry(H, D, R, t, dh_H, dh_D, de_H, de_D)


    segments = [H, D, H, D]
    θ = [π/2, π, -π/2, 0.0]
    r = [R, R, R, R]
    n = [4, 4, 5, 4]
    n_r = [3, 3, 3, 3]

    section_geometry = CrossSection.Geometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D, 0.0))

    x = [section_geometry.center[i][1] for i in eachindex(section_geometry.center)]
    y = [section_geometry.center[i][2] for i in eachindex(section_geometry.center)]

    nodes = [x y zeros(Float64, length(x))]


    #D hole minus
    D_flat_minus = segments[4] - 2*R 
    xloc = t/2 + (R-t/2) + D_flat_minus/n[4] * 1.5
    yloc = t/2
    zloc = 0.0
    atol_x = D_flat_minus/n[4] * 0.55
    atol_y = 0.0
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_D_minus = hole_node_index[2] - 1

    x[hole_node_index[1]] = de_D - dh_D/2
    x[hole_node_index[2]] = de_D + dh_D/2

    #D hole plus
    D_flat_plus = segments[2] - 2*R 
    xloc = t/2 + (R-t/2) + D_flat_plus/n[2] * 1.5
    yloc = H - t/2
    zloc = 0.0
    atol_x = D_flat_plus/n[2] * 0.55
    atol_y = 0.0
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_D_plus = hole_node_index[2] - 1

    x[hole_node_index[1]] = de_D + dh_D/2
    x[hole_node_index[2]] = de_D - dh_D/2

    D_hole_element_index = [hole_element_index_D_plus; hole_element_index_D_minus]


    #H hole plus
    H_flat = segments[3] - 2*R 
    xloc = t/2
    yloc = H - t/2 - (R-t/2) - H_flat/n[3] * 1.5
    zloc = 0.0
    atol_x = 0.0
    atol_y = H_flat/n[3] * 0.55
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_H_plus = hole_node_index[2] - 1

    y[hole_node_index[1]] = H - de_H + dh_H/2
    y[hole_node_index[2]] = H - de_H - dh_H/2

    #H hole minus
    xloc = t/2
    yloc = H - t/2 - (R-t/2) - H_flat/n[3] * 3.5
    zloc = 0.0
    atol_x = 0.0
    atol_y = H_flat/n[3] * 0.55
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_H_minus = hole_node_index[2] - 1

    y[hole_node_index[1]] = de_H + dh_H/2
    y[hole_node_index[2]] = de_H - dh_H/2

    H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]

    section_info = (geometry = section_geometry, x=x, y=y, D_hole_element_index = D_hole_element_index, H_hole_element_index=H_hole_element_index)

    return section_info

end



function rectangular_tube(H, D, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D)

    input = RectangularTubeInput(H, D, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D)

    section_info = rectangular_tube_geometry(H, D, R, t, dh_H, dh_D, de_H, de_D)

    #gross section properties 
    gross_section_properties = CrossSection.Properties.closed_thin_walled(section_info.geometry.center, fill(t, length(section_info.x))) 

    #calculate reduced thickness at hole
    Lnp_H = hole_pitch_H - hole_length_H
    Lnp_D = hole_pitch_D - hole_length_D

    kg = 0.60

    tg_H = RMI.v2021.eq8_2__1(kg, t, Lnp_H, hole_pitch_H)
    tg_D = RMI.v2021.eq8_2__1(kg, t, Lnp_D, hole_pitch_D)

    #define cross-section element thicknesses
    num_elem = length(section_info.x)
    tg = fill(t, num_elem)
    tg[section_info.D_hole_element_index] .= tg_D
    tg[section_info.H_hole_element_index] .= tg_H

    #net section properties 
    net_section_properties = CrossSection.Properties.closed_thin_walled(section_info.geometry.center, tg) 

    #local buckling, compression
    P = 1.0
    Mxx = 0.0
    Mzz = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 5)
    local_buckling_P = CUFSM.Tools.closed_section_analysis(section_info.x, section_info.y, tg, lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs)
    Pcrℓ = minimum(CUFSM.Tools.get_load_factor(local_buckling_P))

    #local buckling, Mxx 
    P = 0.0
    Mxx = 1.0
    Mzz = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 5)
    local_buckling_Mxx = CUFSM.Tools.closed_section_analysis(section_info.x, section_info.y, tg, lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs)
    Mcrℓ_xx  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Mxx))

    #local buckling, Mzz_neg 
    P = 0.0
    Mxx = 0.0
    Mzz = 1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 5)
    local_buckling_Mzz_neg = CUFSM.Tools.closed_section_analysis(section_info.x, section_info.y, tg, lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs)
    Mcrℓ_zz_neg  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Mzz_neg))

    #local buckling, Mzz_pos
    P = 0.0
    Mxx = 0.0
    Mzz = -1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.25*minimum([H, D]), 1.25*minimum([H,D]), 5)
    local_buckling_Mzz_pos = CUFSM.Tools.closed_section_analysis(section_info.x, section_info.y, tg, lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs)
    Mcrℓ_zz_pos  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Mzz_pos))

    #gather up everything 
    properties = RectangularTube(input, section_info, gross_section_properties, net_section_properties, Lnp_H, Lnp_D, tg_H, tg_D, tg, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Mzz_neg, Mcrℓ_zz_neg, local_buckling_Mzz_pos, Mcrℓ_zz_pos)

    return properties 

end





function step_beam_geometry(H, D, W, L, R, t)

    segments = [H-L, D-W, L, W, H, D]
    θ = [π/2, π, π/2,π, -π/2, 0.0]
    r = [R, R-t, R, R, R, R]
    n = [4, 4, 4, 4, 4, 4]
    n_r = [3, 3, 3, 3, 3, 3, 3]

    section_geometry = CrossSection.Geometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D, 0.0))

    x = [section_geometry.center[i][1] for i in eachindex(section_geometry.center)]
    y = [section_geometry.center[i][2] for i in eachindex(section_geometry.center)]

    section_info = (geometry = section_geometry, x=x, y=y)

    return section_info

end




function step_beam(H, D, W, L, R, t, E, ν)

    input = StepBeamInput(H, D, W, L, R, t, E, ν)

    section_info = step_beam_geometry(H, D, W, L, R, t)

    #gross section properties 
    gross_section_properties = CrossSection.Properties.closed_thin_walled(section_info.geometry.center, fill(t, length(section_info.x))) 

    
    #local buckling, compression
    P = 1.0
    Mxx = 0.0
    Mzz = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.5*H, 1.25*H, 5)
    local_buckling_P = CUFSM.Tools.closed_section_analysis(section_info.x, section_info.y, fill(t, length(section_info.x)), lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs)
    Pcrℓ = minimum(CUFSM.Tools.get_load_factor(local_buckling_P))

    #local buckling, Mxx 
    P = 0.0
    Mxx = 1.0
    Mzz = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.25*H, 1.0*H, 5)
    local_buckling_Mxx = CUFSM.Tools.closed_section_analysis(section_info.x, section_info.y, fill(t, length(section_info.x)), lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs)
    Mcrℓ_xx  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Mxx))

    #local buckling, Mzz_neg 
    P = 0.0
    Mxx = 0.0
    Mzz = 1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.5*H, 1.0*H, 5)
    local_buckling_Mzz_neg = CUFSM.Tools.closed_section_analysis(section_info.x, section_info.y, fill(t, length(section_info.x)), lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs)
    Mcrℓ_zz_neg  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Mzz_neg))

    #local buckling, Mzz_pos
    P = 0.0
    Mxx = 0.0
    Mzz = -1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.5*(H-L), 1.0*(H-L), 5)
    local_buckling_Mzz_pos = CUFSM.Tools.closed_section_analysis(section_info.x, section_info.y, fill(t, length(section_info.x)), lengths, E, ν, P, Mxx, Mzz, M11, M22, constraints, springs)
    Mcrℓ_zz_pos  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Mzz_pos))

    #gather up everything 
    properties = StepBeam(input, section_info, gross_section_properties, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Mzz_neg, Mcrℓ_zz_neg, local_buckling_Mzz_pos, Mcrℓ_zz_pos)

    return properties 

end


end # module RackSections
