module Columns 

using CrossSection, RMI, LinesCurvesNodes, Parameters, CUFSM, AISIS100, LinearAlgebra

@with_kw struct CeeLipsInput

    member_type::String
    section_type::String
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
    geometry::@NamedTuple{coordinates::@NamedTuple{center::Vector{Vector{Float64}}, left::Vector{Vector{Float64}}, right::Vector{Vector{Float64}}}, x::Vector{Float64}, y::Vector{Float64}, D_hole_element_index::Vector{Int64}, H_hole_element_index::Vector{Int64}}
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
    local_buckling_Myy_neg::CUFSM.Model
    Mcrℓ_yy_neg::Float64
    local_buckling_Myy_pos::CUFSM.Model
    Mcrℓ_yy_pos::Float64 
    distortional_buckling_P::CUFSM.Model
    Pcrd::Float64
    distortional_buckling_Mxx::CUFSM.Model
    Mcrd::Float64

end



@with_kw struct CeeLipsRibInput

    member_type::String
    section_type::String
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
    rib_depth::Float64
    rib_length::Float64
    rib_radius_start::Float64
    rib_radius_peak::Float64

end



@with_kw struct CeeLipsRib

    input::CeeLipsRibInput
    geometry::@NamedTuple{coordinates::@NamedTuple{center::Vector{Vector{Float64}}, left::Vector{Vector{Float64}}, right::Vector{Vector{Float64}}}, x::Vector{Float64}, y::Vector{Float64}, D_hole_element_index::Vector{Int64}, H_hole_element_index::Vector{Int64}}
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
    local_buckling_Myy_neg::CUFSM.Model
    Mcrℓ_yy_neg::Float64
    local_buckling_Myy_pos::CUFSM.Model
    Mcrℓ_yy_pos::Float64 
    distortional_buckling_P::CUFSM.Model
    Pcrd::Float64
    distortional_buckling_Mxx::CUFSM.Model
    Mcrd::Float64

end


@with_kw struct RectangularTubeInput

    member_type::String
    section_type::String
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
    geometry::@NamedTuple{coordinates::@NamedTuple{center::Vector{Vector{Float64}}, left::Vector{Vector{Float64}}, right::Vector{Vector{Float64}}}, x::Vector{Float64}, y::Vector{Float64}, D_hole_element_index::Vector{Int64}, H_hole_element_index::Vector{Int64}}
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
    local_buckling_Myy_neg::CUFSM.Model
    Mcrℓ_yy_neg::Float64
    local_buckling_Myy_pos::CUFSM.Model
    Mcrℓ_yy_pos::Float64 

end

@with_kw struct HatLipsRibInput

    member_type::String
    section_type::String
    H::Float64
    D1::Float64
    D2::Float64
    D3::Float64
    A::Float64
    X::Float64
    L::Float64
    R::Float64
    t::Float64
    E::Float64
    ν::Float64
    dh_H::Float64
    dh_D1::Float64
    dh_D2::Float64
    de_H::Float64
    de_D1::Float64
    de_D2::Float64
    hole_pitch_H::Float64
    hole_pitch_D1::Float64
    hole_pitch_D2::Float64
    hole_length_H::Float64
    hole_length_D1::Float64
    hole_length_D2::Float64
    rib_depth::Float64
    rib_length::Float64
    rib_radius_start::Float64
    rib_radius_peak::Float64

end


@with_kw struct HatLipsRib

    input::HatLipsRibInput
    geometry::@NamedTuple{coordinates::@NamedTuple{center::Vector{Vector{Float64}}, left::Vector{Vector{Float64}}, right::Vector{Vector{Float64}}}, x::Vector{Float64}, y::Vector{Float64}, D1_hole_element_index::Vector{Int64}, D2_hole_element_index::Vector{Int64}, H_hole_element_index::Vector{Int64}}
    properties::CUFSM.SectionPropertiesObject
    net_properties::CUFSM.SectionPropertiesObject
    Lnp_H::Float64
    Lnp_D1::Float64
    Lnp_D2::Float64
    tg_H::Float64
    tg_D1::Float64
    tg_D2::Float64
    tg::Vector{Float64}
    td_H::Float64
    td_D1::Float64
    td_D2::Float64
    td::Vector{Float64}
    local_buckling_P::CUFSM.Model
    Pcrℓ::Float64
    local_buckling_Mxx::CUFSM.Model
    Mcrℓ_xx::Float64
    local_buckling_Myy_neg::CUFSM.Model
    Mcrℓ_yy_neg::Float64
    distortional_buckling_Myy_pos::CUFSM.Model
    Mcrd_yy_pos::Float64 
    distortional_buckling_P::CUFSM.Model
    Pcrd::Float64
    distortional_buckling_Mxx::CUFSM.Model
    Mcrd::Float64

end


@with_kw struct HatRibInput

    member_type::String
    section_type::String
    H::Float64
    D1::Float64
    D2::Float64
    D3::Float64
    A::Float64
    X::Float64
    R::Float64
    t::Float64
    E::Float64
    ν::Float64
    dh_H::Float64
    dh_D1::Float64
    dh_D2::Float64
    de_H::Float64
    de_D1::Float64
    de_D2::Float64
    hole_pitch_H::Float64
    hole_pitch_D1::Float64
    hole_pitch_D2::Float64
    hole_length_H::Float64
    hole_length_D1::Float64
    hole_length_D2::Float64
    rib_depth::Float64
    rib_length::Float64
    rib_radius_start::Float64
    rib_radius_peak::Float64

end


@with_kw struct HatRib

    input::HatRibInput
    geometry::@NamedTuple{coordinates::@NamedTuple{center::Vector{Vector{Float64}}, left::Vector{Vector{Float64}}, right::Vector{Vector{Float64}}}, x::Vector{Float64}, y::Vector{Float64}, D1_hole_element_index::Vector{Int64}, D2_hole_element_index::Vector{Int64}, H_hole_element_index::Vector{Int64}}
    properties::CUFSM.SectionPropertiesObject
    net_properties::CUFSM.SectionPropertiesObject
    Lnp_H::Float64
    Lnp_D1::Float64
    Lnp_D2::Float64
    tg_H::Float64
    tg_D1::Float64
    tg_D2::Float64
    tg::Vector{Float64}
    td_H::Float64
    td_D1::Float64
    td_D2::Float64
    td::Vector{Float64}
    local_buckling_P::CUFSM.Model
    Pcrℓ::Float64
    local_buckling_Mxx::CUFSM.Model
    Mcrℓ_xx::Float64
    local_buckling_Myy_neg::CUFSM.Model
    Mcrℓ_yy_neg::Float64
    distortional_buckling_Myy_pos::CUFSM.Model
    Mcrd_yy_pos::Float64 
    distortional_buckling_P::CUFSM.Model
    Pcrd::Float64
    distortional_buckling_Mxx::CUFSM.Model
    Mcrd::Float64

end


@with_kw struct HatLipsTrapezoidalRibInput

    member_type::String
    section_type::String
    H::Float64
    D1::Float64
    D2::Float64
    D3::Float64
    A1::Float64
    X::Float64
    L::Float64
    R::Float64
    t::Float64
    E::Float64
    ν::Float64
    dh_H::Float64
    dh_D1::Float64
    dh_D2::Float64
    de_H::Float64
    de_D1::Float64
    de_D2::Float64
    hole_pitch_H::Float64
    hole_pitch_D1::Float64
    hole_pitch_D2::Float64
    hole_length_H::Float64
    hole_length_D1::Float64
    hole_length_D2::Float64
    A2::Float64
    hr::Float64
    wr::Float64
    Rr::Float64

end



@with_kw struct HatLipsTrapezoidalRib

    input::HatLipsTrapezoidalRibInput
    geometry::@NamedTuple{coordinates::@NamedTuple{center::Vector{Vector{Float64}}, left::Vector{Vector{Float64}}, right::Vector{Vector{Float64}}}, x::Vector{Float64}, y::Vector{Float64}, D1_hole_element_index::Vector{Int64}, D2_hole_element_index::Vector{Int64}, H_hole_element_index::Vector{Int64}}
    properties::CUFSM.SectionPropertiesObject
    net_properties::CUFSM.SectionPropertiesObject
    Lnp_H::Float64
    Lnp_D1::Float64
    Lnp_D2::Float64
    tg_H::Float64
    tg_D1::Float64
    tg_D2::Float64
    tg::Vector{Float64}
    td_H::Float64
    td_D1::Float64
    td_D2::Float64
    td::Vector{Float64}
    local_buckling_P::CUFSM.Model
    Pcrℓ::Float64
    local_buckling_Mxx::CUFSM.Model
    Mcrℓ_xx::Float64
    local_buckling_Myy_neg::CUFSM.Model
    Mcrℓ_yy_neg::Float64
    distortional_buckling_Myy_pos::CUFSM.Model
    Mcrd_yy_pos::Float64 
    distortional_buckling_P::CUFSM.Model
    Pcrd::Float64
    distortional_buckling_Mxx::CUFSM.Model
    Mcrd::Float64

end



@with_kw struct UniStrutInput

    member_type::String
    section_type::String
    H::Float64
    D::Float64
    L1::Float64
    L2::Float64
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
    rib_depth::Float64
    rib_length::Float64
    rib_radius_start::Float64
    rib_radius_peak::Float64

end


@with_kw struct UniStrut

    input::UniStrutInput
    geometry::@NamedTuple{coordinates::@NamedTuple{center::Vector{Vector{Float64}}, left::Vector{Vector{Float64}}, right::Vector{Vector{Float64}}}, x::Vector{Float64}, y::Vector{Float64}, D_hole_element_index::Vector{Int64}, H_hole_element_index::Vector{Int64}}
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
    local_buckling_Myy_neg::CUFSM.Model
    Mcrℓ_yy_neg::Float64
    local_buckling_Myy_pos::CUFSM.Model
    Mcrℓ_yy_pos::Float64 
    distortional_buckling_P::CUFSM.Model
    Pcrd::Float64
    distortional_buckling_Mxx::CUFSM.Model
    Mcrd::Float64

end


function cee_with_lips(section_inputs)

    @unpack member_type, section_type, H, D, L, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D = section_inputs

    geometry = cee_with_lips_geometry(H, D, L, R, t, dh_H, dh_D, de_H, de_D)

    #gross section properties 
    gross_section_properties = CrossSection.Properties.open_thin_walled(geometry.coordinates.center, fill(t, length(geometry.x)-1)) 

    #calculate reduced thickness at hole
    Lnp_H = hole_pitch_H - hole_length_H
    Lnp_D = hole_pitch_D - hole_length_D

    kg = 0.60
    kd = 0.80

    tg_H = RMI.v2021.eq8_2__1(kg, t, Lnp_H, hole_pitch_H)
    tg_D = RMI.v2021.eq8_2__1(kg, t, Lnp_D, hole_pitch_D)

    td_H = RMI.v2021.eq8_2__2(kd, t, Lnp_H, hole_pitch_H)
    td_D = RMI.v2021.eq8_2__2(kd, t, Lnp_D, hole_pitch_D)

    #define cross-section element thicknesses
    num_elem = length(geometry.x) - 1
    tg = fill(t, num_elem)
    tg[geometry.D_hole_element_index] .= tg_D
    tg[geometry.H_hole_element_index] .= tg_H

    td = fill(t, num_elem)
    td[geometry.D_hole_element_index] .= td_D
    td[geometry.H_hole_element_index] .= td_H

    #net section properties 
    xy_coords_with_holes = [[geometry.x[i], geometry.y[i]] for i in eachindex(geometry.x)]
    net_section_properties = CrossSection.Properties.open_thin_walled(xy_coords_with_holes, tg) 


    #elastic buckling properties 

    #local buckling, compression
    P = 1.0
    Mxx = 0.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []
    neigs = 1
    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 5)
    local_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    eig = 1
    Pcrℓ = minimum(CUFSM.Tools.get_load_factor(local_buckling_P, eig))

    #local buckling, Mxx 
    P = 0.0
    Mxx = 1.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 5)
    local_buckling_Mxx = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrℓ_xx  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Mxx, eig))


    #local buckling, Myy_neg 
    P = 0.0
    Mxx = 0.0
    Myy = 1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 5)
    local_buckling_Myy_neg = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrℓ_yy_neg  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Myy_neg, eig))

    #local buckling, Myy_pos
    P = 0.0
    Mxx = 0.0
    Myy = -1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.25*minimum([H, D]), 1.25*minimum([H,D]), 5)
    local_buckling_Myy_pos = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrℓ_yy_pos  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Myy_pos, eig))

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
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*Lcrd, 1.5*Lcrd, 9)
    distortional_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Pcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_P, eig))


    #distortional buckling, Mxx
    P = 0.0
    Mxx = 1.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*Lcrd, 2.0*Lcrd, 9)
    distortional_buckling_Mxx = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_Mxx, eig))


    #gather up everything 
    properties = CeeLips(section_inputs, geometry, gross_section_properties, net_section_properties, Lnp_H, Lnp_D, tg_H, tg_D, tg, td_H, td_D, td, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Myy_neg, Mcrℓ_yy_neg, local_buckling_Myy_pos, Mcrℓ_yy_pos, distortional_buckling_P, Pcrd, distortional_buckling_Mxx, Mcrd)

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

    geometry = (coordinates = section_geometry, x=x, y=y, D_hole_element_index = D_hole_element_index, H_hole_element_index=H_hole_element_index)

    return geometry

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

    geometry = (coordinates = section_geometry, x=x, y=y, D_hole_element_index = D_hole_element_index, H_hole_element_index=H_hole_element_index)

    return geometry

end



function rectangular_tube(section_inputs)

    @unpack member_type, section_type, H, D, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D = section_inputs

    # input = RectangularTubeInput(H, D, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D)

    geometry = rectangular_tube_geometry(H, D, R, t, dh_H, dh_D, de_H, de_D)

    #gross section properties 
    gross_section_properties = CrossSection.Properties.closed_thin_walled(geometry.coordinates.center, fill(t, length(geometry.x))) 

    #remove NaNs to allow for writing to JSON
    gross_section_properties.xs = -1
    gross_section_properties.ys = -1
    gross_section_properties.Cw = -1
    gross_section_properties.B1 = -1
    gross_section_properties.B2 = -1
    gross_section_properties.wn = [-1, -1]

    #calculate reduced thickness at hole
    Lnp_H = hole_pitch_H - hole_length_H
    Lnp_D = hole_pitch_D - hole_length_D

    kg = 0.60

    tg_H = RMI.v2021.eq8_2__1(kg, t, Lnp_H, hole_pitch_H)
    tg_D = RMI.v2021.eq8_2__1(kg, t, Lnp_D, hole_pitch_D)

    #define cross-section element thicknesses
    num_elem = length(geometry.x)
    tg = fill(t, num_elem)
    tg[geometry.D_hole_element_index] .= tg_D
    tg[geometry.H_hole_element_index] .= tg_H

    #net section properties 
    xy_coords_with_holes = [[geometry.x[i], geometry.y[i]] for i in eachindex(geometry.x)]
    net_section_properties = CrossSection.Properties.closed_thin_walled(xy_coords_with_holes, tg) 
    # net_section_properties = CrossSection.Properties.closed_thin_walled(geometry.coordinates.center, tg) 

    net_section_properties.xs = -1
    net_section_properties.ys = -1
    net_section_properties.Cw = -1
    net_section_properties.B1 = -1
    net_section_properties.B2 = -1
    net_section_properties.wn = [-1, -1]

    #local buckling, compression
    P = 1.0
    Mxx = 0.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []
    neigs = 1
    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 5)
    local_buckling_P = CUFSM.Tools.closed_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs)
    eig = 1
    Pcrℓ = minimum(CUFSM.Tools.get_load_factor(local_buckling_P, eig))

    #local buckling, Mxx 
    P = 0.0
    Mxx = 1.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 5)
    local_buckling_Mxx = CUFSM.Tools.closed_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs)
    Mcrℓ_xx  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Mxx, eig))

    #local buckling, Myy_neg 
    P = 0.0
    Mxx = 0.0
    Myy = 1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 5)
    local_buckling_Myy_neg = CUFSM.Tools.closed_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs)
    Mcrℓ_yy_neg  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Myy_neg, eig))

    #local buckling, Myy_pos
    P = 0.0
    Mxx = 0.0
    Myy = -1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.25*minimum([H, D]), 1.25*minimum([H,D]), 5)
    local_buckling_Myy_pos = CUFSM.Tools.closed_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs)
    Mcrℓ_yy_pos  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Myy_pos, eig))

    #gather up everything 
    properties = RectangularTube(section_inputs, geometry, gross_section_properties, net_section_properties, Lnp_H, Lnp_D, tg_H, tg_D, tg, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Myy_neg, Mcrℓ_yy_neg, local_buckling_Myy_pos, Mcrℓ_yy_pos)

    return properties 

end


function cee_with_lips_rib_geometry(H, D, L, R, t, dh_H, dh_D, de_H, de_D, rib_depth, rib_length, rib_radius_start, rib_radius_peak)

    θ_rib = atan(rib_depth/(rib_length/4))

    segments = [L, D, H/2-rib_length/2 + rib_length/4, rib_depth/sin(θ_rib), rib_depth/sin(θ_rib), H/2-rib_length/2 + rib_length/4, D, L]

    θ = [π/2, π, -π/2, -π/4, -3π/4, -π/2, 0.0, π/2]
    r = [R, R, rib_radius_start, rib_radius_peak, rib_radius_start, R, R]
    n = [4, 4, 3, 4, 4, 3, 4, 4]
    n_r = [3, 3, 3, 3, 3, 3, 3]

    section_geometry = CrossSection.Geometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D, H-L))

    x = [section_geometry.center[i][1] for i in eachindex(section_geometry.center)]
    y = [section_geometry.center[i][2] for i in eachindex(section_geometry.center)]

    # x = [section_geometry.right[i][1] for i in eachindex(section_geometry.center)]
    # y = [section_geometry.right[i][2] for i in eachindex(section_geometry.center)]

    nodes = [x y zeros(Float64, length(x))]

    #flange hole minus
    D_flat_minus = segments[7] - 2*R 
    xloc = t/2 + (R-t/2) + D_flat_minus/n[7] * 1.5
    yloc = t/2
    zloc = 0.0
    atol_x = D_flat_minus/n[7] * 0.55
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

    index_start = n[1] + n_r[1] + n[2] + n_r[2] + 1
    index_end = sum(n) + sum(n_r) - n[end] - n_r[end] - n[end-1] - n_r[end-1] + 1 
    H_flat = y[index_start] - y[index_end]

    index_start = n[1] + n_r[1] + n[2] + n_r[2] + n[3] + 1
    index_end = sum(n) + sum(n_r) - n[end] - n_r[end] - n[end-1] - n_r[end-1] - n[end-2] + 1 

    centerline_rib_length = y[index_start] - y[index_end]

    xloc = t/2
    H_flat_from_rib = H_flat/2 - centerline_rib_length/2
    yloc = H/2 + centerline_rib_length/2 + H_flat_from_rib/n[3] * 1.5
    zloc = 0.0
    atol_x = 0.0
    atol_y = H_flat_from_rib/n[3] * 0.55
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_H_plus = hole_node_index[2] - 1

    y[hole_node_index[1]] = H - de_H + dh_H/2
    y[hole_node_index[2]] = H - de_H - dh_H/2


    #web hole minus
    xloc = t/2
    yloc = H/2 - centerline_rib_length/2 - H_flat_from_rib/n[3] * 1.5
    zloc = 0.0
    atol_x = 0.0
    atol_y = H_flat/n[3] * 0.55
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_H_minus = hole_node_index[2] - 1

    y[hole_node_index[1]] = de_H + dh_H/2
    y[hole_node_index[2]] = de_H - dh_H/2

    H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]

    geometry = (coordinates = section_geometry, x=x, y=y, D_hole_element_index = D_hole_element_index, H_hole_element_index=H_hole_element_index)

    return geometry

end




function cee_with_lips_rib(section_inputs)

    @unpack member_type, section_type, H, D, L, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D, rib_depth, rib_length, rib_radius_start, rib_radius_peak = section_inputs

    # input = CeeLipsRibInput(H, D, L, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D, rib_depth, rib_length, rib_radius_start, rib_radius_peak)

    geometry = cee_with_lips_rib_geometry(H, D, L, R, t, dh_H, dh_D, de_H, de_D, rib_depth, rib_length, rib_radius_start, rib_radius_peak)

    #gross section properties 
    gross_section_properties = CrossSection.Properties.open_thin_walled(geometry.coordinates.center, fill(t, length(geometry.x)-1)) 

    #calculate reduced thickness at hole
    Lnp_H = hole_pitch_H - hole_length_H
    Lnp_D = hole_pitch_D - hole_length_D

    kg = 0.60
    kd = 0.80

    tg_H = RMI.v2021.eq8_2__1(kg, t, Lnp_H, hole_pitch_H)
    tg_D = RMI.v2021.eq8_2__1(kg, t, Lnp_D, hole_pitch_D)

    td_H = RMI.v2021.eq8_2__2(kd, t, Lnp_H, hole_pitch_H)
    td_D = RMI.v2021.eq8_2__2(kd, t, Lnp_D, hole_pitch_D)

    #define cross-section element thicknesses
    num_elem = length(geometry.x) - 1
    tg = fill(t, num_elem)
    tg[geometry.D_hole_element_index] .= tg_D
    tg[geometry.H_hole_element_index] .= tg_H

    td = fill(t, num_elem)
    td[geometry.D_hole_element_index] .= td_D
    td[geometry.H_hole_element_index] .= td_H

    #net section properties 
    xy_coords_with_holes = [[geometry.x[i], geometry.y[i]] for i in eachindex(geometry.x)]
    net_section_properties = CrossSection.Properties.open_thin_walled(xy_coords_with_holes, tg) 
    # net_section_properties = CrossSection.Properties.open_thin_walled(geometry.coordinates.center, tg) 


    #elastic buckling properties 

    #local buckling, compression
    P = 1.0
    Mxx = 0.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []
    neigs = 1
    lengths = range(1.0*minimum([H, D]), 1.75*minimum([H,D]), 5)
    local_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    eig = 1
    Pcrℓ = minimum(CUFSM.Tools.get_load_factor(local_buckling_P, eig))

    #local buckling, Mxx 
    P = 0.0
    Mxx = 1.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 5)
    local_buckling_Mxx = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrℓ_xx  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Mxx, eig))


    #local buckling, Myy_neg 
    P = 0.0
    Mxx = 0.0
    Myy = 1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*minimum([H, D]), 1.5*minimum([H,D]), 5)
    local_buckling_Myy_neg = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrℓ_yy_neg  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Myy_neg, eig))

    #local buckling, Myy_pos
    P = 0.0
    Mxx = 0.0
    Myy = -1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.25*minimum([H, D]), 1.25*minimum([H,D]), 5)
    local_buckling_Myy_pos = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrℓ_yy_pos  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Myy_pos, eig))

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
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*Lcrd, 1.5*Lcrd, 9)
    distortional_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Pcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_P, eig))


    #distortional buckling, Mxx
    P = 0.0
    Mxx = 1.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*Lcrd, 2.0*Lcrd, 9)
    distortional_buckling_Mxx = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_Mxx, eig))


    #gather up everything 
    properties = CeeLipsRib(section_inputs, geometry, gross_section_properties, net_section_properties, Lnp_H, Lnp_D, tg_H, tg_D, tg, td_H, td_D, td, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Myy_neg, Mcrℓ_yy_neg, local_buckling_Myy_pos, Mcrℓ_yy_pos, distortional_buckling_P, Pcrd, distortional_buckling_Mxx, Mcrd)

    return properties 

end




function hat_with_lips_rib_geometry(H, D1, D2, D3, A, X, L, R, t, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2, rib_depth, rib_length, rib_radius_start, rib_radius_peak)

    D = D1 + D2 + D3
    α = A - π/2
    yc = t * tan(α)
    Δ = π - A
    T = R * tan(Δ/2)

    θ_rib = atan(rib_depth/(rib_length/4))

    segments = [L-t, D2-t-yc, (X-t)/sin(π-A), D1, H/2-rib_length/4, rib_depth/sin(θ_rib), rib_depth/sin(θ_rib), H/2-rib_length/4, D1, (X-t)/sin(π-A), D2-t-yc, L-t]

    θ = [-π/2, π, A, π, -π/2, -π/4, -3π/4, -π/2, 0.0, π-A, 0.0, -π/2]
    r = [R-t, R-t, R, R, rib_radius_start, rib_radius_peak, rib_radius_start, R, R, R-t, R-t]
    n = [4, 3, 4, 3, 3, 4, 4, 3, 3, 4, 3, 4]
    n_r = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]

    section_geometry = CrossSection.Geometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D-t, H - (X - L)))

    x = [section_geometry.center[i][1] for i in eachindex(section_geometry.center)]
    y = [section_geometry.center[i][2] for i in eachindex(section_geometry.center)]

    nodes = [x y zeros(Float64, length(x))]

    #D1 flange hole minus
    D1_flat_minus = segments[9] - R - T 
    xloc = t/2 + (R-t/2) + D1_flat_minus/n[9] * 1.5
    yloc = t/2
    zloc = 0.0
    atol_x = D1_flat_minus/n[9] * 0.55
    atol_y = 0.0
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_D1_minus = hole_node_index[2] - 1

    x[hole_node_index[1]] = de_D1 - dh_D1/2
    x[hole_node_index[2]] = de_D1 + dh_D1/2


    #D1 flange hole plus
    D1_flat_plus = segments[4] - R - T 
    xloc = t/2 + (R-t/2) + D1_flat_plus/n[4] * 1.5
    yloc = H-t/2
    zloc = 0.0
    atol_x = D1_flat_plus/n[4] * 0.55
    atol_y = 0.0
    atol_z = 0.0 
    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_D1_plus = hole_node_index[2] - 1

    x[hole_node_index[2]] = de_D1 - dh_D1/2
    x[hole_node_index[1]] = de_D1 + dh_D1/2

    D1_hole_element_index = [hole_element_index_D1_plus; hole_element_index_D1_minus]


    #D2 flange hole minus
    D2_flat_minus = segments[11] - (R-t) - (R-t)*tan(Δ/2) 
    xloc = D - t - (R-t) - D2_flat_minus/2
    yloc = X - t/2
    zloc = 0.0
    atol_x = D2_flat_minus/n[11] * 0.55
    atol_y = 0.001
    atol_z = 0.0 
    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_D2_minus = hole_node_index[2] - 1

    x[hole_node_index[1]] = D - de_D2 - dh_D2/2
    x[hole_node_index[2]] = D - de_D2 + dh_D2/2



    #D2 flange hole plus
    D2_flat_plus = segments[2] - (R-t) - (R-t)*tan(Δ/2) 
    xloc = D - t - (R-t) - D2_flat_plus/2
    yloc = H - (X - t/2)
    zloc = 0.0
    atol_x = D2_flat_plus/n[2] * 0.55
    atol_y = 0.001
    atol_z = 0.0 
    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_D2_plus = hole_node_index[2] - 1

    x[hole_node_index[2]] = D - de_D2 - dh_D2/2
    x[hole_node_index[1]] = D - de_D2 + dh_D2/2

    D2_hole_element_index = [hole_element_index_D2_plus; hole_element_index_D2_minus]


    #web hole plus

    # index_start = n[1] + n_r[1] + n[2] + n_r[2] + n[3] + n_r[3] + n[4] + n_r[4] + 1
    # index_end = sum(n) + sum(n_r) - n[end] - n_r[end] - n[end-1] - n_r[end-1] + 1 
    H_flat = H - 2 * R

    # index_start = n[1] + n_r[1] + n[2] + n_r[2] + n[3] + 1
    index_start = sum(n[1:5]) + sum(n_r[1:4]) + 1
    # index_end = sum(n) + sum(n_r) - n[end] - n_r[end] - n[end-1] - n_r[end-1] - n[end-2] + 1 
    index_end = sum(n[1:7]) + sum(n_r[1:7]) + 1

    centerline_rib_length = y[index_start] - y[index_end]

    xloc = t/2
    H_flat_from_rib = H_flat/2 - centerline_rib_length/2
    yloc = H/2 + centerline_rib_length/2 + H_flat_from_rib/n[5] * 1.5
    zloc = 0.0
    atol_x = 0.001
    atol_y = H_flat_from_rib/n[5] * 0.55
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_H_plus = hole_node_index[2] - 1

    y[hole_node_index[1]] = H - de_H + dh_H/2
    y[hole_node_index[2]] = H - de_H - dh_H/2



    #web hole minus

    # index_start = n[1] + n_r[1] + n[2] + n_r[2] + n[3] + n_r[3] + n[4] + n_r[4] + 1
    # index_end = sum(n) + sum(n_r) - n[end] - n_r[end] - n[end-1] - n_r[end-1] + 1 
    H_flat = H - 2 * R

    # index_start = n[1] + n_r[1] + n[2] + n_r[2] + n[3] + 1
    index_start = sum(n[1:5]) + sum(n_r[1:4]) + 1
    # index_end = sum(n) + sum(n_r) - n[end] - n_r[end] - n[end-1] - n_r[end-1] - n[end-2] + 1 
    index_end = sum(n[1:7]) + sum(n_r[1:7]) + 1

    centerline_rib_length = y[index_start] - y[index_end]

    xloc = t/2
    H_flat_from_rib = H_flat/2 - centerline_rib_length/2
    yloc = H/2 - centerline_rib_length/2 - H_flat_from_rib/n[8] * 1.5
    zloc = 0.0
    atol_x = 0.001
    atol_y = H_flat_from_rib/n[8] * 0.55
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_H_minus = hole_node_index[2] - 1

    H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]


    y[hole_node_index[1]] = de_H + dh_H/2
    y[hole_node_index[2]] = de_H - dh_H/2

    geometry = (coordinates = section_geometry, x=x, y=y, D1_hole_element_index = D1_hole_element_index, D2_hole_element_index = D2_hole_element_index, H_hole_element_index=H_hole_element_index)

    return geometry

end




function hat_with_lips_rib(section_inputs)

    @unpack member_type, section_type, H, D1, D2, D3, A, X, L, R, t, E, ν, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2, hole_pitch_H, hole_pitch_D1, hole_pitch_D2, hole_length_H, hole_length_D1, hole_length_D2, rib_depth, rib_length, rib_radius_start, rib_radius_peak = section_inputs

    D = D1 + D2 + D3

    # input = HatLipsRibInput(H, D1, D2, D3, A, X, L, R, t, E, ν, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2, hole_pitch_H, hole_pitch_D1, hole_pitch_D2, hole_length_H, hole_length_D1, hole_length_D2, rib_depth, rib_length, rib_radius_start, rib_radius_peak)

    geometry = hat_with_lips_rib_geometry(H, D1, D2, D3, A, X, L, R, t, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2, rib_depth, rib_length, rib_radius_start, rib_radius_peak)
    #gross section properties 
    gross_section_properties = CrossSection.Properties.open_thin_walled(geometry.coordinates.center, fill(t, length(geometry.x)-1)) 

    #calculate reduced thickness at hole
    Lnp_H = hole_pitch_H - hole_length_H
    Lnp_D1 = hole_pitch_D1 - hole_length_D1
    Lnp_D2 = hole_pitch_D2 - hole_length_D2

    kg = 0.60
    kd = 0.80

    tg_H = RMI.v2021.eq8_2__1(kg, t, Lnp_H, hole_pitch_H)
    tg_D1 = RMI.v2021.eq8_2__1(kg, t, Lnp_D1, hole_pitch_D1)
    tg_D2 = RMI.v2021.eq8_2__1(kg, t, Lnp_D2, hole_pitch_D2)

    td_H = RMI.v2021.eq8_2__2(kd, t, Lnp_H, hole_pitch_H)
    td_D1 = RMI.v2021.eq8_2__2(kd, t, Lnp_D1, hole_pitch_D1)
    td_D2 = RMI.v2021.eq8_2__2(kd, t, Lnp_D2, hole_pitch_D2)

    #define cross-section element thicknesses
    num_elem = length(geometry.x) - 1
    tg = fill(t, num_elem)
    tg[geometry.D1_hole_element_index] .= tg_D1
    tg[geometry.D2_hole_element_index] .= tg_D2
    tg[geometry.H_hole_element_index] .= tg_H

    td = fill(t, num_elem)
    td[geometry.D1_hole_element_index] .= td_D1
    td[geometry.D2_hole_element_index] .= td_D2
    td[geometry.H_hole_element_index] .= td_H

    #net section properties 
    xy_coords_with_holes = [[geometry.x[i], geometry.y[i]] for i in eachindex(geometry.x)]
    net_section_properties = CrossSection.Properties.open_thin_walled(xy_coords_with_holes, tg) 
    # net_section_properties = CrossSection.Properties.open_thin_walled(geometry.coordinates.center, tg) 


    #elastic buckling properties 

    #local buckling, compression
    P = 1.0
    Mxx = 0.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []
    neigs = 1
    lengths = range(1.0*minimum([H, D]), 1.75*minimum([H,D]), 5)
    local_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    eig = 1
    Pcrℓ = minimum(CUFSM.Tools.get_load_factor(local_buckling_P, eig))

    #local buckling, Mxx 
    P = 0.0
    Mxx = 1.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.25*minimum([H, D]), 0.6*minimum([H,D]), 5)
    local_buckling_Mxx = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrℓ_xx  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Mxx, eig))


    #local buckling, Myy_neg 
    P = 0.0
    Mxx = 0.0
    Myy = 1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(1.0*minimum([H, D]), 2.0*minimum([H,D]), 5)
    local_buckling_Myy_neg = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrℓ_yy_neg  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Myy_neg, eig))

    #distortional buckling, Myy_pos
    P = 0.0
    Mxx = 0.0
    Myy = -1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(4.0*minimum([H, D]), 8.0*minimum([H,D]), 5)
    distortional_buckling_Myy_pos = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrd_yy_pos  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_Myy_pos, eig))

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
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*Lcrd, 2.0*Lcrd, 9)
    distortional_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Pcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_P, eig))


    #distortional buckling, Mxx
    P = 0.0
    Mxx = 1.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*Lcrd, 2.0*Lcrd, 9)
    distortional_buckling_Mxx = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_Mxx, eig))


    #gather up everything 
    properties = HatLipsRib(section_inputs, geometry, gross_section_properties, net_section_properties, Lnp_H, Lnp_D1, Lnp_D2, tg_H, tg_D1, tg_D2, tg, td_H, td_D1, td_D2, td, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Myy_neg, Mcrℓ_yy_neg, distortional_buckling_Myy_pos, Mcrd_yy_pos, distortional_buckling_P, Pcrd, distortional_buckling_Mxx, Mcrd)

    return properties 

end


function hat_with_rib_geometry(H, D1, D2, D3, A, X, R, t, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2, rib_depth, rib_length, rib_radius_start, rib_radius_peak)

    D = D1 + D2 + D3
    α = A - π/2
    yc = t * tan(α)
    Δ = π - A
    T = R * tan(Δ/2)

    θ_rib = atan(rib_depth/(rib_length/4))

    segments = [D2-yc, (X-t)/sin(π-A), D1, H/2-rib_length/4, rib_depth/sin(θ_rib), rib_depth/sin(θ_rib), H/2-rib_length/4, D1, (X-t)/sin(π-A), D2-yc]

    θ = [π, A, π, -π/2, -π/4, -3π/4, -π/2, 0.0, π-A, 0.0]
    r = [R-t, R, R, rib_radius_start, rib_radius_peak, rib_radius_start, R, R, R-t]
    n = [3, 4, 3, 3, 4, 4, 3, 3, 4, 3]
    n_r = [3, 3, 3, 3, 3, 3, 3, 3, 3]

    section_geometry = CrossSection.Geometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D, H - (X-t)))

    x = [section_geometry.center[i][1] for i in eachindex(section_geometry.center)]
    y = [section_geometry.center[i][2] for i in eachindex(section_geometry.center)]

    nodes = [x y zeros(Float64, length(x))]

    #D1 flange hole minus
    D1_flat_minus = segments[8] - R - T 
    xloc = t/2 + (R-t/2) + D1_flat_minus/n[8] * 1.5
    yloc = t/2
    zloc = 0.0
    atol_x = D1_flat_minus/n[8] * 0.55
    atol_y = 0.0
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_D1_minus = hole_node_index[2] - 1

    x[hole_node_index[1]] = de_D1 - dh_D1/2
    x[hole_node_index[2]] = de_D1 + dh_D1/2


    #D1 flange hole plus
    D1_flat_plus = segments[3] - R - T 
    xloc = t/2 + (R-t/2) + D1_flat_plus/n[3] * 1.5
    yloc = H-t/2
    zloc = 0.0
    atol_x = D1_flat_plus/n[3] * 0.55
    atol_y = 0.0
    atol_z = 0.0 
    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_D1_plus = hole_node_index[2] - 1

    x[hole_node_index[2]] = de_D1 - dh_D1/2
    x[hole_node_index[1]] = de_D1 + dh_D1/2

    D1_hole_element_index = [hole_element_index_D1_plus; hole_element_index_D1_minus]


    #D2 flange hole minus
    D2_flat_minus = segments[10] - (R-t)*tan(Δ/2) 
    xloc = D - D2_flat_minus/2
    yloc = X - t/2
    zloc = 0.0
    atol_x = D2_flat_minus/n[10] * 0.55
    atol_y = 0.001
    atol_z = 0.0 
    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_D2_minus = hole_node_index[2] - 1

    x[hole_node_index[1]] = D - de_D2 - dh_D2/2
    x[hole_node_index[2]] = D - de_D2 + dh_D2/2



    #D2 flange hole plus
    D2_flat_plus = segments[1] - (R-t)*tan(Δ/2) 
    xloc = D - D2_flat_plus/2
    yloc = H - (X - t/2)
    zloc = 0.0
    atol_x = D2_flat_plus/n[1] * 0.55
    atol_y = 0.001
    atol_z = 0.0 
    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_D2_plus = hole_node_index[2] - 1

    x[hole_node_index[2]] = D - de_D2 - dh_D2/2
    x[hole_node_index[1]] = D - de_D2 + dh_D2/2

    D2_hole_element_index = [hole_element_index_D2_plus; hole_element_index_D2_minus]


    #web hole plus

    # index_start = n[1] + n_r[1] + n[2] + n_r[2] + n[3] + n_r[3] + n[4] + n_r[4] + 1
    # index_end = sum(n) + sum(n_r) - n[end] - n_r[end] - n[end-1] - n_r[end-1] + 1 
    H_flat = H - 2 * R

    # index_start = n[1] + n_r[1] + n[2] + n_r[2] + n[3] + 1
    index_start = sum(n[1:4]) + sum(n_r[1:3]) + 1
    # index_end = sum(n) + sum(n_r) - n[end] - n_r[end] - n[end-1] - n_r[end-1] - n[end-2] + 1 
    index_end = sum(n[1:6]) + sum(n_r[1:6]) + 1

    centerline_rib_length = y[index_start] - y[index_end]

    xloc = t/2
    H_flat_from_rib = H_flat/2 - centerline_rib_length/2
    yloc = H/2 + centerline_rib_length/2 + H_flat_from_rib/n[4] * 1.5
    zloc = 0.0
    atol_x = 0.001
    atol_y = H_flat_from_rib/n[4] * 0.55
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_H_plus = hole_node_index[2] - 1

    y[hole_node_index[1]] = H - de_H + dh_H/2
    y[hole_node_index[2]] = H - de_H - dh_H/2



    #web hole minus

    # index_start = n[1] + n_r[1] + n[2] + n_r[2] + n[3] + n_r[3] + n[4] + n_r[4] + 1
    # index_end = sum(n) + sum(n_r) - n[end] - n_r[end] - n[end-1] - n_r[end-1] + 1 
    H_flat = H - 2 * R

    # index_start = n[1] + n_r[1] + n[2] + n_r[2] + n[3] + 1
    index_start = sum(n[1:4]) + sum(n_r[1:3]) + 1
    # index_end = sum(n) + sum(n_r) - n[end] - n_r[end] - n[end-1] - n_r[end-1] - n[end-2] + 1 
    index_end = sum(n[1:6]) + sum(n_r[1:6]) + 1

    centerline_rib_length = y[index_start] - y[index_end]

    xloc = t/2
    H_flat_from_rib = H_flat/2 - centerline_rib_length/2
    yloc = H/2 - centerline_rib_length/2 - H_flat_from_rib/n[7] * 1.5
    zloc = 0.0
    atol_x = 0.001
    atol_y = H_flat_from_rib/n[7] * 0.55
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_H_minus = hole_node_index[2] - 1

    H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]


    y[hole_node_index[1]] = de_H + dh_H/2
    y[hole_node_index[2]] = de_H - dh_H/2

    geometry = (coordinates = section_geometry, x=x, y=y, D1_hole_element_index = D1_hole_element_index, D2_hole_element_index = D2_hole_element_index, H_hole_element_index=H_hole_element_index)

    return geometry

end



function hat_with_rib(section_inputs)


    @unpack member_type, section_type, H, D1, D2, D3, A, X, R, t, E, ν, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2, hole_pitch_H, hole_pitch_D1, hole_pitch_D2, hole_length_H, hole_length_D1, hole_length_D2, rib_depth, rib_length, rib_radius_start, rib_radius_peak = section_inputs

    D = D1 + D2 + D3

    # input = HatRibInput(H, D1, D2, D3, A, X, R, t, E, ν, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2, hole_pitch_H, hole_pitch_D1, hole_pitch_D2, hole_length_H, hole_length_D1, hole_length_D2, rib_depth, rib_length, rib_radius_start, rib_radius_peak)


    geometry = hat_with_rib_geometry(H, D1, D2, D3, A, X, R, t, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2, rib_depth, rib_length, rib_radius_start, rib_radius_peak)
    #gross section properties 
    gross_section_properties = CrossSection.Properties.open_thin_walled(geometry.coordinates.center, fill(t, length(geometry.x)-1)) 

    #calculate reduced thickness at hole
    Lnp_H = hole_pitch_H - hole_length_H
    Lnp_D1 = hole_pitch_D1 - hole_length_D1
    Lnp_D2 = hole_pitch_D2 - hole_length_D2

    kg = 0.60
    kd = 0.80

    tg_H = RMI.v2021.eq8_2__1(kg, t, Lnp_H, hole_pitch_H)
    tg_D1 = RMI.v2021.eq8_2__1(kg, t, Lnp_D1, hole_pitch_D1)
    tg_D2 = RMI.v2021.eq8_2__1(kg, t, Lnp_D2, hole_pitch_D2)

    td_H = RMI.v2021.eq8_2__2(kd, t, Lnp_H, hole_pitch_H)
    td_D1 = RMI.v2021.eq8_2__2(kd, t, Lnp_D1, hole_pitch_D1)
    td_D2 = RMI.v2021.eq8_2__2(kd, t, Lnp_D2, hole_pitch_D2)

    #define cross-section element thicknesses
    num_elem = length(geometry.x) - 1
    tg = fill(t, num_elem)
    tg[geometry.D1_hole_element_index] .= tg_D1
    tg[geometry.D2_hole_element_index] .= tg_D2
    tg[geometry.H_hole_element_index] .= tg_H

    td = fill(t, num_elem)
    td[geometry.D1_hole_element_index] .= td_D1
    td[geometry.D2_hole_element_index] .= td_D2
    td[geometry.H_hole_element_index] .= td_H

    #net section properties 
    xy_coords_with_holes = [[geometry.x[i], geometry.y[i]] for i in eachindex(geometry.x)]
    net_section_properties = CrossSection.Properties.open_thin_walled(xy_coords_with_holes, tg) 
    # net_section_properties = CrossSection.Properties.open_thin_walled(geometry.coordinates.center, tg) 


    #elastic buckling properties 

    #local buckling, compression
    P = 1.0
    Mxx = 0.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []
    neigs = 1
    lengths = range(0.75*maximum([H, D]), 1.75*maximum([H,D]), 5)
    local_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    eig = 1
    Pcrℓ = minimum(CUFSM.Tools.get_load_factor(local_buckling_P, eig))

    #local buckling, Mxx 
    P = 0.0
    Mxx = 1.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 5)
    local_buckling_Mxx = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrℓ_xx  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Mxx, eig))


    #local buckling, Myy_neg 
    P = 0.0
    Mxx = 0.0
    Myy = 1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(1.0*minimum([H, D]), 2.0*minimum([H,D]), 5)
    local_buckling_Myy_neg = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrℓ_yy_neg  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Myy_neg, eig))

    #distortional buckling, Myy_pos
    P = 0.0
    Mxx = 0.0
    Myy = -1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(3.0*minimum([H, D]), 7.0*minimum([H,D]), 5)
    distortional_buckling_Myy_pos = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrd_yy_pos  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_Myy_pos, eig))

    # #distortional buckling 

    #find approximate distortional buckling half-wavelength
    CorZ = 0
    b = D
    d = 0.0  #try to use equation without lips 
    θ = 90.0
    Af,Jf,Ixf,Iyf,Ixyf,Cwf,xof,hxf,hyf,yof = AISIS100.v16.table23131(CorZ,t,b,d,θ)

    ho = H
    μ = ν
    Lm = 999999999.0
    Lcrd, not_used = AISIS100.v16.app23334(ho, μ, t, Ixf, xof, hxf, Cwf, Ixyf, Iyf, Lm)


    #distortional buckling, compression 
    P = 1.0
    Mxx = 0.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(1.5*Lcrd, 2.5*Lcrd, 9)
    distortional_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Pcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_P, eig))


    #distortional buckling, Mxx
    P = 0.0
    Mxx = 1.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(1.25*Lcrd, 2.75*Lcrd, 9)
    distortional_buckling_Mxx = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_Mxx, eig))


    #gather up everything 
    properties = HatRib(section_inputs, geometry, gross_section_properties, net_section_properties, Lnp_H, Lnp_D1, Lnp_D2, tg_H, tg_D1, tg_D2, tg, td_H, td_D1, td_D2, td, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Myy_neg, Mcrℓ_yy_neg, distortional_buckling_Myy_pos, Mcrd_yy_pos, distortional_buckling_P, Pcrd, distortional_buckling_Mxx, Mcrd)

    return properties 

end


function hat_with_lips_trapezoidal_rib_geometry(H, D1, D2, D3, A1, X, L, R, t, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2, A2, hr, wr, Rr)

    D = D1 + D2 + D3
    α = A1 - π/2
    yc = t * tan(α)
    Δ = π - A1
    T = R * tan(Δ/2)

    α2 = A2 - π/2
    yc2 = t * tan(α2)
    Δ2 = π - A2
    T2= Rr * tan(Δ2/2)
    
    hr1 = (wr-t) / tan(Δ2) - yc2
    
    
    segments = [L-t, D2-t-yc, (X-t)/sin(π-A1), D1, H/2-hr/2-hr1, (wr-t)/sin(Δ2), hr - 2*yc2, (wr-t)/sin(Δ2), H/2-hr/2 - hr1, D1, (X-t)/sin(π-A1), D2-t-yc, L-t]
    
    θ = [-π/2, π, A1, π, -π/2, -π/2 + Δ2, -π/2, -π/2 - Δ2, -π/2, 0.0, π-A1, 0.0, -π/2]
    r = [R-t, R-t, R, R, Rr, Rr-t, Rr-t, Rr, R, R, R-t, R-t]
    n = [4, 3, 4, 3, 3, 4, 4, 4, 3, 3, 4, 3, 4]
    n_r = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]

    section_geometry = CrossSection.Geometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D-t, H - (X-L)))

    x = [section_geometry.center[i][1] for i in eachindex(section_geometry.center)]
    y = [section_geometry.center[i][2] for i in eachindex(section_geometry.center)]

    nodes = [x y zeros(Float64, length(x))]

    #D1 flange hole minus
    D1_flat_minus = segments[10] - R - T 
    xloc = t/2 + (R-t/2) + D1_flat_minus/n[10] * 1.5
    yloc = t/2
    zloc = 0.0
    atol_x = D1_flat_minus/n[10] * 0.55
    atol_y = 0.0
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_D1_minus = hole_node_index[2] - 1

    x[hole_node_index[1]] = de_D1 - dh_D1/2
    x[hole_node_index[2]] = de_D1 + dh_D1/2


    #D1 flange hole plus
    D1_flat_plus = segments[4] - R - T 
    xloc = t/2 + (R-t/2) + D1_flat_plus/n[4] * 1.5
    yloc = H-t/2
    zloc = 0.0
    atol_x = D1_flat_plus/n[4] * 0.55
    atol_y = 0.0
    atol_z = 0.0 
    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_D1_plus = hole_node_index[2] - 1

    x[hole_node_index[2]] = de_D1 - dh_D1/2
    x[hole_node_index[1]] = de_D1 + dh_D1/2

    D1_hole_element_index = [hole_element_index_D1_plus; hole_element_index_D1_minus]


    #D2 flange hole minus
    D2_flat_minus = segments[12] - (R-t) - (R-t)*tan(Δ/2) 
    xloc = D - t - (R-t) - D2_flat_minus/2
    yloc = X - t/2
    zloc = 0.0
    atol_x = D2_flat_minus/n[12] * 0.55
    atol_y = 0.001
    atol_z = 0.0 
    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_D2_minus = hole_node_index[2] - 1

    x[hole_node_index[1]] = D - de_D2 - dh_D2/2
    x[hole_node_index[2]] = D - de_D2 + dh_D2/2



    #D2 flange hole plus
    D2_flat_plus = segments[2] - (R-t) - (R-t)*tan(Δ/2) 
    xloc = D - t - (R-t) - D2_flat_plus/2
    yloc = H - (X - t/2)
    zloc = 0.0
    atol_x = D2_flat_plus/n[2] * 0.55
    atol_y = 0.001
    atol_z = 0.0 
    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_D2_plus = hole_node_index[2] - 1

    x[hole_node_index[2]] = D - de_D2 - dh_D2/2
    x[hole_node_index[1]] = D - de_D2 + dh_D2/2

    D2_hole_element_index = [hole_element_index_D2_plus; hole_element_index_D2_minus]


    #web hole plus

    # index_start = n[1] + n_r[1] + n[2] + n_r[2] + n[3] + n_r[3] + n[4] + n_r[4] + 1
    # index_end = sum(n) + sum(n_r) - n[end] - n_r[end] - n[end-1] - n_r[end-1] + 1 
    H_flat = H - 2 * R

    # index_start = n[1] + n_r[1] + n[2] + n_r[2] + n[3] + 1
    index_start = sum(n[1:5]) + sum(n_r[1:4]) + 1
    # index_end = sum(n) + sum(n_r) - n[end] - n_r[end] - n[end-1] - n_r[end-1] - n[end-2] + 1 
    index_end = sum(n[1:8]) + sum(n_r[1:8]) + 1

    centerline_rib_length = y[index_start] - y[index_end]

    xloc = t/2
    H_flat_from_rib = H_flat/2 - centerline_rib_length/2
    yloc = H/2 + centerline_rib_length/2 + H_flat_from_rib/n[5] * 1.5
    zloc = 0.0
    atol_x = 0.001
    atol_y = H_flat_from_rib/n[5] * 0.55
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_H_plus = hole_node_index[2] - 1

    y[hole_node_index[1]] = H - de_H + dh_H/2
    y[hole_node_index[2]] = H - de_H - dh_H/2



    #web hole minus

    # index_start = n[1] + n_r[1] + n[2] + n_r[2] + n[3] + n_r[3] + n[4] + n_r[4] + 1
    # index_end = sum(n) + sum(n_r) - n[end] - n_r[end] - n[end-1] - n_r[end-1] + 1 
    H_flat = H - 2 * R

    # index_start = n[1] + n_r[1] + n[2] + n_r[2] + n[3] + 1
    # index_start = sum(n[1:5]) + sum(n_r[1:4]) + 1
    # # index_end = sum(n) + sum(n_r) - n[end] - n_r[end] - n[end-1] - n_r[end-1] - n[end-2] + 1 
    # index_end = sum(n[1:7]) + sum(n_r[1:7]) + 1

    # centerline_rib_length = y[index_start] - y[index_end]

    xloc = t/2
    H_flat_from_rib = H_flat/2 - centerline_rib_length/2
    yloc = H/2 - centerline_rib_length/2 - H_flat_from_rib/n[9] * 1.5
    zloc = 0.0
    atol_x = 0.001
    atol_y = H_flat_from_rib/n[9] * 0.55
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_H_minus = hole_node_index[2] - 1

    H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]


    y[hole_node_index[1]] = de_H + dh_H/2
    y[hole_node_index[2]] = de_H - dh_H/2

    geometry = (coordinates = section_geometry, x=x, y=y, D1_hole_element_index = D1_hole_element_index, D2_hole_element_index = D2_hole_element_index, H_hole_element_index=H_hole_element_index)

    return geometry

end



function hat_with_lips_trapezoidal_rib(section_inputs)


    @unpack member_type, section_type, H, D1, D2, D3, A1, X, L, R, t, E, ν, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2, hole_pitch_H, hole_pitch_D1, hole_pitch_D2, hole_length_H, hole_length_D1, hole_length_D2, A2, hr, wr, Rr = section_inputs

    D = D1 + D2 + D3

    # input = HatLipsTrapezoidalRibInput(H, D1, D2, D3, A1, X, L, R, t, E, ν, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2, hole_pitch_H, hole_pitch_D1, hole_pitch_D2, hole_length_H, hole_length_D1, hole_length_D2, A2, hr, wr, Rr)

    geometry = hat_with_lips_trapezoidal_rib_geometry(H, D1, D2, D3, A1, X, L, R, t, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2, A2, hr, wr, Rr)
    #gross section properties 
    gross_section_properties = CrossSection.Properties.open_thin_walled(geometry.coordinates.center, fill(t, length(geometry.x)-1)) 

    #calculate reduced thickness at hole
    Lnp_H = hole_pitch_H - hole_length_H
    Lnp_D1 = hole_pitch_D1 - hole_length_D1
    Lnp_D2 = hole_pitch_D2 - hole_length_D2

    kg = 0.60
    kd = 0.80

    tg_H = RMI.v2021.eq8_2__1(kg, t, Lnp_H, hole_pitch_H)
    tg_D1 = RMI.v2021.eq8_2__1(kg, t, Lnp_D1, hole_pitch_D1)
    tg_D2 = RMI.v2021.eq8_2__1(kg, t, Lnp_D2, hole_pitch_D2)

    td_H = RMI.v2021.eq8_2__2(kd, t, Lnp_H, hole_pitch_H)
    td_D1 = RMI.v2021.eq8_2__2(kd, t, Lnp_D1, hole_pitch_D1)
    td_D2 = RMI.v2021.eq8_2__2(kd, t, Lnp_D2, hole_pitch_D2)

    #define cross-section element thicknesses
    num_elem = length(geometry.x) - 1
    tg = fill(t, num_elem)
    tg[geometry.D1_hole_element_index] .= tg_D1
    tg[geometry.D2_hole_element_index] .= tg_D2
    tg[geometry.H_hole_element_index] .= tg_H

    td = fill(t, num_elem)
    td[geometry.D1_hole_element_index] .= td_D1
    td[geometry.D2_hole_element_index] .= td_D2
    td[geometry.H_hole_element_index] .= td_H

    #net section properties 
    xy_coords_with_holes = [[geometry.x[i], geometry.y[i]] for i in eachindex(geometry.x)]
    net_section_properties = CrossSection.Properties.open_thin_walled(xy_coords_with_holes, tg) 
    # net_section_properties = CrossSection.Properties.open_thin_walled(geometry.coordinates.center, tg) 


    #elastic buckling properties 

    #local buckling, compression
    P = 1.0
    Mxx = 0.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []
    neigs = 1
    lengths = range(0.5*maximum([H, D]), 1.8*maximum([H,D]), 5)
    local_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    eig = 1
    Pcrℓ = minimum(CUFSM.Tools.get_load_factor(local_buckling_P, eig))

    #local buckling, Mxx 
    P = 0.0
    Mxx = 1.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.25*minimum([H, D]), 0.6*minimum([H,D]), 5)
    local_buckling_Mxx = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrℓ_xx  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Mxx, eig))


    #local buckling, Myy_neg 
    P = 0.0
    Mxx = 0.0
    Myy = 1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(1.0*maximum([H, D]), 2.25*maximum([H,D]), 5)
    local_buckling_Myy_neg = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrℓ_yy_neg  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Myy_neg, eig))

    #distortional buckling, Myy_pos
    P = 0.0
    Mxx = 0.0
    Myy = -1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(4.0*minimum([H, D]), 8.0*minimum([H,D]), 5)
    distortional_buckling_Myy_pos = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrd_yy_pos  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_Myy_pos, eig))

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
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*Lcrd, 2.0*Lcrd, 9)
    distortional_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Pcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_P, eig))


    #distortional buckling, Mxx
    P = 0.0
    Mxx = 1.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*Lcrd, 2.0*Lcrd, 9)
    distortional_buckling_Mxx = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_Mxx, eig))


    #gather up everything 
    properties = HatLipsTrapezoidalRib(section_inputs, geometry, gross_section_properties, net_section_properties, Lnp_H, Lnp_D1, Lnp_D2, tg_H, tg_D1, tg_D2, tg, td_H, td_D1, td_D2, td, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Myy_neg, Mcrℓ_yy_neg, distortional_buckling_Myy_pos, Mcrd_yy_pos, distortional_buckling_P, Pcrd, distortional_buckling_Mxx, Mcrd)

    return properties 

end



function unistrut_in_geometry(H, D, L1, L2, R, t, dh_H, dh_D, de_H, de_D, rib_depth, rib_length, rib_radius_start, rib_radius_peak)

    θ_rib = atan(rib_depth/(rib_length/4))

    segments = [L2, L1, D, H/2-rib_length/2 + rib_length/4, rib_depth/sin(θ_rib), rib_depth/sin(θ_rib), H/2-rib_length/2 + rib_length/4, D, L1, L2]

    θ = [0.0, π/2, π, -π/2, -π/4, -3π/4, -π/2, 0.0, π/2, π]
    r = [R, R, R, rib_radius_start, rib_radius_peak, rib_radius_start, R, R, R]
    n = [4, 4, 4, 3, 4, 4, 3, 4, 4, 4]
    n_r = [3, 3, 3, 3, 3, 3, 3, 3, 3]

    section_geometry = CrossSection.Geometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D-L2, H-L1))

    x = [section_geometry.center[i][1] for i in eachindex(section_geometry.center)]
    y = [section_geometry.center[i][2] for i in eachindex(section_geometry.center)]

    nodes = [x y zeros(Float64, length(x))]

    #flange hole minus
    D_flat_minus = segments[8] - 2*R 
    xloc = t/2 + (R-t/2) + D_flat_minus/n[8] * 1.5
    yloc = t/2
    zloc = 0.0
    atol_x = D_flat_minus/n[8] * 0.55
    atol_y = 0.0
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_D_minus = hole_node_index[2] - 1

    x[hole_node_index[1]] = de_D - dh_D/2
    x[hole_node_index[2]] = de_D + dh_D/2

    #flange hole plus
    D_flat_plus = segments[3] - 2*R 
    xloc = t/2 + (R-t/2) + D_flat_plus/n[3] * 1.5
    yloc = H - t/2
    zloc = 0.0
    atol_x = D_flat_plus/n[3] * 0.55
    atol_y = 0.0
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_D_plus = hole_node_index[2] - 1

    x[hole_node_index[1]] = de_D + dh_D/2
    x[hole_node_index[2]] = de_D - dh_D/2

    D_hole_element_index = [hole_element_index_D_plus; hole_element_index_D_minus]

    #web hole plus

    # index_start = n[1] + n_r[1] + n[2] + n_r[2] + n[3] + n_r[3] + 1
    # index_end = sum(n) + sum(n_r) - n[end] - n_r[end] - n[end-1] - n_r[end-1] + 1 
    # H_flat = y[index_start] - y[index_end]
    H_flat = H - 2*R

    index_start = n[1] + n_r[1] + n[2] + n_r[2] + n[3] + n_r[3] + n[4] + 1
    index_end = sum(n) + sum(n_r) - n[end] - n_r[end] - n[end-1] - n_r[end-1] - n[end-2] - n_r[end-2]- n[end-3] + 1 

    centerline_rib_length = y[index_start] - y[index_end]

    xloc = t/2
    H_flat_from_rib = H_flat/2 - centerline_rib_length/2
    yloc = H/2 + centerline_rib_length/2 + H_flat_from_rib/n[4] * 1.5
    zloc = 0.0
    atol_x = 0.0
    atol_y = H_flat_from_rib/n[4] * 0.55
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_H_plus = hole_node_index[2] - 1

    y[hole_node_index[1]] = H - de_H + dh_H/2
    y[hole_node_index[2]] = H - de_H - dh_H/2


    #web hole minus
    xloc = t/2
    yloc = H/2 - centerline_rib_length/2 - H_flat_from_rib/n[4] * 1.5
    zloc = 0.0
    atol_x = 0.0
    atol_y = H_flat/n[4] * 0.55
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_H_minus = hole_node_index[2] - 1

    y[hole_node_index[1]] = de_H + dh_H/2
    y[hole_node_index[2]] = de_H - dh_H/2

    H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]

    geometry = (coordinates = section_geometry, x=x, y=y, D_hole_element_index = D_hole_element_index, H_hole_element_index=H_hole_element_index)

    return geometry

end




function unistrut_in(section_inputs)

    # input = UniStrutInput(H, D, L1, L2, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D, rib_depth, rib_length, rib_radius_start, rib_radius_peak)

    @unpack member_type, section_type, H, D, L1, L2, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D, rib_depth, rib_length, rib_radius_start, rib_radius_peak = section_inputs

    geometry = unistrut_in_geometry(H, D, L1, L2, R, t, dh_H, dh_D, de_H, de_D, rib_depth, rib_length, rib_radius_start, rib_radius_peak)

    #gross section properties 
    gross_section_properties = CrossSection.Properties.open_thin_walled(geometry.coordinates.center, fill(t, length(geometry.x)-1)) 

    #calculate reduced thickness at hole
    Lnp_H = hole_pitch_H - hole_length_H
    Lnp_D = hole_pitch_D - hole_length_D

    kg = 0.60
    kd = 0.80

    tg_H = RMI.v2021.eq8_2__1(kg, t, Lnp_H, hole_pitch_H)
    tg_D = RMI.v2021.eq8_2__1(kg, t, Lnp_D, hole_pitch_D)

    td_H = RMI.v2021.eq8_2__2(kd, t, Lnp_H, hole_pitch_H)
    td_D = RMI.v2021.eq8_2__2(kd, t, Lnp_D, hole_pitch_D)

    #define cross-section element thicknesses
    num_elem = length(geometry.x) - 1
    tg = fill(t, num_elem)
    tg[geometry.D_hole_element_index] .= tg_D
    tg[geometry.H_hole_element_index] .= tg_H

    td = fill(t, num_elem)
    td[geometry.D_hole_element_index] .= td_D
    td[geometry.H_hole_element_index] .= td_H

    #net section properties 
    xy_coords_with_holes = [[geometry.x[i], geometry.y[i]] for i in eachindex(geometry.x)]
    net_section_properties = CrossSection.Properties.open_thin_walled(xy_coords_with_holes, tg) 
    # net_section_properties = CrossSection.Properties.open_thin_walled(geometry.coordinates.center, tg) 


    #elastic buckling properties 

    #local buckling, compression
    P = 1.0
    Mxx = 0.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []
    neigs = 1
    lengths = range(1.0*minimum([H, D]), 1.75*minimum([H,D]), 5)
    local_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    eig = 1
    Pcrℓ = minimum(CUFSM.Tools.get_load_factor(local_buckling_P, eig))

    #local buckling, Mxx 
    P = 0.0
    Mxx = 1.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 5)
    local_buckling_Mxx = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrℓ_xx  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Mxx, eig))


    #local buckling, Myy_neg 
    P = 0.0
    Mxx = 0.0
    Myy = 1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*minimum([H, D]), 1.5*minimum([H,D]), 5)
    local_buckling_Myy_neg = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrℓ_yy_neg  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Myy_neg, eig))

    #local buckling, Myy_pos
    P = 0.0
    Mxx = 0.0
    Myy = -1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.25 * D, 0.75 * D, 5)
    local_buckling_Myy_pos = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrℓ_yy_pos  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Myy_pos, eig))

    #distortional buckling 

    #find approximate distortional buckling half-wavelength
    CorZ = 0
    b = D
    d = L1
    θ = 90.0
    Af,Jf,Ixf,Iyf,Ixyf,Cwf,xof,hxf,hyf,yof = AISIS100.v16.table23131(CorZ,t,b,d,θ)

    ho = H
    μ = ν
    Lm = 999999999.0
    Lcrd, not_used = AISIS100.v16.app23334(ho, μ, t, Ixf, xof, hxf, Cwf, Ixyf, Iyf, Lm)


    #distortional buckling, compression 
    P = 1.0
    Mxx = 0.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(1.0*Lcrd, 2.0*Lcrd, 9)
    distortional_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Pcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_P, eig))


    #distortional buckling, Mxx
    P = 0.0
    Mxx = 1.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.5*Lcrd, 5.0*Lcrd, 9)
    distortional_buckling_Mxx = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_Mxx, eig))


    #gather up everything 
    properties = UniStrut(section_inputs, geometry, gross_section_properties, net_section_properties, Lnp_H, Lnp_D, tg_H, tg_D, tg, td_H, td_D, td, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Myy_neg, Mcrℓ_yy_neg, local_buckling_Myy_pos, Mcrℓ_yy_pos, distortional_buckling_P, Pcrd, distortional_buckling_Mxx, Mcrd)

    return properties 

end




function unistrut_out_geometry(H, D, L1, L2, R, t, dh_H, dh_D, de_H, de_D, rib_depth, rib_length, rib_radius_start, rib_radius_peak)

    θ_rib = atan(rib_depth/(rib_length/4))

    # segments = [L2, L1, D, H/2-rib_length/2 + rib_length/4, rib_depth/sin(θ_rib), rib_depth/sin(θ_rib), H/2-rib_length/2 + rib_length/4, D, L1, L2]

    # θ = [0.0, -π/2, π, -π/2, -π/4, -3π/4, -π/2, 0.0, -π/2, π]
    # r = [R, R, R, rib_radius_start, rib_radius_peak, rib_radius_start, R, R, R]
    # n = [4, 4, 4, 3, 4, 4, 3, 4, 4, 4]
    # n_r = [3, 3, 3, 3, 3, 3, 3, 3, 3]

    segments = [L2, L1, D - t, H/2 - t -rib_length/2 + rib_length/4, rib_depth/sin(θ_rib), rib_depth/sin(θ_rib), H/2 - t - rib_length/2 + rib_length/4, D - t, L1, L2]

    θ = [0.0, -π/2, π, -π/2, -π/4, -3π/4, -π/2, 0.0, -π/2, π]
    r = [R, R, R-t, rib_radius_start, rib_radius_peak, rib_radius_start, R-t, R, R]
    n = [4, 4, 4, 3, 4, 4, 3, 4, 4, 4]
    n_r = [3, 3, 3, 3, 3, 3, 3, 3, 3]


    section_geometry = CrossSection.Geometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to right", offset = (D-L2, H+2*L1 - 2*t))

    x = [section_geometry.center[i][1] for i in eachindex(section_geometry.center)]
    y = [section_geometry.center[i][2] for i in eachindex(section_geometry.center)]

    nodes = [x y zeros(Float64, length(x))]

    #flange hole minus
    D_flat_minus = segments[8] + t/2 - 2*R 
    xloc = t/2 + (R-t/2) + D_flat_minus/n[8] * 1.5
    yloc = L1 - t/2
    zloc = 0.0
    atol_x = D_flat_minus/n[8] * 0.55
    atol_y = 0.0
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_D_minus = hole_node_index[2] - 1

    x[hole_node_index[1]] = de_D - dh_D/2
    x[hole_node_index[2]] = de_D + dh_D/2

    #flange hole plus
    D_flat_plus = segments[3] + t/2 - 2*R 
    xloc = t/2 + (R-t/2) + D_flat_plus/n[3] * 1.5
    yloc = L1 - t + H - t/2
    zloc = 0.0
    atol_x = D_flat_plus/n[3] * 0.55
    atol_y = 0.0
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_D_plus = hole_node_index[2] - 1

    x[hole_node_index[1]] = de_D + dh_D/2
    x[hole_node_index[2]] = de_D - dh_D/2

    D_hole_element_index = [hole_element_index_D_plus; hole_element_index_D_minus]

    #web hole plus

    # index_start = n[1] + n_r[1] + n[2] + n_r[2] + n[3] + n_r[3] + 1
    # index_end = sum(n) + sum(n_r) - n[end] - n_r[end] - n[end-1] - n_r[end-1] + 1 
    # H_flat = y[index_start] - y[index_end]
    H_flat = H - 2*R

    index_start = n[1] + n_r[1] + n[2] + n_r[2] + n[3] + n_r[3] + n[4] + 1
    index_end = sum(n) + sum(n_r) - n[end] - n_r[end] - n[end-1] - n_r[end-1] - n[end-2] - n_r[end-2]- n[end-3] + 1 

    centerline_rib_length = y[index_start] - y[index_end]

    xloc = t/2
    H_flat_from_rib = H_flat/2 - centerline_rib_length/2
    yloc = H/2 + (L1 - t) + centerline_rib_length/2 + H_flat_from_rib/n[4] * 1.5
    zloc = 0.0
    atol_x = 0.0
    atol_y = H_flat_from_rib/n[4] * 0.55
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_H_plus = hole_node_index[2] - 1

    y[hole_node_index[1]] = H - de_H + dh_H/2
    y[hole_node_index[2]] = H - de_H - dh_H/2


    #web hole minus
    xloc = t/2
    yloc = H/2 + (L1 - t) - centerline_rib_length/2 - H_flat_from_rib/n[4] * 1.5
    zloc = 0.0
    atol_x = 0.0
    atol_y = H_flat/n[4] * 0.55
    atol_z = 0.0 

    hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
    hole_element_index_H_minus = hole_node_index[2] - 1

    y[hole_node_index[1]] = de_H + dh_H/2
    y[hole_node_index[2]] = de_H - dh_H/2

    H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]

    geometry = (coordinates = section_geometry, x=x, y=y, D_hole_element_index = D_hole_element_index, H_hole_element_index=H_hole_element_index)

    return geometry

end




function unistrut_out(section_inputs)

    # input = UniStrutInput(H, D, L1, L2, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D, rib_depth, rib_length, rib_radius_start, rib_radius_peak)

    @unpack member_type, section_type, H, D, L1, L2, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D, rib_depth, rib_length, rib_radius_start, rib_radius_peak = section_inputs

    geometry = unistrut_out_geometry(H, D, L1, L2, R, t, dh_H, dh_D, de_H, de_D, rib_depth, rib_length, rib_radius_start, rib_radius_peak)

    #gross section properties 
    gross_section_properties = CrossSection.Properties.open_thin_walled(geometry.coordinates.center, fill(t, length(geometry.x)-1)) 

    #calculate reduced thickness at hole
    Lnp_H = hole_pitch_H - hole_length_H
    Lnp_D = hole_pitch_D - hole_length_D

    kg = 0.60
    kd = 0.80

    tg_H = RMI.v2021.eq8_2__1(kg, t, Lnp_H, hole_pitch_H)
    tg_D = RMI.v2021.eq8_2__1(kg, t, Lnp_D, hole_pitch_D)

    td_H = RMI.v2021.eq8_2__2(kd, t, Lnp_H, hole_pitch_H)
    td_D = RMI.v2021.eq8_2__2(kd, t, Lnp_D, hole_pitch_D)

    #define cross-section element thicknesses
    num_elem = length(geometry.x) - 1
    tg = fill(t, num_elem)
    tg[geometry.D_hole_element_index] .= tg_D
    tg[geometry.H_hole_element_index] .= tg_H

    td = fill(t, num_elem)
    td[geometry.D_hole_element_index] .= td_D
    td[geometry.H_hole_element_index] .= td_H

    #net section properties 
    xy_coords_with_holes = [[geometry.x[i], geometry.y[i]] for i in eachindex(geometry.x)]
    net_section_properties = CrossSection.Properties.open_thin_walled(xy_coords_with_holes, tg) 
    # net_section_properties = CrossSection.Properties.open_thin_walled(geometry.coordinates.center, tg) 


    #elastic buckling properties 

    #local buckling, compression
    P = 1.0
    Mxx = 0.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []
    neigs = 1
    lengths = range(1.0*minimum([H, D]), 1.75*minimum([H,D]), 5)
    local_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    eig = 1
    Pcrℓ = minimum(CUFSM.Tools.get_load_factor(local_buckling_P, eig))

    #local buckling, Mxx 
    P = 0.0
    Mxx = 1.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 5)
    local_buckling_Mxx = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrℓ_xx  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Mxx, eig))


    #local buckling, Myy_neg 
    P = 0.0
    Mxx = 0.0
    Myy = 1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*minimum([H, D]), 1.5*minimum([H,D]), 5)
    local_buckling_Myy_neg = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrℓ_yy_neg  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Myy_neg, eig))

    #local buckling, Myy_pos
    P = 0.0
    Mxx = 0.0
    Myy = -1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.25 * D, 0.75 * D, 5)
    local_buckling_Myy_pos = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrℓ_yy_pos  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Myy_pos, eig))

    #distortional buckling 

    #find approximate distortional buckling half-wavelength
    CorZ = 0
    b = D
    d = L1
    θ = 90.0
    Af,Jf,Ixf,Iyf,Ixyf,Cwf,xof,hxf,hyf,yof = AISIS100.v16.table23131(CorZ,t,b,d,θ)

    ho = H
    μ = ν
    Lm = 999999999.0
    Lcrd, not_used = AISIS100.v16.app23334(ho, μ, t, Ixf, xof, hxf, Cwf, Ixyf, Iyf, Lm)


    #distortional buckling, compression 
    P = 1.0
    Mxx = 0.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(1.0*Lcrd, 2.0*Lcrd, 9)
    distortional_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Pcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_P, eig))


    #distortional buckling, Mxx
    P = 0.0
    Mxx = 1.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.5*Lcrd, 5.0*Lcrd, 9)
    distortional_buckling_Mxx = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_Mxx, eig))


    #gather up everything 
    properties = UniStrut(section_inputs, geometry, gross_section_properties, net_section_properties, Lnp_H, Lnp_D, tg_H, tg_D, tg, td_H, td_D, td, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Myy_neg, Mcrℓ_yy_neg, local_buckling_Myy_pos, Mcrℓ_yy_pos, distortional_buckling_P, Pcrd, distortional_buckling_Mxx, Mcrd)

    return properties 

end


end #module