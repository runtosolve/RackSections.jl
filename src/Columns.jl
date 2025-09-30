module Columns 

using CrossSectionGeometry, SectionProperties, RMI, LinesCurvesNodes, CUFSM, AISIS100, LinearAlgebra

struct CeeLipsInput

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

struct CeeLipsRibInput

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
    H_rib::Float64
    R_rib_flat::Float64
    R_rib_peak::Float64
    
end


struct CeeLips

    input::CeeLipsInput
    geometry::@NamedTuple{coordinates::CrossSectionGeometry.ThinWalled, x::Vector{Float64}, y::Vector{Float64}, D_hole_element_index::Vector{Int64}, H_hole_element_index::Vector{Int64}}
    properties::SectionProperties.SectionPropertiesObject
    net_properties::SectionProperties.SectionPropertiesObject
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
    distortional_buckling_Myy_pos::CUFSM.Model
    Mcrd_yy_pos::Float64

end


struct CeeLipsRib

    input::CeeLipsRibInput
    geometry::@NamedTuple{coordinates::CrossSectionGeometry.ThinWalled, x::Vector{Float64}, y::Vector{Float64}, D_hole_element_index::Vector{Int64}, H_hole_element_index::Vector{Int64}}
    properties::SectionProperties.SectionPropertiesObject
    net_properties::SectionProperties.SectionPropertiesObject
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
    distortional_buckling_Myy_pos::CUFSM.Model
    Mcrd_yy_pos::Float64

end


struct RectangularTubeInput

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

struct RectangularTube

    input::RectangularTubeInput
    geometry::@NamedTuple{coordinates::CrossSectionGeometry.ThinWalled, x::Vector{Float64}, y::Vector{Float64}, D_hole_element_index::Vector{Int64}, H_hole_element_index::Vector{Int64}}
    properties::SectionProperties.SectionPropertiesObject
    net_properties::SectionProperties.SectionPropertiesObject
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




struct HatLipsInput

    H::Float64
    D1::Float64
    D2::Float64
    D::Float64
    A::Float64
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

end


struct HatLipsRibInput

    H::Float64
    D1::Float64
    D2::Float64
    D::Float64
    A::Float64
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
    H_rib::Float64
    R_rib_peak::Float64
    R_rib_flat::Float64

end


struct HatLipsTrapezoidalRibInput

    H::Float64
    D1::Float64
    D2::Float64
    D::Float64
    A1::Float64
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
    H_rib::Float64
    W_rib::Float64
end



struct HatLips

    input::HatLipsInput
    geometry::@NamedTuple{coordinates::CrossSectionGeometry.ThinWalled, x::Vector{Float64}, y::Vector{Float64}, D1_hole_element_index::Vector{Int64}, D2_hole_element_index::Vector{Int64}, H_hole_element_index::Vector{Int64}}
    properties::SectionProperties.SectionPropertiesObject
    net_properties::SectionProperties.SectionPropertiesObject
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


struct HatLipsRib

    input::HatLipsRibInput
    geometry::@NamedTuple{coordinates::CrossSectionGeometry.ThinWalled, x::Vector{Float64}, y::Vector{Float64}, D1_hole_element_index::Vector{Int64}, D2_hole_element_index::Vector{Int64}, H_hole_element_index::Vector{Int64}}
    properties::SectionProperties.SectionPropertiesObject
    net_properties::SectionProperties.SectionPropertiesObject
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




struct HatLipsTrapezoidalRib

    input::HatLipsTrapezoidalRibInput
    geometry::@NamedTuple{coordinates::CrossSectionGeometry.ThinWalled, x::Vector{Float64}, y::Vector{Float64}, D1_hole_element_index::Vector{Int64}, D2_hole_element_index::Vector{Int64}, H_hole_element_index::Vector{Int64}}
    properties::SectionProperties.SectionPropertiesObject
    net_properties::SectionProperties.SectionPropertiesObject
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



struct HatInput

    H::Float64
    D1::Float64
    D2::Float64
    D::Float64
    A::Float64
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

end



struct HatRibInput

    H::Float64
    D1::Float64
    D2::Float64
    D::Float64
    A::Float64
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
    H_rib::Float64
    R_rib_peak::Float64
    R_rib_flat::Float64

end







struct Hat

    input::HatInput
    geometry::@NamedTuple{coordinates::CrossSectionGeometry.ThinWalled, x::Vector{Float64}, y::Vector{Float64}, D1_hole_element_index::Vector{Int64}, D2_hole_element_index::Vector{Int64}, H_hole_element_index::Vector{Int64}}
    properties::SectionProperties.SectionPropertiesObject
    net_properties::SectionProperties.SectionPropertiesObject
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




struct HatRib

    input::HatRibInput
    geometry::@NamedTuple{coordinates::CrossSectionGeometry.ThinWalled, x::Vector{Float64}, y::Vector{Float64}, D1_hole_element_index::Vector{Int64}, D2_hole_element_index::Vector{Int64}, H_hole_element_index::Vector{Int64}}
    properties::SectionProperties.SectionPropertiesObject
    net_properties::SectionProperties.SectionPropertiesObject
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



struct UniStrutInput

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

end


struct UniStrutRibInput

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
    H_rib::Float64
    R_rib_peak::Float64
    R_rib_flat::Float64

end



struct UniStrutIn

    input::UniStrutInput
    geometry::@NamedTuple{coordinates::CrossSectionGeometry.ThinWalled, x::Vector{Float64}, y::Vector{Float64}, D_hole_element_index::Vector{Int64}, H_hole_element_index::Vector{Int64}}
    properties::SectionProperties.SectionPropertiesObject
    net_properties::SectionProperties.SectionPropertiesObject
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
    Mcrd_xx::Float64
    distortional_buckling_Myy_pos::CUFSM.Model
    Mcrd_yy_pos::Float64


end


struct UniStrutInRib

    input::UniStrutRibInput
    geometry::@NamedTuple{coordinates::CrossSectionGeometry.ThinWalled, x::Vector{Float64}, y::Vector{Float64}, D_hole_element_index::Vector{Int64}, H_hole_element_index::Vector{Int64}}
    properties::SectionProperties.SectionPropertiesObject
    net_properties::SectionProperties.SectionPropertiesObject
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
    Mcrd_xx::Float64
    distortional_buckling_Myy_pos::CUFSM.Model
    Mcrd_yy_pos::Float64

end



struct UniStrutOut

    input::UniStrutInput
    geometry::@NamedTuple{coordinates::CrossSectionGeometry.ThinWalled, x::Vector{Float64}, y::Vector{Float64}, D_hole_element_index::Vector{Int64}, H_hole_element_index::Vector{Int64}}
    properties::SectionProperties.SectionPropertiesObject
    net_properties::SectionProperties.SectionPropertiesObject
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

end



struct UniStrutOutRib

    input::UniStrutRibInput
    geometry::@NamedTuple{coordinates::CrossSectionGeometry.ThinWalled, x::Vector{Float64}, y::Vector{Float64}, D_hole_element_index::Vector{Int64}, H_hole_element_index::Vector{Int64}}
    properties::SectionProperties.SectionPropertiesObject
    net_properties::SectionProperties.SectionPropertiesObject
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



function cee_with_lips(section_inputs)

   (; H, D, L, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D) = section_inputs


    geometry = cee_with_lips_geometry(H, D, L, R, t, dh_H, dh_D, de_H, de_D)

    #gross section properties 
    gross_section_properties = SectionProperties.open_thin_walled(geometry.coordinates.centerline_node_XY, fill(t, length(geometry.x)-1)) 

    #calculate reduced thickness at hole

    #initialize element thicknesses 
    num_elem = length(geometry.x) - 1
    tg = fill(t, num_elem) 
    td = fill(t, num_elem)

    kg = 0.60
    kd = 0.80

    if (hole_pitch_H != 0.0) & (hole_length_H != 0.0) & (dh_H != 0.0) & (de_H != 0.0)

        #calculate reduced thickness at hole
        Lnp_H = hole_pitch_H - hole_length_H
        tg_H = RMI.v2021.eq8_2__1(kg, t, Lnp_H, hole_pitch_H)
        td_H = RMI.v2021.eq8_2__2(kd, t, Lnp_H, hole_pitch_H)

        #define cross-section element thicknesses
        tg[geometry.H_hole_element_index] .= tg_H
        td[geometry.H_hole_element_index] .= td_H

    else 
        Lnp_H = 0.0
        tg_H = t
        td_H = t
    end

    if (hole_pitch_D != 0.0) & (hole_length_D != 0.0) & (dh_D != 0.0) & (de_D != 0.0)

        Lnp_D = hole_pitch_D - hole_length_D
        tg_D = RMI.v2021.eq8_2__1(kg, t, Lnp_D, hole_pitch_D)
        td_D = RMI.v2021.eq8_2__2(kd, t, Lnp_D, hole_pitch_D)
        
        #define cross-section element thicknesses
        tg[geometry.D_hole_element_index] .= tg_D
        td[geometry.D_hole_element_index] .= td_D

    else
        Lnp_D = 0.0
        td_D = t
        tg_D = t
    end



    #net section properties 
    xy_coords_with_holes = [[geometry.x[i], geometry.y[i]] for i in eachindex(geometry.x)]
    net_section_properties = SectionProperties.open_thin_walled(xy_coords_with_holes, tg) 


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
    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 7)
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

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 7)
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

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 7)
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

    lengths = range(0.25*minimum([H, D]), 1.25*minimum([H,D]), 7)
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


    #distortional buckling, Mcrd_yy_pos
    P = 0.0
    Mxx = 0.0
    Myy = -1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.5 * Lcrd, 2.0 * Lcrd, 9)
    distortional_buckling_Myy_pos = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrd_yy_pos  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_Mxx, eig))

    
    # lengths = range(0.5 * Lcrd, 1.5 * Lcrd, 9)

    # material = CeeSectionBuckling.Material(E, ν)
    # dimensions = CeeSectionBuckling.Dimensions(t, L, D, H, R-t)
    # section = CeeSectionBuckling.calculate_Mcrd_yy_pos(dimensions, material)
    # distortional_buckling_Myy_pos = section.results.model
    # Mcrd_yy_pos = section.results.Rcr


    #gather up everything 
    properties = CeeLips(section_inputs, geometry, gross_section_properties, net_section_properties, Lnp_H, Lnp_D, tg_H, tg_D, tg, td_H, td_D, td, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Myy_neg, Mcrℓ_yy_neg, local_buckling_Myy_pos, Mcrℓ_yy_pos, distortional_buckling_P, Pcrd, distortional_buckling_Mxx, Mcrd, distortional_buckling_Myy_pos, Mcrd_yy_pos)

    return properties 

end

function cee_with_lips_geometry(H, D, L, R, t, dh_H, dh_D, de_H, de_D)

    segments = [L, D, H, D, L]
    θ = [π/2, π, -π/2, 0.0, π/2]
    r = [R, R, R, R]
    n = [4, 4, 5, 4, 4]
    n_r = [3, 3, 3, 3] .* 3

    section_geometry = CrossSectionGeometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D, H-L))

    x = [section_geometry.centerline_node_XY[i][1] for i in eachindex(section_geometry.centerline_node_XY)]
    y = [section_geometry.centerline_node_XY[i][2] for i in eachindex(section_geometry.centerline_node_XY)]

    nodes = [x y zeros(Float64, length(x))]

    if (dh_D != 0) & (de_D != 0)

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

    else 

        hole_element_index_D_plus = 0
        hole_element_index_D_minus = 0

    end


    D_hole_element_index = [hole_element_index_D_plus; hole_element_index_D_minus]


    if (dh_H != 0) & (de_H != 0)

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

    else

        hole_element_index_H_plus = 0
        hole_element_index_H_minus = 0

    end

    H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]

    geometry = (coordinates = section_geometry, x=x, y=y, D_hole_element_index = D_hole_element_index, H_hole_element_index=H_hole_element_index)

    return geometry

end




function rectangular_tube_geometry(H, D, R, t, dh_H, dh_D, de_H, de_D)


    segments = [H, D, H, D]
    θ = [π/2, π, -π/2, 0.0]
    r = [R, R, R, R]
    n = [4, 4, 5, 4]
    n_r = [9, 9, 9, 9]

    section_geometry = CrossSectionGeometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D, 0.0))

    x = [section_geometry.centerline_node_XY[i][1] for i in eachindex(section_geometry.centerline_node_XY)]
    y = [section_geometry.centerline_node_XY[i][2] for i in eachindex(section_geometry.centerline_node_XY)]

    nodes = [x y zeros(Float64, length(x))]

    if (dh_D != 0.0) & (de_D != 0.0)

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
    
    else

        hole_element_index_D_plus = 0
        hole_element_index_D_minus = 0

    end


    D_hole_element_index = [hole_element_index_D_plus; hole_element_index_D_minus]

    if (dh_H != 0.0) & (de_H != 0.0)

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

    else

        hole_element_index_H_plus = 0
        hole_element_index_H_minus = 0

    end


    H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]

    geometry = (coordinates = section_geometry, x=x, y=y, D_hole_element_index = D_hole_element_index, H_hole_element_index=H_hole_element_index)

    return geometry

end



function rectangular_tube(section_inputs)

     (; H, D, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D) = section_inputs

    # input = RectangularTubeInput(H, D, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D)

    geometry = rectangular_tube_geometry(H, D, R, t, dh_H, dh_D, de_H, de_D)

    #gross section properties 
    gross_section_properties = SectionProperties.closed_thin_walled(geometry.coordinates.centerline_node_XY, fill(t, length(geometry.x))) 

    #remove NaNs to allow for writing to JSON
    gross_section_properties.xs = -1
    gross_section_properties.ys = -1
    gross_section_properties.Cw = -1
    gross_section_properties.B1 = -1
    gross_section_properties.B2 = -1
    gross_section_properties.wn = [-1, -1]

    ####

    num_elem = length(geometry.x)
    tg = fill(t, num_elem)
    kg = 0.60

    if (hole_pitch_H != 0.0) & (hole_length_H != 0.0) & (dh_H != 0.0) & (de_H != 0.0)

        #calculate reduced thickness at hole
        Lnp_H = hole_pitch_H - hole_length_H
        tg_H = RMI.v2021.eq8_2__1(kg, t, Lnp_H, hole_pitch_H)

        #define cross-section element thicknesses
        tg[geometry.H_hole_element_index] .= tg_H

    else 
        Lnp_H = 0.0
        tg_H = t
    end

    if (hole_pitch_D != 0.0) & (hole_length_D != 0.0) & (dh_D != 0.0) & ((de_D != 0.0))

        Lnp_D = hole_pitch_D - hole_length_D
        tg_D = RMI.v2021.eq8_2__1(kg, t, Lnp_D, hole_pitch_D)

       #define cross-section element thicknesses
        tg[geometry.D_hole_element_index] .= tg_D
 

    else
        Lnp_D = 0.0
        tg_D = 0.0
    end


    #net section properties 
    xy_coords_with_holes = [[geometry.x[i], geometry.y[i]] for i in eachindex(geometry.x)]
    net_section_properties = SectionProperties.closed_thin_walled(xy_coords_with_holes, tg) 

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
    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 7)
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

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 7)
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

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 7)
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

    lengths = range(0.25*minimum([H, D]), 1.25*minimum([H,D]), 7)
    local_buckling_Myy_pos = CUFSM.Tools.closed_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs)
    Mcrℓ_yy_pos  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Myy_pos, eig))

    #gather up everything 
    properties = RectangularTube(section_inputs, geometry, gross_section_properties, net_section_properties, Lnp_H, Lnp_D, tg_H, tg_D, tg, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Myy_neg, Mcrℓ_yy_neg, local_buckling_Myy_pos, Mcrℓ_yy_pos)

    return properties 

end


function cee_with_lips_rib_geometry(H, D, L, R, t, dh_H, dh_D, de_H, de_D, H_rib, R_rib_flat, R_rib_peak)


    #rib geometry, bottom of section  
    I = π / 2   #assume this is how long the rib arc is, hard coded, always 
    Tr = ((H_rib - t) - (R_rib_peak - t) * (1 - cos(I/2))) / sin(I/2)
    T = (R_rib_peak - t) * tan(I/2)
    T_total = Tr + T
    ΔTx = T_total * cos(I/2)


    segments = [L, D, H/2-ΔTx, T_total, T_total, H/2-ΔTx, D, L]

    θ = [π/2, π, -π/2, -π/2 + I/2, -π/2 - I/2, -π/2, 0.0, π/2]
    r = [R, R, R_rib_flat, R_rib_peak-t, R_rib_flat, R, R]
    n = [4, 4, 3, 4, 4, 3, 4, 4]
    n_r = [3, 3, 3, 3, 3, 3, 3] .* 3

    section_geometry = CrossSectionGeometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D, H-L))

    x = [section_geometry.centerline_node_XY[i][1] for i in eachindex(section_geometry.centerline_node_XY)]
    y = [section_geometry.centerline_node_XY[i][2] for i in eachindex(section_geometry.centerline_node_XY)]

    nodes = [x y zeros(Float64, length(x))]

    if (dh_D != 0) & (de_D != 0)

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

    else

        hole_element_index_D_plus = 0
        hole_element_index_D_minus = 0

    end

    D_hole_element_index = [hole_element_index_D_plus; hole_element_index_D_minus]

    if (dh_H != 0) & (de_H != 0)

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

    else 

        hole_element_index_H_plus = 0
        hole_element_index_H_minus = 0

    end

    H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]

    geometry = (coordinates = section_geometry, x=x, y=y, D_hole_element_index = D_hole_element_index, H_hole_element_index=H_hole_element_index)

    return geometry

end




function cee_with_lips_rib(section_inputs)

     (; H, D, L, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D, H_rib, R_rib_flat, R_rib_peak) = section_inputs


    geometry = cee_with_lips_rib_geometry(H, D, L, R, t, dh_H, dh_D, de_H, de_D, H_rib, R_rib_flat, R_rib_peak)

    #gross section properties 
    gross_section_properties = SectionProperties.open_thin_walled(geometry.coordinates.centerline_node_XY, fill(t, length(geometry.x)-1)) 

    # #calculate reduced thickness at hole

    num_elem = length(geometry.x) - 1
    tg = fill(t, num_elem) 
    td = fill(t, num_elem)

    kg = 0.60
    kd = 0.80

    if (hole_pitch_H != 0.0) & (hole_length_H != 0.0) & (dh_H != 0.0) & (de_H != 0.0)

        #calculate reduced thickness at hole
        Lnp_H = hole_pitch_H - hole_length_H
        tg_H = RMI.v2021.eq8_2__1(kg, t, Lnp_H, hole_pitch_H)
        td_H = RMI.v2021.eq8_2__2(kd, t, Lnp_H, hole_pitch_H)

        #define cross-section element thicknesses
        tg[geometry.H_hole_element_index] .= tg_H
        td[geometry.H_hole_element_index] .= td_H

    else 
        Lnp_H = 0.0
        tg_H = t
        td_H = t
    end

    if (hole_pitch_D != 0.0) & (hole_length_D != 0.0) & (dh_D != 0.0) & ((de_D != 0.0))

        Lnp_D = hole_pitch_D - hole_length_D
        tg_D = RMI.v2021.eq8_2__1(kg, t, Lnp_D, hole_pitch_D)
        td_D = RMI.v2021.eq8_2__2(kd, t, Lnp_D, hole_pitch_D)
        
        #define cross-section element thicknesses
        tg[geometry.D_hole_element_index] .= tg_D
        td[geometry.D_hole_element_index] .= td_D

    else
        Lnp_D = 0.0
        td_D = t
        tg_D = t
    end
  
    #net section properties 
    xy_coords_with_holes = [[geometry.x[i], geometry.y[i]] for i in eachindex(geometry.x)]
    net_section_properties = SectionProperties.open_thin_walled(xy_coords_with_holes, tg) 
    # net_section_properties = CrossSection.Properties.open_thin_walled(geometry.coordinates.centerline_node_XY, tg) 


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
    lengths = range(1.25*minimum([H, D]), 2.25*minimum([H,D]), 7)
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

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 7)
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

    lengths = range(1.0*minimum([H, D]), 2.0*minimum([H,D]), 7)
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

    lengths = range(0.25*minimum([H, D]), 1.25*minimum([H,D]), 7)
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



    #distortional buckling, Mcrd_yy_pos
    P = 0.0
    Mxx = 0.0
    Myy = -1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.5 * Lcrd, 2.0 * Lcrd, 9)
    distortional_buckling_Myy_pos = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrd_yy_pos  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_Mxx, eig))

    
 

    #gather up everything 
    properties = CeeLipsRib(section_inputs, geometry, gross_section_properties, net_section_properties, Lnp_H, Lnp_D, tg_H, tg_D, tg, td_H, td_D, td, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Myy_neg, Mcrℓ_yy_neg, local_buckling_Myy_pos, Mcrℓ_yy_pos, distortional_buckling_P, Pcrd, distortional_buckling_Mxx, Mcrd, distortional_buckling_Myy_pos, Mcrd_yy_pos)



    # #gather up everything 
    # properties = CeeLipsRib(section_inputs, geometry, gross_section_properties, net_section_properties, Lnp_H, Lnp_D, tg_H, tg_D, tg, td_H, td_D, td, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Myy_neg, Mcrℓ_yy_neg, local_buckling_Myy_pos, Mcrℓ_yy_pos, distortional_buckling_P, Pcrd, distortional_buckling_Mxx, Mcrd)

    return properties 

end




function hat_with_lips_rib_geometry(H, D1, D2, D, A, L, R, t, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2, H_rib, R_rib_peak, R_rib_flat)

    A = deg2rad(A)  #convert 

    yc = t * cos(A/2)


    #rib geometry, bottom of section  
    I = π / 2   #assume this is how long the rib arc is, hard coded, always 

    Tr = ((H_rib - t) - (R_rib_peak - t) * (1 - cos(I/2))) / sin(I/2)
    T = (R_rib_peak - t) * tan(I/2)
    T_total = Tr + T
    ΔTx = T_total * cos(I/2)

    #calculate X 
    D3 = D - D2 - D1 
    X = (D3 + yc) * tan(π - A) + t

    segments = [L-t, D2-t-yc, (X-t)/sin(π-A), D1, H/2-ΔTx, T_total, T_total, H/2-ΔTx, D1, (X-t)/sin(π-A), D2-t-yc, L-t]

    θ = [-π/2, π, A, π, -π/2, -π/2 + I/2, -π/2 - I/2, -π/2, 0.0, π-A, 0.0, -π/2]
    r = [R-t, R-t, R, R, R_rib_flat, R_rib_peak - t, R_rib_flat, R, R, R-t, R-t]
    n = [4, 3, 4, 3, 3, 4, 4, 3, 3, 4, 3, 4]
    n_r = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3] .* 3

    section_geometry = CrossSectionGeometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D-t, H - (X - L)))

    x = [section_geometry.centerline_node_XY[i][1] for i in eachindex(section_geometry.centerline_node_XY)]
    y = [section_geometry.centerline_node_XY[i][2] for i in eachindex(section_geometry.centerline_node_XY)]


    if (de_D1 !=0) & (dh_D1 !=0)

        #D1 flange hole minus
        hole_node_index_start = sum(n[1:8]) + 1 + sum(n_r[1:8]) + 1
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

        #web hole plus
        hole_node_index_start = sum(n[1:4]) + 1 + sum(n_r[1:4]) + 1
        hole_node_index = [hole_node_index_start, hole_node_index_start+1]
        hole_element_index_H_plus = hole_node_index[2] - 1
        y[hole_node_index[1]] = H - de_H + dh_H/2
        y[hole_node_index[2]] = H - de_H - dh_H/2


        #web hole minus
        hole_node_index_start = sum(n[1:7]) + 1 + sum(n_r[1:7]) + 1
        hole_node_index = [hole_node_index_start, hole_node_index_start+1]
        hole_element_index_H_minus = hole_node_index[2] - 1
        y[hole_node_index[1]] = de_H + dh_H/2
        y[hole_node_index[2]] = de_H - dh_H/2

    else

        hole_element_index_H_plus = 0
        hole_element_index_H_minus = 0

    end

    H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]


    geometry = (coordinates = section_geometry, x=x, y=y, D1_hole_element_index = D1_hole_element_index, D2_hole_element_index = D2_hole_element_index, H_hole_element_index=H_hole_element_index)

    return geometry

end




function hat_with_lips_rib(section_inputs)

     (; H, D1, D2, D, A, L, R, t, E, ν, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2, hole_pitch_H, hole_pitch_D1, hole_pitch_D2, hole_length_H, hole_length_D1, hole_length_D2, H_rib, R_rib_peak, R_rib_flat) = section_inputs


    geometry = hat_with_lips_rib_geometry(H, D1, D2, D, A, L, R, t, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2, H_rib, R_rib_peak, R_rib_flat)
    #gross section properties 
    gross_section_properties = SectionProperties.open_thin_walled(geometry.coordinates.centerline_node_XY, fill(t, length(geometry.x)-1)) 

    #calculate reduced thickness at hole

    num_elem = length(geometry.x) - 1
    tg = fill(t, num_elem)
    td = fill(t, num_elem)

    kg = 0.60
    kd = 0.80

    if (hole_pitch_H != 0.0) & (hole_length_H != 0.0) & (dh_H != 0.0) & (de_H != 0.0)

        Lnp_H = hole_pitch_H - hole_length_H
        tg_H = RMI.v2021.eq8_2__1(kg, t, Lnp_H, hole_pitch_H)
        td_H = RMI.v2021.eq8_2__2(kd, t, Lnp_H, hole_pitch_H)
        tg[geometry.H_hole_element_index] .= tg_H
        td[geometry.H_hole_element_index] .= td_H

    else
        Lnp_H = 0.0
        tg_H = t
        td_H = t
    end

    if (hole_pitch_D1 != 0.0) & (hole_length_D1 != 0.0) & (dh_D1 != 0.0) & (de_D1 != 0.0)

        Lnp_D1 = hole_pitch_D1 - hole_length_D1
        tg_D1 = RMI.v2021.eq8_2__1(kg, t, Lnp_D1, hole_pitch_D1)
        td_D1 = RMI.v2021.eq8_2__2(kd, t, Lnp_D1, hole_pitch_D1)
        tg[geometry.D1_hole_element_index] .= tg_D1
        td[geometry.D1_hole_element_index] .= td_D1

    else
        Lnp_D1 = 0.0
        td_D1 = t 
        tg_D1 = t
    end

    if (hole_pitch_D2 != 0.0) & (hole_length_D2 != 0.0) & (dh_D2 != 0.0) & (de_D2 != 0.0)
        Lnp_D2 = hole_pitch_D2 - hole_length_D2
        tg_D2 = RMI.v2021.eq8_2__1(kg, t, Lnp_D2, hole_pitch_D2)
        td_D2 = RMI.v2021.eq8_2__2(kd, t, Lnp_D2, hole_pitch_D2)
        tg[geometry.D2_hole_element_index] .= tg_D2
        td[geometry.D2_hole_element_index] .= td_D2
    else
        Lnp_D2 = 0.0
        td_D2 = t 
        tg_D2 = t
    end


    #net section properties 
    xy_coords_with_holes = [[geometry.x[i], geometry.y[i]] for i in eachindex(geometry.x)]
    net_section_properties = SectionProperties.open_thin_walled(xy_coords_with_holes, tg) 


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
    lengths = range(1.0*minimum([H, D]), 3.0*minimum([H,D]), 7)
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

    lengths = range(0.25*minimum([H, D]), 0.6*minimum([H,D]), 7)
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

    lengths = range(1.0*minimum([H, D]), 2.0*minimum([H,D]), 7)
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

    lengths = range(4.0*minimum([H, D]), 8.0*minimum([H,D]), 7)
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


function hat_with_rib_geometry(H, D1, D2, D, A, R, t, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2, H_rib, R_rib_peak, R_rib_flat)

        A = deg2rad(A)  #convert 

        yc = t * cos(A/2)
        
        
        #rib geometry, bottom of section  
        I = π / 2   #assume this is how long the rib arc is, hard coded, always 
        
        Tr = ((H_rib - t) - (R_rib_peak - t) * (1 - cos(I/2))) / sin(I/2)
        T = (R_rib_peak - t) * tan(I/2)
        T_total = Tr + T
        ΔTx = T_total * cos(I/2)
        
        #calculate X 
        D3 = D - D2 - D1 
        X = (D3 + yc) * tan(π - A) + t
        
        #define straight line segments 
        segments = [D2-yc, (X-t)/sin(π-A), D1, H/2-ΔTx, T_total, T_total, H/2-ΔTx, D1, (X-t)/sin(π-A), D2-yc]
        
        θ = [π, A, π, -π/2, -π/2 + I/2, -π/2 - I/2, -π/2, 0.0, π-A, 0.0]
        r = [R-t, R, R, R_rib_flat, R_rib_peak - t, R_rib_flat, R, R, R-t]
        n = [3, 4, 3, 3, 4, 4, 3, 3, 4, 3]
        n_r = [5, 5, 5, 5, 5, 5, 5, 5, 5]
        
        section_geometry = CrossSectionGeometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D, H - (X-t)))
        
        x = [section_geometry.centerline_node_XY[i][1] for i in eachindex(section_geometry.centerline_node_XY)]
        y = [section_geometry.centerline_node_XY[i][2] for i in eachindex(section_geometry.centerline_node_XY)]
        
        nodes = [x y zeros(Float64, length(x))]
        
        
        if (de_D1 != 0.0) & (dh_D1 != 0.0)

            #D1 flange hole minus
            hole_node_index_start = sum(n[1:7]) + 1 + sum(n_r[1:7]) + 1
            hole_node_index = [hole_node_index_start, hole_node_index_start + 1]
            hole_element_index_D1_minus = hole_node_index[2] - 1
            x[hole_node_index[1]] = de_D1 - dh_D1/2
            x[hole_node_index[2]] = de_D1 + dh_D1/2
            
            
            #D1 flange hole plus
            hole_node_index_start = sum(n[1:2]) + 1 + sum(n_r[1:2]) + 1
            hole_node_index = [hole_node_index_start, hole_node_index_start + 1]
            hole_element_index_D1_plus = hole_node_index[2] - 1
            x[hole_node_index[2]] = de_D1 - dh_D1/2
            x[hole_node_index[1]] = de_D1 + dh_D1/2

        else 

            hole_element_index_D1_plus = 0
            hole_element_index_D1_minus = 0

        end

        
        D1_hole_element_index = [hole_element_index_D1_plus; hole_element_index_D1_minus]
        
        
        if (de_D2 != 0.0) & (dh_D2 != 0.0)

            #D2 flange hole minus
            hole_node_index_start = sum(n[1:end-1]) + 1 + sum(n_r[1:end]) + 1
            hole_node_index = [hole_node_index_start, hole_node_index_start + 1]
            hole_element_index_D2_minus = hole_node_index[2] - 1
            x[hole_node_index[1]] = D - de_D2 - dh_D2/2
            x[hole_node_index[2]] = D - de_D2 + dh_D2/2
            
            #D2 flange hole plus
            hole_node_index = [2, 3]
            hole_element_index_D2_plus = hole_node_index[2] - 1
            x[hole_node_index[2]] = D - de_D2 - dh_D2/2
            x[hole_node_index[1]] = D - de_D2 + dh_D2/2

        else

            hole_element_index_D2_plus = 0
            hole_element_index_D2_minus = 0

        end
        
        D2_hole_element_index = [hole_element_index_D2_plus; hole_element_index_D2_minus]
        

        if (de_H != 0.0) & (dh_H != 0.0)

            #web hole plus
            hole_node_index_start = sum(n[1:3]) + 1 + sum(n_r[1:3]) + 1
            hole_node_index = [hole_node_index_start, hole_node_index_start+1]
            hole_element_index_H_plus = hole_node_index[2] - 1
            y[hole_node_index[1]] = H - de_H + dh_H/2
            y[hole_node_index[2]] = H - de_H - dh_H/2
            
            #web hole minus
            hole_node_index_start = sum(n[1:6]) + 1 + sum(n_r[1:6]) + 1
            hole_node_index = [hole_node_index_start, hole_node_index_start+1]
            hole_element_index_H_minus = hole_node_index[2] - 1
            y[hole_node_index[1]] = de_H + dh_H/2
            y[hole_node_index[2]] = de_H - dh_H/2

        else

            hole_element_index_H_plus = 0
            hole_element_index_H_minus = 0

        end

        
        H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]
        
        
        geometry = (coordinates = section_geometry, x=x, y=y, D1_hole_element_index = D1_hole_element_index, D2_hole_element_index = D2_hole_element_index, H_hole_element_index=H_hole_element_index)
        

    return geometry

end

#####


function hat_geometry(H, D1, D2, D, A, R, t, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2)

    #no rib, no lips!

        A = deg2rad(A)  #convert 

        yc = t * cos(A/2)
        

        #calculate X 
        D3 = D - D2 - D1 
        X = (D3 + yc) * tan(π - A) + t
        
        #define straight line segments 
        segments = [D2-yc, (X-t)/sin(π-A), D1, H, D1, (X-t)/sin(π-A), D2-yc]
        
        θ = [π, A, π, -π/2, 0.0, π-A, 0.0]
        r = [R-t, R, R, R, R, R-t]
        n = [3, 4, 3, 6, 3, 4, 3]
        n_r = [5, 5, 5, 5, 5, 5]
        
        section_geometry = CrossSectionGeometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D, H - (X-t)))
        
        x = [section_geometry.centerline_node_XY[i][1] for i in eachindex(section_geometry.centerline_node_XY)]
        y = [section_geometry.centerline_node_XY[i][2] for i in eachindex(section_geometry.centerline_node_XY)]
        
  
        if (de_D1 != 0.0) & (dh_D1 != 0.0)

            #D1 flange hole minus
            hole_node_index_start = sum(n[1:7]) + 1 + sum(n_r[1:7]) + 1
            hole_node_index = [hole_node_index_start, hole_node_index_start + 1]
            hole_element_index_D1_minus = hole_node_index[2] - 1
            x[hole_node_index[1]] = de_D1 - dh_D1/2
            x[hole_node_index[2]] = de_D1 + dh_D1/2
            
            
            #D1 flange hole plus
            hole_node_index_start = sum(n[1:2]) + 1 + sum(n_r[1:2]) + 1
            hole_node_index = [hole_node_index_start, hole_node_index_start + 1]
            hole_element_index_D1_plus = hole_node_index[2] - 1
            x[hole_node_index[2]] = de_D1 - dh_D1/2
            x[hole_node_index[1]] = de_D1 + dh_D1/2

        else 

            hole_element_index_D1_plus = 0
            hole_element_index_D1_minus = 0

        end

        
        D1_hole_element_index = [hole_element_index_D1_plus; hole_element_index_D1_minus]
        
        
        if (de_D2 != 0.0) & (dh_D2 != 0.0)

            #D2 flange hole minus
            hole_node_index_start = sum(n[1:end-1]) + 1 + sum(n_r[1:end]) + 1
            hole_node_index = [hole_node_index_start, hole_node_index_start + 1]
            hole_element_index_D2_minus = hole_node_index[2] - 1
            x[hole_node_index[1]] = D - de_D2 - dh_D2/2
            x[hole_node_index[2]] = D - de_D2 + dh_D2/2
            
            #D2 flange hole plus
            hole_node_index = [2, 3]
            hole_element_index_D2_plus = hole_node_index[2] - 1
            x[hole_node_index[2]] = D - de_D2 - dh_D2/2
            x[hole_node_index[1]] = D - de_D2 + dh_D2/2

        else

            hole_element_index_D2_plus = 0
            hole_element_index_D2_minus = 0

        end
        
        D2_hole_element_index = [hole_element_index_D2_plus; hole_element_index_D2_minus]
        

        if (de_H != 0.0) & (dh_H != 0.0)

            #web hole plus
            hole_node_index_start = sum(n[1:3]) + 1 + sum(n_r[1:3]) + 1
            hole_node_index = [hole_node_index_start, hole_node_index_start+1]
            hole_element_index_H_plus = hole_node_index[2] - 1
            y[hole_node_index[1]] = H - de_H + dh_H/2
            y[hole_node_index[2]] = H - de_H - dh_H/2
            
            #web hole minus
            hole_node_index_start = sum(n[1:4]) + sum(n_r[1:3]) - 1
            hole_node_index = [hole_node_index_start, hole_node_index_start+1]
            
            hole_element_index_H_minus = hole_node_index[2] - 1
            y[hole_node_index[1]] = de_H + dh_H/2
            y[hole_node_index[2]] = de_H - dh_H/2

        else

            hole_element_index_H_plus = 0
            hole_element_index_H_minus = 0

        end

        
        H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]
      
       
        
        geometry = (coordinates = section_geometry, x=x, y=y, D1_hole_element_index = D1_hole_element_index, D2_hole_element_index = D2_hole_element_index, H_hole_element_index=H_hole_element_index)
        

    return geometry

end


####


function hat(section_inputs)


     (; H, D1, D2, D, A, R, t, E, ν, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2,
     
     hole_pitch_H,
    hole_pitch_D1,
    hole_pitch_D2,
    hole_length_H,
    hole_length_D1,
    hole_length_D2) = section_inputs

    geometry = hat_geometry(H, D1, D2, D, A, R, t, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2)
    #gross section properties 
    gross_section_properties = SectionProperties.open_thin_walled(geometry.coordinates.centerline_node_XY, fill(t, length(geometry.x)-1)) 

 

    num_elem = length(geometry.x) - 1
    tg = fill(t, num_elem)
    td = fill(t, num_elem)

    kg = 0.60
    kd = 0.80

   

    if (hole_pitch_H != 0.0) & (hole_length_H != 0.0) & (dh_H != 0.0) & (de_H != 0.0)

        Lnp_H = hole_pitch_H - hole_length_H
        tg_H = RMI.v2021.eq8_2__1(kg, t, Lnp_H, hole_pitch_H)
        td_H = RMI.v2021.eq8_2__2(kd, t, Lnp_H, hole_pitch_H)
        tg[geometry.H_hole_element_index] .= tg_H
        td[geometry.H_hole_element_index] .= td_H

        # print(geometry.H_hole_element_index)

    else
        Lnp_H = 0.0
        td_H = t 
        tg_H = t
    end

    if (hole_pitch_D1 != 0.0) & (hole_length_D1 != 0.0) & (dh_D1 != 0.0) & (de_D1 != 0.0)

        Lnp_D1 = hole_pitch_D1 - hole_length_D1
        tg_D1 = RMI.v2021.eq8_2__1(kg, t, Lnp_D1, hole_pitch_D1)
        td_D1 = RMI.v2021.eq8_2__2(kd, t, Lnp_D1, hole_pitch_D1)
        tg[geometry.D1_hole_element_index] .= tg_D1
        td[geometry.D1_hole_element_index] .= td_D1

    else
        Lnp_D1 = 0.0
        tg_D1 = t 
        td_D1 = t
    end

    if (hole_pitch_D2 != 0.0) & (hole_length_D2 != 0.0) & (dh_D2 != 0.0) & (de_D2 != 0.0)
        Lnp_D2 = hole_pitch_D2 - hole_length_D2
        tg_D2 = RMI.v2021.eq8_2__1(kg, t, Lnp_D2, hole_pitch_D2)
        td_D2 = RMI.v2021.eq8_2__2(kd, t, Lnp_D2, hole_pitch_D2)
        tg[geometry.D2_hole_element_index] .= tg_D2
        td[geometry.D2_hole_element_index] .= td_D2
    else
        Lnp_D2 = 0.0
        tg_D2 = t 
        td_D2 = t 
    end

    #net section properties 
    xy_coords_with_holes = [[geometry.x[i], geometry.y[i]] for i in eachindex(geometry.x)]
    net_section_properties = SectionProperties.open_thin_walled(xy_coords_with_holes, tg) 


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
    # lengths = range(0.75*maximum([H, D]), 1.75*maximum([H,D]), 7)
    lengths = range(0.5*maximum([H, D]), 1.25*maximum([H,D]), 7)

    supports = [[1, 1, 1, 1, 1], [length(geometry.x), 1, 1, 1, 1]]

    local_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, supports, neigs)
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

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 7)
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

    lengths = range(1.0*minimum([H, D]), 2.0*minimum([H,D]), 7)
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

    lengths = range(3.0*minimum([H, D]), 7.0*minimum([H,D]), 7)
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
    properties = Hat(section_inputs, geometry, gross_section_properties, net_section_properties, Lnp_H, Lnp_D1, Lnp_D2, tg_H, tg_D1, tg_D2, tg, td_H, td_D1, td_D2, td, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Myy_neg, Mcrℓ_yy_neg, distortional_buckling_Myy_pos, Mcrd_yy_pos, distortional_buckling_P, Pcrd, distortional_buckling_Mxx, Mcrd)

    return properties 

end







#####

function hat_with_rib(section_inputs)


     (; H, D1, D2, D, A, R, t, E, ν, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2,
     
     hole_pitch_H,
    hole_pitch_D1,
    hole_pitch_D2,
    hole_length_H,
    hole_length_D1,
    hole_length_D2,
     
     H_rib, R_rib_peak, R_rib_flat) = section_inputs


    geometry = hat_with_rib_geometry(H, D1, D2, D, A, R, t, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2, H_rib, R_rib_peak, R_rib_flat)
    #gross section properties 
    gross_section_properties = SectionProperties.open_thin_walled(geometry.coordinates.centerline_node_XY, fill(t, length(geometry.x)-1)) 


    num_elem = length(geometry.x) - 1
    tg = fill(t, num_elem)
    td = fill(t, num_elem)

    kg = 0.60
    kd = 0.80

    if (hole_pitch_H != 0.0) & (hole_length_H != 0.0) & (dh_H != 0.0) & (de_H != 0.0)

        Lnp_H = hole_pitch_H - hole_length_H
        tg_H = RMI.v2021.eq8_2__1(kg, t, Lnp_H, hole_pitch_H)
        td_H = RMI.v2021.eq8_2__2(kd, t, Lnp_H, hole_pitch_H)
        tg[geometry.H_hole_element_index] .= tg_H
        td[geometry.H_hole_element_index] .= td_H

    else
        Lnp_H = 0.0
        tg_H = 0.0
        td_H = 0.0
    end

    if (hole_pitch_D1 != 0.0) & (hole_length_D1 != 0.0) & (dh_D1 != 0.0) & (de_D1 != 0.0)

        Lnp_D1 = hole_pitch_D1 - hole_length_D1
        tg_D1 = RMI.v2021.eq8_2__1(kg, t, Lnp_D1, hole_pitch_D1)
        td_D1 = RMI.v2021.eq8_2__2(kd, t, Lnp_D1, hole_pitch_D1)
        tg[geometry.D1_hole_element_index] .= tg_D1
        td[geometry.D1_hole_element_index] .= td_D1

    else
        Lnp_D1 = 0.0
        tg_D1 = 0.0
        td_D1 = 0.0
    end

    if (hole_pitch_D2 != 0.0) & (hole_length_D2 != 0.0) & (dh_D2 != 0.0) & (de_D2 != 0.0)
        Lnp_D2 = hole_pitch_D2 - hole_length_D2
        tg_D2 = RMI.v2021.eq8_2__1(kg, t, Lnp_D2, hole_pitch_D2)
        td_D2 = RMI.v2021.eq8_2__2(kd, t, Lnp_D2, hole_pitch_D2)
        tg[geometry.D2_hole_element_index] .= tg_D2
        td[geometry.D2_hole_element_index] .= td_D2
    else
        Lnp_D2 = 0.0
        tg_D2 = 0.0
        td_D2 = 0.0
    end

    #net section properties 
    xy_coords_with_holes = [[geometry.x[i], geometry.y[i]] for i in eachindex(geometry.x)]
    net_section_properties = SectionProperties.open_thin_walled(xy_coords_with_holes, tg) 
    # net_section_properties = CrossSection.Properties.open_thin_walled(geometry.coordinates.centerline_node_XY, tg) 


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
    # lengths = range(0.75*maximum([H, D]), 1.75*maximum([H,D]), 7)
    lengths = range(0.25*maximum([H, D]), 0.75*maximum([H,D]), 7)

    supports = [[1, 1, 1, 1, 1], [length(geometry.x), 1, 1, 1, 1]]

    local_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, supports, neigs)
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

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 7)
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

    lengths = range(1.0*minimum([H, D]), 2.0*minimum([H,D]), 7)
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

    lengths = range(3.0*minimum([H, D]), 7.0*minimum([H,D]), 7)
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


function hat_with_lips_trapezoidal_rib_geometry(H, D1, D2, D, A1, L, R, t, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2, A2, H_rib, W_rib)

    A1 = deg2rad(A1)  #convert 
    A2 = deg2rad(A2)  #convert 

    yc1 = t * cos(A1/2)
    yc2 = t * cos(A2/2)

    #calculate X 
    D3 = D - D2 - D1 
    X = (D3 + yc1) * tan(π - A1) + t

    
    H2 = (H_rib - t) / tan(π - A2)

    H1 = (H - (W_rib - yc2 - yc2) - H2 - H2)/2


    segments = [L-t, D2-t-yc1, (X-t)/sin(π-A1), D1, H1, (H_rib - t) / sin(π - A2), W_rib - 2 * yc2, (H_rib - t) / sin(π - A2), H1, D1, (X-t)/sin(π-A1), D2-t-yc1, L-t]
    

    θ = [-π/2, π, A1, π, -π/2, -π/2 + (π - A2), -π/2, -π/2 - (π - A2), -π/2, 0.0, π-A1, 0.0, -π/2]
    r = [R-t, R-t, R, R, R, R-t, R-t, R, R, R, R-t, R-t]
    n = [4, 3, 4, 3, 3, 4, 4, 4, 3, 3, 4, 3, 4]
    n_r = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]

    section_geometry = CrossSectionGeometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D-t, H - (X-L)))

    x = [section_geometry.centerline_node_XY[i][1] for i in eachindex(section_geometry.centerline_node_XY)]
    y = [section_geometry.centerline_node_XY[i][2] for i in eachindex(section_geometry.centerline_node_XY)]

    nodes = [x y zeros(Float64, length(x))]

    #D1 flange hole minus


    if (de_D1 != 0.0) & (dh_D1 != 0.0)

        hole_node_index_start = sum(n[1:9]) + 1 + sum(n_r[1:9]) + 1
        hole_node_index = [hole_node_index_start, hole_node_index_start + 1]

        hole_element_index_D1_minus = hole_node_index[2] - 1

        x[hole_node_index[1]] = de_D1 - dh_D1/2
        x[hole_node_index[2]] = de_D1 + dh_D1/2

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


    #D2 flange hole minus

    if (de_D2 != 0.0) & (dh_D2 != 0.0)

        hole_node_index_start = sum(n[1:11]) + 1 + sum(n_r[1:11]) + 1
        hole_node_index = [hole_node_index_start, hole_node_index_start + 1]
        hole_element_index_D2_minus = hole_node_index[2] - 1

        x[hole_node_index[1]] = D - de_D2 - dh_D2/2
        x[hole_node_index[2]] = D - de_D2 + dh_D2/2

        hole_node_index_start = sum(n[1:1]) + 1 + sum(n_r[1:1]) + 1
        hole_node_index = [hole_node_index_start, hole_node_index_start + 1]
        hole_element_index_D2_plus = hole_node_index[2] - 1

        x[hole_node_index[2]] = D - de_D2 - dh_D2/2
        x[hole_node_index[1]] = D - de_D2 + dh_D2/2

    else

        hole_element_index_D2_plus = 0
        hole_element_index_D2_minus = 0

    end

    D2_hole_element_index = [hole_element_index_D2_plus; hole_element_index_D2_minus]



    if (de_H != 0.0) & (dh_H != 0.0)


        #web hole plus

        hole_node_index_start = sum(n[1:4]) + 1 + sum(n_r[1:4]) + 1
        hole_node_index = [hole_node_index_start, hole_node_index_start + 1]
        hole_element_index_H_plus = hole_node_index[2] - 1

        y[hole_node_index[1]] = H - de_H + dh_H/2
        y[hole_node_index[2]] = H - de_H - dh_H/2



        #web hole minus

        hole_node_index_start = sum(n[1:8]) + 1 + sum(n_r[1:8]) + 1
        hole_node_index = [hole_node_index_start, hole_node_index_start + 1]
        hole_element_index_H_minus = hole_node_index[2] - 1


        y[hole_node_index[1]] = de_H + dh_H/2
        y[hole_node_index[2]] = de_H - dh_H/2

    else

        hole_element_index_H_plus = 0
        hole_element_index_H_minus = 0

    end

    H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]

    geometry = (coordinates = section_geometry, x=x, y=y, D1_hole_element_index = D1_hole_element_index, D2_hole_element_index = D2_hole_element_index, H_hole_element_index=H_hole_element_index)

    return geometry

end



function hat_with_lips_trapezoidal_rib(section_inputs)


     (; H, D1, D2, D, A1, L, R, t, E, ν, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2, hole_pitch_H, hole_pitch_D1, hole_pitch_D2, hole_length_H, hole_length_D1, hole_length_D2, H_rib, A2, W_rib) = section_inputs

    geometry = hat_with_lips_trapezoidal_rib_geometry(H, D1, D2, D, A1, L, R, t, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2, A2, H_rib, W_rib)
    #gross section properties 
    gross_section_properties = SectionProperties.open_thin_walled(geometry.coordinates.centerline_node_XY, fill(t, length(geometry.x)-1)) 


    num_elem = length(geometry.x) - 1
    tg = fill(t, num_elem)
    td = fill(t, num_elem)

    kg = 0.60
    kd = 0.80

    if (hole_pitch_H != 0.0) & (hole_length_H != 0.0) & (dh_H != 0.0) & (de_H != 0.0)

        Lnp_H = hole_pitch_H - hole_length_H
        tg_H = RMI.v2021.eq8_2__1(kg, t, Lnp_H, hole_pitch_H)
        td_H = RMI.v2021.eq8_2__2(kd, t, Lnp_H, hole_pitch_H)
        tg[geometry.H_hole_element_index] .= tg_H
        td[geometry.H_hole_element_index] .= td_H

    else
        Lnp_H = 0.0
        tg_H = t 
        td_H = t
    end

    if (hole_pitch_D1 != 0.0) & (hole_length_D1 != 0.0) & (dh_D1 != 0.0) & (de_D1 != 0.0)

        Lnp_D1 = hole_pitch_D1 - hole_length_D1
        tg_D1 = RMI.v2021.eq8_2__1(kg, t, Lnp_D1, hole_pitch_D1)
        td_D1 = RMI.v2021.eq8_2__2(kd, t, Lnp_D1, hole_pitch_D1)
        tg[geometry.D1_hole_element_index] .= tg_D1
        td[geometry.D1_hole_element_index] .= td_D1

    else
        Lnp_D1 = 0.0
        tg_D1 = t 
        td_D1 = t
    end

    if (hole_pitch_D2 != 0.0) & (hole_length_D2 != 0.0) & (dh_D2 != 0.0) & (de_D2 != 0.0)
        Lnp_D2 = hole_pitch_D2 - hole_length_D2
        tg_D2 = RMI.v2021.eq8_2__1(kg, t, Lnp_D2, hole_pitch_D2)
        td_D2 = RMI.v2021.eq8_2__2(kd, t, Lnp_D2, hole_pitch_D2)
        tg[geometry.D2_hole_element_index] .= tg_D2
        td[geometry.D2_hole_element_index] .= td_D2
    else
        Lnp_D2 = 0.0
        tg_D2 = t 
        td_D2 = t
    end

    #net section properties 
    xy_coords_with_holes = [[geometry.x[i], geometry.y[i]] for i in eachindex(geometry.x)]
    net_section_properties = SectionProperties.open_thin_walled(xy_coords_with_holes, tg) 


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
    lengths = range(0.5*maximum([H, D]), 1.8*maximum([H,D]), 7)
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

    lengths = range(0.25*minimum([H, D]), 0.6*minimum([H,D]), 7)
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

    lengths = range(1.0*maximum([H, D]), 2.25*maximum([H,D]), 7)
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

    lengths = range(4.0*minimum([H, D]), 8.0*minimum([H,D]), 7)
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

    lengths = range(0.75*Lcrd, 2.0*Lcrd, 7)
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

    lengths = range(0.75*Lcrd, 2.0*Lcrd, 7)
    distortional_buckling_Mxx = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_Mxx, eig))


    #gather up everything 
    properties = HatLipsTrapezoidalRib(section_inputs, geometry, gross_section_properties, net_section_properties, Lnp_H, Lnp_D1, Lnp_D2, tg_H, tg_D1, tg_D2, tg, td_H, td_D1, td_D2, td, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Myy_neg, Mcrℓ_yy_neg, distortional_buckling_Myy_pos, Mcrd_yy_pos, distortional_buckling_P, Pcrd, distortional_buckling_Mxx, Mcrd)

    return properties 

end

#####



function hat_with_lips_geometry(H, D1, D2, D, A, L, R, t, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2)

    A = deg2rad(A)  #convert 

    yc = t * cos(A/2)

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


    if (de_D1 !=0) & (dh_D1 !=0)

        #D1 flange hole minus
        hole_node_index_start = sum(n[1:8]) + 1 + sum(n_r[1:8]) + 1
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

        #web hole plus
        hole_node_index_start = sum(n[1:4]) + 1 + sum(n_r[1:4]) + 1
        hole_node_index = [hole_node_index_start, hole_node_index_start+1]
        hole_element_index_H_plus = hole_node_index[2] - 1
        y[hole_node_index[1]] = H - de_H + dh_H/2
        y[hole_node_index[2]] = H - de_H - dh_H/2


        # #web hole minus
        # hole_node_index_start = sum(n[1:7]) + 1 + sum(n_r[1:7]) + 1
        hole_node_index_start = sum(n[1:5]) + 1 + sum(n_r[1:4]) - 2
        hole_node_index = [hole_node_index_start, hole_node_index_start+1]
        hole_element_index_H_minus = hole_node_index[2] - 1
        y[hole_node_index[1]] = de_H + dh_H/2
        y[hole_node_index[2]] = de_H - dh_H/2

    else

        hole_element_index_H_plus = 0
        hole_element_index_H_minus = 0

    end

    H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]


    geometry = (coordinates = section_geometry, x=x, y=y, D1_hole_element_index = D1_hole_element_index, D2_hole_element_index = D2_hole_element_index, H_hole_element_index=H_hole_element_index)

    return geometry

end



#####



function hat_with_lips(section_inputs)

     (; H, D1, D2, D, A, L, R, t, E, ν, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2, hole_pitch_H, hole_pitch_D1, hole_pitch_D2, hole_length_H, hole_length_D1, hole_length_D2) = section_inputs

    geometry = hat_with_lips_geometry(H, D1, D2, D, A, L, R, t, dh_H, dh_D1, dh_D2, de_H, de_D1, de_D2)
    #gross section properties 
    gross_section_properties = SectionProperties.open_thin_walled(geometry.coordinates.centerline_node_XY, fill(t, length(geometry.x)-1)) 


    num_elem = length(geometry.x) - 1
    tg = fill(t, num_elem)
    td = fill(t, num_elem)

    kg = 0.60
    kd = 0.80

    if (hole_pitch_H != 0.0) & (hole_length_H != 0.0) & (dh_H != 0.0) & (de_H != 0.0)

        Lnp_H = hole_pitch_H - hole_length_H
        tg_H = RMI.v2021.eq8_2__1(kg, t, Lnp_H, hole_pitch_H)
        td_H = RMI.v2021.eq8_2__2(kd, t, Lnp_H, hole_pitch_H)
        tg[geometry.H_hole_element_index] .= tg_H
        td[geometry.H_hole_element_index] .= td_H

    else
        Lnp_H = 0.0
        tg_H = t 
        td_H = t
    end

    if (hole_pitch_D1 != 0.0) & (hole_length_D1 != 0.0) & (dh_D1 != 0.0) & (de_D1 != 0.0)

        Lnp_D1 = hole_pitch_D1 - hole_length_D1
        tg_D1 = RMI.v2021.eq8_2__1(kg, t, Lnp_D1, hole_pitch_D1)
        td_D1 = RMI.v2021.eq8_2__2(kd, t, Lnp_D1, hole_pitch_D1)
        tg[geometry.D1_hole_element_index] .= tg_D1
        td[geometry.D1_hole_element_index] .= td_D1

    else
        Lnp_D1 = 0.0
        tg_D1 = t 
        td_D1 = t
    end

    if (hole_pitch_D2 != 0.0) & (hole_length_D2 != 0.0) & (dh_D2 != 0.0) & (de_D2 != 0.0)
        Lnp_D2 = hole_pitch_D2 - hole_length_D2
        tg_D2 = RMI.v2021.eq8_2__1(kg, t, Lnp_D2, hole_pitch_D2)
        td_D2 = RMI.v2021.eq8_2__2(kd, t, Lnp_D2, hole_pitch_D2)
        tg[geometry.D2_hole_element_index] .= tg_D2
        td[geometry.D2_hole_element_index] .= td_D2
    else
        Lnp_D2 = 0.0
        tg_D2 = t 
        td_D2 = t
    end

    #net section properties 
    xy_coords_with_holes = [[geometry.x[i], geometry.y[i]] for i in eachindex(geometry.x)]
    net_section_properties = SectionProperties.open_thin_walled(xy_coords_with_holes, tg) 


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
    lengths = range(1.0*minimum([H, D]), 3.0*minimum([H,D]), 7)
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

    lengths = range(0.25*minimum([H, D]), 0.6*minimum([H,D]), 7)
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

    lengths = range(1.0*minimum([H, D]), 2.0*minimum([H,D]), 7)
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

    lengths = range(4.0*minimum([H, D]), 8.0*minimum([H,D]), 7)
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
    properties = HatLips(section_inputs, geometry, gross_section_properties, net_section_properties, Lnp_H, Lnp_D1, Lnp_D2, tg_H, tg_D1, tg_D2, tg, td_H, td_D1, td_D2, td, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Myy_neg, Mcrℓ_yy_neg, distortional_buckling_Myy_pos, Mcrd_yy_pos, distortional_buckling_P, Pcrd, distortional_buckling_Mxx, Mcrd)

    return properties 

end



#####

function unistrut_in_with_rib_geometry(H, D, L1, L2, R, t, dh_H, dh_D, de_H, de_D, H_rib, R_rib_flat, R_rib_peak)

    #rib geometry, bottom of section  
    I = π / 2   #assume this is how long the rib arc is, hard coded, always 

    Tr = ((H_rib - t) - (R_rib_peak - t) * (1 - cos(I/2))) / sin(I/2)
    T = (R_rib_peak - t) * tan(I/2)
    T_total = Tr + T
    ΔTx = T_total * cos(I/2)


    segments = [L2, L1, D, H/2-ΔTx,T_total, T_total, H/2-ΔTx, D, L1, L2]

    θ = [0.0, π/2, π, -π/2, -π/2 + I/2, -π/2 - I/2, -π/2, 0.0, π/2, π]
    r = [R, R, R, R_rib_flat, R_rib_peak - t, R_rib_flat, R, R, R]
    n = [4, 4, 4, 3, 4, 4, 3, 4, 4, 4]
    n_r = [5, 5, 5, 5, 5, 5, 5, 5, 5]

    section_geometry = CrossSectionGeometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D-L2, H-L1))

    x = [section_geometry.centerline_node_XY[i][1] for i in eachindex(section_geometry.centerline_node_XY)]
    y = [section_geometry.centerline_node_XY[i][2] for i in eachindex(section_geometry.centerline_node_XY)]

    nodes = [x y zeros(Float64, length(x))]

    if (de_D != 0.0) & (dh_D != 0)


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

    else 

        hole_element_index_D_plus = 0
        hole_element_index_D_minus = 0

    end


    D_hole_element_index = [hole_element_index_D_plus; hole_element_index_D_minus]

    #web hole plus

    if (de_H != 0.0) & (dh_H != 0)

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

    else 

        hole_element_index_H_plus = 0
        hole_element_index_H_minus = 0

    end

    H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]

    geometry = (coordinates = section_geometry, x=x, y=y, D_hole_element_index = D_hole_element_index, H_hole_element_index=H_hole_element_index)

    return geometry

end


#####


function unistrut_in_geometry(H, D, L1, L2, R, t, dh_H, dh_D, de_H, de_D)


    segments = [L2, L1, D, H, D, L1, L2]

    θ = [0.0, π/2, π, -π/2, 0.0, π/2, π]
    r = [R, R, R, R, R, R]
    n = [4, 4, 3, 6, 3, 4, 4]
    n_r = [5, 5, 5, 5, 5, 5]

    section_geometry = CrossSectionGeometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D-L2, H-L1))

    x = [section_geometry.centerline_node_XY[i][1] for i in eachindex(section_geometry.centerline_node_XY)]
    y = [section_geometry.centerline_node_XY[i][2] for i in eachindex(section_geometry.centerline_node_XY)]

    nodes = [x y zeros(Float64, length(x))]

    if (de_D != 0.0) & (dh_D != 0)


        #flange hole minus
        D_flat = D - 2*R 
        xloc = D/2
        yloc = t/2
        zloc = 0.0
        atol_x = D_flat/3/2 * 1.05
        atol_y = 0.0
        atol_z = 0.0 

        hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
        hole_element_index_D_minus = hole_node_index[2] - 1

        x[hole_node_index[1]] = de_D - dh_D/2
        x[hole_node_index[2]] = de_D + dh_D/2

        #flange hole plus
        # D_flat_plus = D - 2*R 
        xloc = D/2
        yloc = H - t/2
        zloc = 0.0
        atol_x = D_flat/3/2 * 1.05
        atol_y = 0.0
        atol_z = 0.0 

        hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
        hole_element_index_D_plus = hole_node_index[2] - 1

        x[hole_node_index[1]] = de_D + dh_D/2
        x[hole_node_index[2]] = de_D - dh_D/2

    else 

        hole_element_index_D_plus = 0
        hole_element_index_D_minus = 0

    end


    D_hole_element_index = [hole_element_index_D_plus; hole_element_index_D_minus]

    #web hole plus

    if (de_H != 0.0) & (dh_H != 0)

        H_flat = H - 2*R

        xloc = t/2
        
        yloc = R + 3/4 * H_flat 
        zloc = 0.0
        atol_x = 0.0
        atol_y = H_flat / 2 / 3 / 2 * 1.05
        atol_z = 0.0 

        hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
        hole_element_index_H_plus = hole_node_index[2] - 1

        y[hole_node_index[1]] = H - de_H + dh_H/2
        y[hole_node_index[2]] = H - de_H - dh_H/2


        #web hole minus
        xloc = t/2
        yloc = R + 1/4 * H_flat
        zloc = 0.0
        atol_x = 0.0
        atol_y = H_flat / 2 / 3 / 2 * 1.05
        atol_z = 0.0 

        hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
        hole_element_index_H_minus = hole_node_index[2] - 1

        y[hole_node_index[1]] = de_H + dh_H/2
        y[hole_node_index[2]] = de_H - dh_H/2

    else 

        hole_element_index_H_plus = 0
        hole_element_index_H_minus = 0

    end

    H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]

    geometry = (coordinates = section_geometry, x=x, y=y, D_hole_element_index = D_hole_element_index, H_hole_element_index=H_hole_element_index)

    return geometry

end


#####



function unistrut_in_with_rib(section_inputs)

     (; H, D, L1, L2, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D, H_rib, R_rib_flat, R_rib_peak) = section_inputs

    geometry = unistrut_in_with_rib_geometry(H, D, L1, L2, R, t, dh_H, dh_D, de_H, de_D, H_rib, R_rib_flat, R_rib_peak)

    #gross section properties 
    gross_section_properties = SectionProperties.open_thin_walled(geometry.coordinates.centerline_node_XY, fill(t, length(geometry.x)-1)) 


    #calculate reduced thickness at hole

    #initialize element thicknesses 
    num_elem = length(geometry.x) - 1
    tg = fill(t, num_elem) 
    td = fill(t, num_elem)

    kg = 0.60
    kd = 0.80

    if (hole_pitch_H != 0.0) & (hole_length_H != 0.0) & (dh_H != 0.0) & (de_H != 0.0)

        #calculate reduced thickness at hole
        Lnp_H = hole_pitch_H - hole_length_H
        tg_H = RMI.v2021.eq8_2__1(kg, t, Lnp_H, hole_pitch_H)
        td_H = RMI.v2021.eq8_2__2(kd, t, Lnp_H, hole_pitch_H)

        #define cross-section element thicknesses
        tg[geometry.H_hole_element_index] .= tg_H
        td[geometry.H_hole_element_index] .= td_H

    else 
        Lnp_H = 0.0
        tg_H = t 
        td_H = t
    end

    if (hole_pitch_D != 0.0) & (hole_length_D != 0.0) & (dh_D != 0.0) & ((de_D != 0.0))

        Lnp_D = hole_pitch_D - hole_length_D
        tg_D = RMI.v2021.eq8_2__1(kg, t, Lnp_D, hole_pitch_D)
        td_D = RMI.v2021.eq8_2__2(kd, t, Lnp_D, hole_pitch_D)
        
        #define cross-section element thicknesses
        tg[geometry.D_hole_element_index] .= tg_D
        td[geometry.D_hole_element_index] .= td_D

    else
        Lnp_D = 0.0
        tg_D = t 
        td_D = t
    end


    #net section properties 
    xy_coords_with_holes = [[geometry.x[i], geometry.y[i]] for i in eachindex(geometry.x)]
    net_section_properties = SectionProperties.open_thin_walled(xy_coords_with_holes, tg) 
    # net_section_properties = CrossSection.Properties.open_thin_walled(geometry.coordinates.centerline_node_XY, tg) 


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
    lengths = range(1.5*minimum([H, D]), 2.5*minimum([H,D]), 7)
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

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 9)
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

    lengths = range(0.75*minimum([H, D]), 1.5*minimum([H,D]), 7)
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

    lengths = range(0.25 * D, 0.75 * D, 7)
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

    lengths = range(1.25*Lcrd, 2.0*Lcrd, 7)
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

    lengths = range(0.5*Lcrd, 5.0*Lcrd, 7)
    distortional_buckling_Mxx = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrd_xx  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_Mxx, eig))




    #distortional buckling, Mcrd_yy_pos
    P = 0.0
    Mxx = 0.0
    Myy = -1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75 * Lcrd, 2.2 * Lcrd, 9)
    distortional_buckling_Myy_pos = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrd_yy_pos  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_Mxx, eig))




    #gather up everything 
    properties = UniStrutInRib(section_inputs, geometry, gross_section_properties, net_section_properties, Lnp_H, Lnp_D, tg_H, tg_D, tg, td_H, td_D, td, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Myy_neg, Mcrℓ_yy_neg, local_buckling_Myy_pos, Mcrℓ_yy_pos, distortional_buckling_P, Pcrd, distortional_buckling_Mxx, Mcrd_xx, distortional_buckling_Myy_pos, Mcrd_yy_pos)

    return properties 

end





function unistrut_in(section_inputs)

    # input = UniStrutInput(H, D, L1, L2, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D, rib_depth, rib_length, rib_radius_start, rib_radius_peak)

     (; H, D, L1, L2, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D) = section_inputs

    geometry = unistrut_in_geometry(H, D, L1, L2, R, t, dh_H, dh_D, de_H, de_D)

    #gross section properties 
    gross_section_properties = SectionProperties.open_thin_walled(geometry.coordinates.centerline_node_XY, fill(t, length(geometry.x)-1)) 


    #calculate reduced thickness at hole

    #initialize element thicknesses 
    num_elem = length(geometry.x) - 1
    tg = fill(t, num_elem) 
    td = fill(t, num_elem)

    kg = 0.60
    kd = 0.80

    if (hole_pitch_H != 0.0) & (hole_length_H != 0.0) & (dh_H != 0.0) & (de_H != 0.0)

        #calculate reduced thickness at hole
        Lnp_H = hole_pitch_H - hole_length_H
        tg_H = RMI.v2021.eq8_2__1(kg, t, Lnp_H, hole_pitch_H)
        td_H = RMI.v2021.eq8_2__2(kd, t, Lnp_H, hole_pitch_H)

        #define cross-section element thicknesses
        tg[geometry.H_hole_element_index] .= tg_H
        td[geometry.H_hole_element_index] .= td_H

    else 
        Lnp_H = 0.0
        tg_H = t 
        td_H = t
    end

    if (hole_pitch_D != 0.0) & (hole_length_D != 0.0) & (dh_D != 0.0) & ((de_D != 0.0))

        Lnp_D = hole_pitch_D - hole_length_D
        tg_D = RMI.v2021.eq8_2__1(kg, t, Lnp_D, hole_pitch_D)
        td_D = RMI.v2021.eq8_2__2(kd, t, Lnp_D, hole_pitch_D)
        
        #define cross-section element thicknesses
        tg[geometry.D_hole_element_index] .= tg_D
        td[geometry.D_hole_element_index] .= td_D

    else
        Lnp_D = 0.0
        tg_D = t 
        td_D = t
    end


    #net section properties 
    xy_coords_with_holes = [[geometry.x[i], geometry.y[i]] for i in eachindex(geometry.x)]
    net_section_properties = SectionProperties.open_thin_walled(xy_coords_with_holes, tg) 


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
    lengths = range(1.5*minimum([H, D]), 2.5*minimum([H,D]), 7)
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

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 9)
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

    lengths = range(0.75*minimum([H, D]), 1.5*minimum([H,D]), 7)
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

    lengths = range(0.25 * D, 0.75 * D, 7)
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

    lengths = range(1.25*Lcrd, 2.0*Lcrd, 7)
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

    lengths = range(0.5*Lcrd, 5.0*Lcrd, 7)
    distortional_buckling_Mxx = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrd_xx  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_Mxx, eig))


    #distortional buckling, Mcrd_yy_pos
    P = 0.0
    Mxx = 0.0
    Myy = -1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.5 * Lcrd, 2.5 * Lcrd, 9)
    distortional_buckling_Myy_pos = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrd_yy_pos  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_Mxx, eig))




    #gather up everything 
    properties = UniStrutIn(section_inputs, geometry, gross_section_properties, net_section_properties, Lnp_H, Lnp_D, tg_H, tg_D, tg, td_H, td_D, td, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Myy_neg, Mcrℓ_yy_neg, local_buckling_Myy_pos, Mcrℓ_yy_pos, distortional_buckling_P, Pcrd, distortional_buckling_Mxx, Mcrd_xx, distortional_buckling_Myy_pos, Mcrd_yy_pos)

    return properties 

end




function unistrut_out_with_rib_geometry(H, D, L1, L2, R, t, dh_H, dh_D, de_H, de_D, H_rib, R_rib_flat, R_rib_peak)

    #rib geometry, bottom of section  
    I = π / 2   #assume this is how long the rib arc is, hard coded, always 

    Tr = ((H_rib - t) - (R_rib_peak - t) * (1 - cos(I/2))) / sin(I/2)
    T = (R_rib_peak - t) * tan(I/2)
    T_total = Tr + T
    ΔTx = T_total * cos(I/2)


    segments = [L2-t, L1 - 2*t, D-t, H/2-ΔTx, T_total, T_total, H/2-ΔTx, D-t, L1-2*t, L2 -t]

    θ = [0.0, -π/2, π, -π/2, -π/2 + I/2, -π/2 - I/2, -π/2, 0.0, -π/2, π]
    r = [R-t, R-t, R, R_rib_flat, R_rib_peak - t, R_rib_flat, R, R-t, R-t]
    n = [4, 4, 4, 3, 4, 4, 3, 4, 4, 4]
    n_r = [5, 5, 5, 5, 5, 5, 5, 5, 5]

    section_geometry = CrossSectionGeometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D-L2, H/2 - t + (L1 - t)))   

    x = [section_geometry.centerline_node_XY[i][1] for i in eachindex(section_geometry.centerline_node_XY)]
    y = [section_geometry.centerline_node_XY[i][2] for i in eachindex(section_geometry.centerline_node_XY)]

    nodes = [x y zeros(Float64, length(x))]

    if (de_D != 0.0) & (dh_D != 0.0)

        #flange hole minus
        D_flat_minus = segments[8] - (R-t) - (R-t) 
        xloc = R + D_flat_minus/n[8] * 1.5
        yloc = -H/2 + t/2
        zloc = 0.0
        atol_x = D_flat_minus/n[8] * 0.55
        atol_y = 0.0
        atol_z = 0.0 

        hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
        hole_element_index_D_minus = hole_node_index[2] - 1

        x[hole_node_index[1]] = de_D - dh_D/2
        x[hole_node_index[2]] = de_D + dh_D/2

        #flange hole plus
        D_flat_plus = segments[3] - (R-t) - (R-t) 
        xloc = R + D_flat_plus/n[3] * 1.5
        yloc = H/2 - t/2
        zloc = 0.0
        atol_x = D_flat_plus/n[3] * 0.55
        atol_y = 0.0
        atol_z = 0.0 

        hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
        hole_element_index_D_plus = hole_node_index[2] - 1

        x[hole_node_index[1]] = de_D + dh_D/2
        x[hole_node_index[2]] = de_D - dh_D/2

    else

        hole_element_index_D_plus = 0
        hole_element_index_D_minus = 0

    end

    D_hole_element_index = [hole_element_index_D_plus; hole_element_index_D_minus]


    if (de_H != 0.0) & (dh_H != 0.0)

        #web hole plus

        H_flat = H - 2*R

        index_start = n[1] + n_r[1] + n[2] + n_r[2] + n[3] + n_r[3] + n[4] + 1
        index_end = sum(n) + sum(n_r) - n[end] - n_r[end] - n[end-1] - n_r[end-1] - n[end-2] - n_r[end-2]- n[end-3] + 1 

        centerline_rib_length = y[index_start] - y[index_end]

        xloc = t/2
        H_flat_from_rib = H_flat/2 - centerline_rib_length/2
        yloc = centerline_rib_length/2 + H_flat_from_rib/n[4] * 1.5
        zloc = 0.0
        atol_x = 0.0
        atol_y = H_flat_from_rib/n[4] * 0.55
        atol_z = 0.0 

        hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
        hole_element_index_H_plus = hole_node_index[2] - 1

        y[hole_node_index[1]] = H/2 - de_H + dh_H/2
        y[hole_node_index[2]] = H/2 - de_H - dh_H/2

        #web hole minus
        xloc = t/2
        yloc =  -centerline_rib_length/2 - H_flat_from_rib/n[4] * 1.5
        zloc = 0.0
        atol_x = 0.0
        atol_y = H_flat_from_rib/n[4] * 0.55
        atol_z = 0.0 

        hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
        hole_element_index_H_minus = hole_node_index[2] - 1

        y[hole_node_index[1]] = -H/2 + de_H + dh_H/2
        y[hole_node_index[2]] = -H/2 + de_H - dh_H/2

    else
        hole_element_index_H_plus = 0
        hole_element_index_H_minus = 0
    end

    H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]

    geometry = (coordinates = section_geometry, x=x, y=y, D_hole_element_index = D_hole_element_index, H_hole_element_index=H_hole_element_index)

    return geometry

end





function unistrut_out_geometry(H, D, L1, L2, R, t, dh_H, dh_D, de_H, de_D)


    segments = [L2-t, L1 - 2*t, D-t, H, D-t, L1-2*t, L2 -t]

    θ = [0.0, -π/2, π, -π/2, 0.0, -π/2, π]
    r = [R-t, R-t, R, R, R-t, R-t]
    n = [4, 4, 3, 6, 3, 4, 4]
    n_r = [5, 5, 5, 5, 5, 5]

    section_geometry = CrossSectionGeometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D-L2, H/2 - t + (L1 - t)))   

    x = [section_geometry.centerline_node_XY[i][1] for i in eachindex(section_geometry.centerline_node_XY)]
    y = [section_geometry.centerline_node_XY[i][2] for i in eachindex(section_geometry.centerline_node_XY)]

    nodes = [x y zeros(Float64, length(x))]

    if (de_D != 0.0) & (dh_D != 0.0)

        #flange hole minus
        D_flat = D - 2 * R
        xloc = D/2
        yloc = -H/2 + t/2
        zloc = 0.0
        atol_x = D_flat/3/2 * 1.05
        atol_y = 0.0
        atol_z = 0.0 

        hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
        hole_element_index_D_minus = hole_node_index[2] - 1

        x[hole_node_index[1]] = de_D - dh_D/2
        x[hole_node_index[2]] = de_D + dh_D/2

        #flange hole plus
        xloc = D/2
        yloc = H/2 - t/2
        zloc = 0.0
        atol_x = D_flat/3/2 * 1.05
        atol_y = 0.0
        atol_z = 0.0 

        hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
        hole_element_index_D_plus = hole_node_index[2] - 1

        x[hole_node_index[1]] = de_D + dh_D/2
        x[hole_node_index[2]] = de_D - dh_D/2

    else

        hole_element_index_D_plus = 0
        hole_element_index_D_minus = 0

    end

    D_hole_element_index = [hole_element_index_D_plus; hole_element_index_D_minus]


    if (de_H != 0.0) & (dh_H != 0.0)

        #web hole plus

        H_flat = H - 2*R

        xloc = t/2
     
        yloc = H_flat / 4
        zloc = 0.0
        atol_x = 0.0
        atol_y = H_flat / 2 / 3 * 1.05
        atol_z = 0.0 

        hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
        hole_element_index_H_plus = hole_node_index[2] - 1

        y[hole_node_index[1]] = H/2 - de_H + dh_H/2
        y[hole_node_index[2]] = H/2 - de_H - dh_H/2

        #web hole minus
        xloc = t/2
        yloc =  -H_flat/4
        zloc = 0.0
        atol_x = 0.0
        atol_y = H_flat / 2 / 3  * 1.05
        atol_z = 0.0 

        hole_node_index = LinesCurvesNodes.find_nodes(nodes, xloc, yloc, zloc, atol_x, atol_y, atol_z)
        hole_element_index_H_minus = hole_node_index[2] - 1

        y[hole_node_index[1]] = -H/2 + de_H + dh_H/2
        y[hole_node_index[2]] = -H/2 + de_H - dh_H/2

    else
        hole_element_index_H_plus = 0
        hole_element_index_H_minus = 0
    end

    H_hole_element_index = [hole_element_index_H_plus; hole_element_index_H_minus]

    geometry = (coordinates = section_geometry, x=x, y=y, D_hole_element_index = D_hole_element_index, H_hole_element_index=H_hole_element_index)

    return geometry

end



function unistrut_out_with_rib(section_inputs)

     (; H, D, L1, L2, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D, H_rib, R_rib_flat, R_rib_peak) = section_inputs

    geometry = unistrut_out_with_rib_geometry(H, D, L1, L2, R, t, dh_H, dh_D, de_H, de_D, H_rib, R_rib_flat, R_rib_peak)

    #gross section properties 
    gross_section_properties = SectionProperties.open_thin_walled(geometry.coordinates.centerline_node_XY, fill(t, length(geometry.x)-1)) 


    #calculate reduced thickness at hole

    #initialize element thicknesses 
    num_elem = length(geometry.x) - 1
    tg = fill(t, num_elem) 
    td = fill(t, num_elem)

    kg = 0.60
    kd = 0.80

    if (hole_pitch_H != 0.0) & (hole_length_H != 0.0) & (dh_H != 0.0) & (de_H != 0.0)

        #calculate reduced thickness at hole
        Lnp_H = hole_pitch_H - hole_length_H
        tg_H = RMI.v2021.eq8_2__1(kg, t, Lnp_H, hole_pitch_H)
        td_H = RMI.v2021.eq8_2__2(kd, t, Lnp_H, hole_pitch_H)

        #define cross-section element thicknesses
        tg[geometry.H_hole_element_index] .= tg_H
        td[geometry.H_hole_element_index] .= td_H

    else 
        Lnp_H = 0.0
        tg_H = t 
        td_H = t
    end

    if (hole_pitch_D != 0.0) & (hole_length_D != 0.0) & (dh_D != 0.0) & ((de_D != 0.0))

        Lnp_D = hole_pitch_D - hole_length_D
        tg_D = RMI.v2021.eq8_2__1(kg, t, Lnp_D, hole_pitch_D)
        td_D = RMI.v2021.eq8_2__2(kd, t, Lnp_D, hole_pitch_D)
        
        #define cross-section element thicknesses
        tg[geometry.D_hole_element_index] .= tg_D
        td[geometry.D_hole_element_index] .= td_D

    else
        Lnp_D = 0.0
        tg_D = t 
        td_D = t
    end


    #net section properties 
    xy_coords_with_holes = [[geometry.x[i], geometry.y[i]] for i in eachindex(geometry.x)]
    net_section_properties = SectionProperties.open_thin_walled(xy_coords_with_holes, tg) 
 

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
    lengths = range(1.0*minimum([H, D]), 2.3*minimum([H,D]), 9)
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

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 7)
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

    lengths = range(0.75*minimum([H, D]), 2.0*minimum([H,D]), 7)
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

    lengths = range(0.25 * D, 0.75 * D, 7)
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
    # P = 1.0
    # Mxx = 0.0
    # Myy = 0.0
    # M11 = 0.0
    # M22 = 0.0
    # constraints = []
    # springs = []

    # lengths = range(1.25*Lcrd, 3.0*Lcrd, 9)
    # distortional_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    # Pcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_P, eig))

    


    # #find approximate distortional buckling half-wavelength
    # CorZ = 0
    # b = D
    # d = L1 + L2
    # θ = 90.0
    # Af,Jf,Ixf,Iyf,Ixyf,Cwf,xof,hxf,hyf,yof = AISIS100.v16.table23131(CorZ,t,b,d,θ)

    # ho = H
    # μ = ν
    # Lm = 999999999.0
    # Lcrd, not_used = AISIS100.v16.app23334(ho, μ, t, Ixf, xof, hxf, Cwf, Ixyf, Iyf, Lm)


    # #distortional buckling, compression 
    # P = 1.0
    # Mxx = 0.0
    # Myy = 0.0
    # M11 = 0.0
    # M22 = 0.0
    # constraints = []
    # springs = []

    # # lengths = range(1.0*Lcrd, 2.0*Lcrd, 9)
    # lengths = range(0.5*Lcrd, 6.0*Lcrd, 15)
    # distortional_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    # Pcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_P, eig))


    # #distortional buckling, Mxx
    # P = 0.0
    # Mxx = 1.0
    # Myy = 0.0
    # M11 = 0.0
    # M22 = 0.0
    # constraints = []
    # springs = []

    # lengths = range(0.5*Lcrd, 5.0*Lcrd, 9)
    # distortional_buckling_Mxx = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    # Mcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_Mxx, eig))


    #gather up everything 
    properties = UniStrutOutRib(section_inputs, geometry, gross_section_properties, net_section_properties, Lnp_H, Lnp_D, tg_H, tg_D, tg, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Myy_neg, Mcrℓ_yy_neg, local_buckling_Myy_pos, Mcrℓ_yy_pos)



    return properties 

end




function unistrut_out(section_inputs)

     (; H, D, L1, L2, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D) = section_inputs

    geometry = unistrut_out_geometry(H, D, L1, L2, R, t, dh_H, dh_D, de_H, de_D)

    #gross section properties 
    gross_section_properties = SectionProperties.open_thin_walled(geometry.coordinates.centerline_node_XY, fill(t, length(geometry.x)-1)) 


    #calculate reduced thickness at hole

    #initialize element thicknesses 
    num_elem = length(geometry.x) - 1
    tg = fill(t, num_elem) 
    td = fill(t, num_elem)

    kg = 0.60
    kd = 0.80

    if (hole_pitch_H != 0.0) & (hole_length_H != 0.0) & (dh_H != 0.0) & (de_H != 0.0)

        #calculate reduced thickness at hole
        Lnp_H = hole_pitch_H - hole_length_H
        tg_H = RMI.v2021.eq8_2__1(kg, t, Lnp_H, hole_pitch_H)
        td_H = RMI.v2021.eq8_2__2(kd, t, Lnp_H, hole_pitch_H)

        #define cross-section element thicknesses
        tg[geometry.H_hole_element_index] .= tg_H
        td[geometry.H_hole_element_index] .= td_H

    else 
        Lnp_H = 0.0
        tg_H = t 
        td_H = t
    end

    if (hole_pitch_D != 0.0) & (hole_length_D != 0.0) & (dh_D != 0.0) & ((de_D != 0.0))

        Lnp_D = hole_pitch_D - hole_length_D
        tg_D = RMI.v2021.eq8_2__1(kg, t, Lnp_D, hole_pitch_D)
        td_D = RMI.v2021.eq8_2__2(kd, t, Lnp_D, hole_pitch_D)
        
        #define cross-section element thicknesses
        tg[geometry.D_hole_element_index] .= tg_D
        td[geometry.D_hole_element_index] .= td_D

    else
        Lnp_D = 0.0
        tg_D = t 
        td_D = t
    end


    #net section properties 
    xy_coords_with_holes = [[geometry.x[i], geometry.y[i]] for i in eachindex(geometry.x)]
    net_section_properties = SectionProperties.open_thin_walled(xy_coords_with_holes, tg) 

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
    lengths = range(1.0*minimum([H, D]), 2.0*minimum([H,D]), 7)
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

    lengths = range(0.75*minimum([H, D]), 1.25*minimum([H,D]), 7)
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

    lengths = range(0.75*minimum([H, D]), 2.0*minimum([H,D]), 7)
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

    lengths = range(0.25 * D, 0.75 * D, 7)
    local_buckling_Myy_pos = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Mcrℓ_yy_pos  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Myy_pos, eig))

    #distortional buckling 

    

    # #find approximate distortional buckling half-wavelength
    # CorZ = 0
    # b = D
    # d = L1 + L2
    # θ = 90.0
    # Af,Jf,Ixf,Iyf,Ixyf,Cwf,xof,hxf,hyf,yof = AISIS100.v16.table23131(CorZ,t,b,d,θ)

    # ho = H
    # μ = ν
    # Lm = 999999999.0
    # Lcrd, not_used = AISIS100.v16.app23334(ho, μ, t, Ixf, xof, hxf, Cwf, Ixyf, Iyf, Lm)


    # #distortional buckling, compression 
    # P = 1.0
    # Mxx = 0.0
    # Myy = 0.0
    # M11 = 0.0
    # M22 = 0.0
    # constraints = []
    # springs = []

    # # lengths = range(1.0*Lcrd, 2.0*Lcrd, 9)
    # lengths = range(0.5*Lcrd, 6.0*Lcrd, 15)
    # distortional_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    # Pcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_P, eig))


    # #distortional buckling, Mxx
    # P = 0.0
    # Mxx = 1.0
    # Myy = 0.0
    # M11 = 0.0
    # M22 = 0.0
    # constraints = []
    # springs = []

    # lengths = range(0.5*Lcrd, 5.0*Lcrd, 9)
    # distortional_buckling_Mxx = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    # Mcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_Mxx, eig))


    #gather up everything 
    # properties = UniStrut(section_inputs, geometry, gross_section_properties, net_section_properties, Lnp_H, Lnp_D, tg_H, tg_D, tg, td_H, td_D, td, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Myy_neg, Mcrℓ_yy_neg, local_buckling_Myy_pos, Mcrℓ_yy_pos, distortional_buckling_P, Pcrd, distortional_buckling_Mxx, Mcrd)
    properties = UniStrutOut(section_inputs, geometry, gross_section_properties, net_section_properties, Lnp_H, Lnp_D, tg_H, tg_D, tg, td_H, td_D, td, local_buckling_P, Pcrℓ, local_buckling_Mxx, Mcrℓ_xx, local_buckling_Myy_neg, Mcrℓ_yy_neg, local_buckling_Myy_pos, Mcrℓ_yy_pos, distortional_buckling_P, Pcrd)



    return properties 

end



end #module