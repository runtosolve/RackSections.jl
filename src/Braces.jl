module Braces

using CrossSectionGeometry, SectionProperties, RMI, LinesCurvesNodes, CUFSM, AISIS100, LinearAlgebra



struct CeeLipsBraceInput

    H::Float64
    D::Float64
    L::Float64
    R::Float64
    t::Float64
    E::Float64
    ν::Float64

end



struct CeeLipsBrace

    input::CeeLipsBraceInput
    geometry::@NamedTuple{coordinates::CrossSectionGeometry.ThinWalled, x::Vector{Float64}, y::Vector{Float64}}
    properties::SectionProperties.SectionPropertiesObject
    local_buckling_P::CUFSM.Model
    Pcrℓ::Float64
    distortional_buckling_P::CUFSM.Model
    Pcrd::Float64

end




struct CeeInput

    H::Float64
    D::Float64
    R::Float64
    t::Float64
    E::Float64
    ν::Float64

end

struct Cee

    input::CeeInput
    geometry::@NamedTuple{coordinates::CrossSectionGeometry.ThinWalled, x::Vector{Float64}, y::Vector{Float64}}
    properties::SectionProperties.SectionPropertiesObject
    local_buckling_P::CUFSM.Model
    Pcrℓ::Float64

end


struct AngleInput

    H::Float64
    D::Float64
    R::Float64
    t::Float64
    E::Float64
    ν::Float64

end

struct Angle

    input::AngleInput
    geometry::@NamedTuple{coordinates::CrossSectionGeometry.ThinWalled, x::Vector{Float64}, y::Vector{Float64}}
    properties::SectionProperties.SectionPropertiesObject
    local_buckling_P::CUFSM.Model
    Pcrℓ::Float64

end

struct PipeInput

    D::Float64
    t::Float64
    E::Float64
    ν::Float64

end

struct Pipe

    input::PipeInput
    geometry::@NamedTuple{coordinates::CrossSectionGeometry.ThinWalled, x::Vector{Float64}, y::Vector{Float64}}
    properties::SectionProperties.SectionPropertiesObject
    local_buckling_P::CUFSM.Model
    Pcrℓ::Float64

end


struct RectangularTubeBraceInput

    H::Float64
    D::Float64
    R::Float64
    t::Float64
    E::Float64
    ν::Float64

end

struct RectangularTubeBrace

    input::RectangularTubeBraceInput
    geometry::@NamedTuple{coordinates::CrossSectionGeometry.ThinWalled, x::Vector{Float64}, y::Vector{Float64}}
    properties::SectionProperties.SectionPropertiesObject
    local_buckling_P::CUFSM.Model
    Pcrℓ::Float64

end




function cee_geometry(H, D, R, t)

    segments = [D, H, D]
    θ = [π, -π/2, 0.0]
    r = [R, R]
    n = [4, 4, 4]
    n_r = [5, 5]

    section_geometry = CrossSectionGeometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D, H))

    x = [section_geometry.centerline_node_XY[i][1] for i in eachindex(section_geometry.centerline_node_XY)]
    y = [section_geometry.centerline_node_XY[i][2] for i in eachindex(section_geometry.centerline_node_XY)]

    geometry = (coordinates = section_geometry, x=x, y=y)

    return geometry

end




function cee(input)

     (; H, D, R, t, E, ν) = input

    # input = CeeInput(H, D, R, t, E, ν)

    geometry = cee_geometry(H, D, R, t)

    #gross section properties 
    gross_section_properties = SectionProperties.open_thin_walled(geometry.coordinates.centerline_node_XY, fill(t, length(geometry.x)-1)) 


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
    local_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, fill(t, length(geometry.x)-1), lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    eig = 1
    Pcrℓ = minimum(CUFSM.Tools.get_load_factor(local_buckling_P, eig))

    #gather up everything 
    properties = Cee(input, geometry, gross_section_properties, local_buckling_P, Pcrℓ)

    return properties 

end




function pipe_geometry(D, t, num_segments)

    # num_segments = 100
    # circumference = 2*π*D/2
    # segments = fill(circumference/num_segments, num_segments-1)
    # θ = range(0.0, 2π, num_segments)[1:end-1]
    # n = ones(Int, num_segments-1)


    # num_segments = 100
    θc = range(0, 2π, num_segments)
    x = D/2 * cos.(θc)
    y = D/2 * sin.(θc)
    

    segments = [norm([x[i+1]-x[i], y[i+1]-y[i]]) for i=1:length(x)-1] 
    θ = [atan((y[i+1]-y[i]), (x[i+1]-x[i])) for i=1:length(x)-1]
    n = ones(Int, num_segments-1)

    section_geometry = CrossSectionGeometry.create_thin_walled_cross_section_geometry(segments, θ, n, t, centerline = "to left", offset = (D, D/2))

    x = [section_geometry.centerline_node_XY[i][1] for i in eachindex(section_geometry.centerline_node_XY)]
    y = [section_geometry.centerline_node_XY[i][2] for i in eachindex(section_geometry.centerline_node_XY)]

    geometry = (coordinates = section_geometry, x=x, y=y)

    return geometry

end


function pipe(input)

     (; D, t, E, ν) = input

    # input = PipeInput(D, t, E, ν)

    geometry = pipe_geometry(D, t, 100)

    #gross section properties 
    gross_section_properties = SectionProperties.closed_thin_walled(geometry.coordinates.centerline_node_XY, fill(t, length(geometry.x))) 

    #remove NaNs to allow for writing to JSON
    gross_section_properties.xs = -1
    gross_section_properties.ys = -1
    gross_section_properties.Cw = -1
    gross_section_properties.B1 = -1
    gross_section_properties.B2 = -1
    gross_section_properties.wn = [-1, -1]
  

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
    lengths = range(2.0*D, 4.0*D, 7)
    local_buckling_P = CUFSM.Tools.closed_section_analysis(geometry.x, geometry.y, fill(t, length(geometry.x)), lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs)
    
    eig = 1
    Pcrℓ = minimum(CUFSM.Tools.get_load_factor(local_buckling_P, eig))

    #gather up everything 
    properties = Pipe(input, geometry, gross_section_properties, local_buckling_P, Pcrℓ)

    return properties 

end


function angle_geometry(H, D, R, t)

    segments = [H, D]
    θ = [-π/2, 0.0]
    r = [R]
    n = [4, 4]
    n_r = [5]

    section_geometry = CrossSectionGeometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (0.0, H))

    x = [section_geometry.centerline_node_XY[i][1] for i in eachindex(section_geometry.centerline_node_XY)]
    y = [section_geometry.centerline_node_XY[i][2] for i in eachindex(section_geometry.centerline_node_XY)]

    geometry = (coordinates = section_geometry, x=x, y=y)

    return geometry

end




function angle(input)

    (; H, D, R, t, E, ν) = input
    # input = AngleInput(H, D, R, t, E, ν)

    geometry = angle_geometry(H, D, R, t)

    #gross section properties 
    gross_section_properties = SectionProperties.open_thin_walled(geometry.coordinates.centerline_node_XY, fill(t, length(geometry.x)-1)) 


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
    lengths = range(3.0*minimum([H, D]), 15.0*minimum([H,D]), 13)
    local_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, fill(t, length(geometry.x)-1), lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    eig = 1
    Pcrℓ = minimum(CUFSM.Tools.get_load_factor(local_buckling_P, eig))

    #gather up everything 
    properties = Angle(input, geometry, gross_section_properties, local_buckling_P, Pcrℓ)

    return properties 

end





function cee_with_lips_brace_geometry(H, D, L, R, t)

    segments = [L, D, H, D, L]
    θ = [π/2, π, -π/2, 0.0, π/2]
    r = [R, R, R, R]
    n = [4, 4, 5, 4, 4]
    n_r = [5, 5, 5, 5]

    section_geometry = CrossSectionGeometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D, H-L))

    x = [section_geometry.centerline_node_XY[i][1] for i in eachindex(section_geometry.centerline_node_XY)]
    y = [section_geometry.centerline_node_XY[i][2] for i in eachindex(section_geometry.centerline_node_XY)]

    geometry = (coordinates = section_geometry, x=x, y=y)

    return geometry

end




function cee_with_lips_brace(inputs)

    #  H, D, L, R, t, E, ν, dh_H, dh_D, de_H, de_D, hole_pitch_H, hole_pitch_D, hole_length_H, hole_length_D = section_inputs

     (; H, D, L, R, t, E, ν) = inputs

    # section_inputs = CeeLipsBraceInput(H, D, L, R, t, E, ν)

    geometry = cee_with_lips_brace_geometry(H, D, L, R, t)

    #gross section properties 
    gross_section_properties = SectionProperties.open_thin_walled(geometry.coordinates.centerline_node_XY, fill(t, length(geometry.x)-1)) 

    num_elem = length(geometry.x) - 1
    tg = fill(t, num_elem)

    
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
    lengths = range(0.75*minimum([H, D]), 2.0*minimum([H,D]), 7)
    # lengths = round.(lengths, digits=2)  #like CUFSM v5.50 GUI 
    local_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, tg, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    #local_buckling_P = CUFSM.Tools.open_section_analysis(round.(geometry.x, digits=4), round.(geometry.y, digits=4), tg, round.(lengths, digits=2), E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    eig = 1
    Pcrℓ = minimum(CUFSM.Tools.get_load_factor(local_buckling_P, eig))

    
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


    td = fill(t, num_elem)


    #distortional buckling, compression 
    P = 1.0
    Mxx = 0.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.75*Lcrd, 1.5*Lcrd, 7)
    distortional_buckling_P = CUFSM.Tools.open_section_analysis(geometry.x, geometry.y, td, lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs, neigs)
    Pcrd  = minimum(CUFSM.Tools.get_load_factor(distortional_buckling_P, eig))


    #gather up everything 
    properties = CeeLipsBrace(inputs, geometry, gross_section_properties, local_buckling_P, Pcrℓ, distortional_buckling_P, Pcrd)

    return properties 

end




function rectangular_tube_brace_geometry(H, D, R, t)


    segments = [H, D, H, D]
    θ = [π/2, π, -π/2, 0.0]
    r = [R, R, R, R]
    n = [4, 4, 5, 4]
    n_r = [5, 5, 5, 5]

    section_geometry = CrossSectionGeometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D, 0.0))

    x = [section_geometry.centerline_node_XY[i][1] for i in eachindex(section_geometry.centerline_node_XY)]
    y = [section_geometry.centerline_node_XY[i][2] for i in eachindex(section_geometry.centerline_node_XY)]

    geometry = (coordinates = section_geometry, x=x, y=y)

    return geometry

end



function rectangular_tube_brace(section_inputs)

     (; H, D, R, t, E, ν) = section_inputs


    geometry = rectangular_tube_brace_geometry(H, D, R, t)

    #gross section properties 
    gross_section_properties = SectionProperties.closed_thin_walled(geometry.coordinates.centerline_node_XY, fill(t, length(geometry.x))) 

    #remove NaNs to allow for writing to JSON
    gross_section_properties.xs = -1
    gross_section_properties.ys = -1
    gross_section_properties.Cw = -1
    gross_section_properties.B1 = -1
    gross_section_properties.B2 = -1
    gross_section_properties.wn = [-1, -1]

    

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
    local_buckling_P = CUFSM.Tools.closed_section_analysis(geometry.x, geometry.y, fill(t, length(geometry.x)), lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs)
    eig = 1
    Pcrℓ = minimum(CUFSM.Tools.get_load_factor(local_buckling_P, eig))

    
    #gather up everything 
    properties = RectangularTubeBrace(section_inputs, geometry, gross_section_properties, local_buckling_P, Pcrℓ)

    return properties 

end






end #module

