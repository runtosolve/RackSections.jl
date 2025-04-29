module Beams

using CrossSectionGeometry, SectionProperties, RMI, LinesCurvesNodes, CUFSM, AISIS100, LinearAlgebra

struct StepBeamInput

    H::Float64
    D::Float64
    W::Float64
    L::Float64
    R::Float64
    t::Float64
    E::Float64
    ν::Float64

end


struct StepBeam

    input::StepBeamInput
    geometry::@NamedTuple{coordinates::CrossSectionGeometry.ThinWalled, x::Vector{Float64}, y::Vector{Float64}}
    properties::SectionProperties.SectionPropertiesObject
    local_buckling_P::CUFSM.Model
    Pcrℓ::Float64
    local_buckling_Mxx_pos::CUFSM.Model
    Mcrℓ_xx_pos::Float64
    local_buckling_Mxx_neg::CUFSM.Model
    Mcrℓ_xx_neg::Float64
    local_buckling_Myy_neg::CUFSM.Model
    Mcrℓ_yy_neg::Float64
    local_buckling_Myy_pos::CUFSM.Model
    Mcrℓ_yy_pos::Float64 

end

struct AngledStepBeamInput

    H::Float64
    D::Float64
    W::Float64
    L::Float64
    A::Float64
    R::Float64
    t::Float64
    E::Float64
    ν::Float64

end


struct AngledStepBeam

    input::AngledStepBeamInput
    geometry::@NamedTuple{coordinates::CrossSectionGeometry.ThinWalled, x::Vector{Float64}, y::Vector{Float64}}
    properties::SectionProperties.SectionPropertiesObject
    local_buckling_P::CUFSM.Model
    Pcrℓ::Float64
    local_buckling_Mxx_pos::CUFSM.Model
    Mcrℓ_xx_pos::Float64
    local_buckling_Mxx_neg::CUFSM.Model
    Mcrℓ_xx_neg::Float64
    local_buckling_Myy_neg::CUFSM.Model
    Mcrℓ_yy_neg::Float64
    local_buckling_Myy_pos::CUFSM.Model
    Mcrℓ_yy_pos::Float64 

end


function step_beam_geometry(H, D, W, L, R, t)

    segments = [H-L, D-W, L, W, H, D]
    θ = [π/2, π, π/2,π, -π/2, 0.0]
    r = [R, R-t, R, R, R, R]
    n = [4, 4, 4, 4, 4, 4]
    # n_r = [3, 3, 3, 3, 3, 3, 3]
    n_r = [5, 5, 5, 5, 5, 5, 5]

    section_geometry = CrossSectionGeometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D, 0.0))

    x = [section_geometry.centerline_node_XY[i][1] for i in eachindex(section_geometry.centerline_node_XY)]
    y = [section_geometry.centerline_node_XY[i][2] for i in eachindex(section_geometry.centerline_node_XY)]

    geometry = (coordinates = section_geometry, x=x, y=y)

    return geometry

end




function step_beam(input)

    # input = StepBeamInput(H, D, W, L, R, t, E, ν)

     (; H, D, W, L, R, t, E, ν) = input

    geometry = step_beam_geometry(H, D, W, L, R, t)

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
    lengths = range(0.5*H, 1.25*H, 7)
    local_buckling_P = CUFSM.Tools.closed_section_analysis(geometry.x, geometry.y, fill(t, length(geometry.x)), lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs)
    eig = 1
    Pcrℓ = minimum(CUFSM.Tools.get_load_factor(local_buckling_P, eig))

    #local buckling, Mxx_pos 
    P = 0.0
    Mxx = 1.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.25*H, 1.0*H, 7)
    local_buckling_Mxx_pos = CUFSM.Tools.closed_section_analysis(geometry.x, geometry.y, fill(t, length(geometry.x)), lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs)
    Mcrℓ_xx_pos  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Mxx_pos, eig))

    #local buckling, Mxx_neg 
    P = 0.0
    Mxx = -1.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.25*H, 1.0*H, 7)
    local_buckling_Mxx_neg = CUFSM.Tools.closed_section_analysis(geometry.x, geometry.y, fill(t, length(geometry.x)), lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs)
    Mcrℓ_xx_neg  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Mxx_neg, eig))
    

    #local buckling, Myy_neg 
    P = 0.0
    Mxx = 0.0
    Myy = 1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.5*H, 1.0*H, 7)
    local_buckling_Myy_neg = CUFSM.Tools.closed_section_analysis(geometry.x, geometry.y, fill(t, length(geometry.x)), lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs)
    Mcrℓ_yy_neg  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Myy_neg, eig))

    #local buckling, Myy_pos
    P = 0.0
    Mxx = 0.0
    Myy = -1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.5*(H-L), 1.0*(H-L), 7)
    local_buckling_Myy_pos = CUFSM.Tools.closed_section_analysis(geometry.x, geometry.y, fill(t, length(geometry.x)), lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs)
    Mcrℓ_yy_pos  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Myy_pos, eig))

    #gather up everything 
    properties = StepBeam(input, geometry, gross_section_properties, local_buckling_P, Pcrℓ, local_buckling_Mxx_pos, Mcrℓ_xx_pos, local_buckling_Mxx_neg, Mcrℓ_xx_neg, local_buckling_Myy_neg, Mcrℓ_yy_neg, local_buckling_Myy_pos, Mcrℓ_yy_pos)

    return properties 

end




function angled_step_beam_geometry(H, D, W, L, A, R, t)

    A = deg2rad(A)

    angled_step_segment = L / sin(A)

    segments = [H-L, D-W+angled_step_segment*cos(A), angled_step_segment, W, H, D]
    θ = [π/2, π, A, π, -π/2, 0.0]
    r = [R, R-t, R, R, R, R]
    n = [4, 4, 4, 4, 4, 4]
    n_r = [5, 5, 5, 5, 5, 5, 5]

    section_geometry = CrossSectionGeometry.create_thin_walled_cross_section_geometry(segments, θ, n, r, n_r, t, centerline = "to left", offset = (D, 0.0))

    x = [section_geometry.centerline_node_XY[i][1] for i in eachindex(section_geometry.centerline_node_XY)]
    y = [section_geometry.centerline_node_XY[i][2] for i in eachindex(section_geometry.centerline_node_XY)]

    geometry = (coordinates = section_geometry, x=x, y=y)

    return geometry

end




function angled_step_beam(input)

     (; H, D, W, L, A, R, t, E, ν) = input

    # input = AngledStepBeamInput(H, D, W, L, A, R, t, E, ν)

    geometry = angled_step_beam_geometry(H, D, W, L, A, R, t)

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
    lengths = range(0.5*H, 1.25*H, 7)
    local_buckling_P = CUFSM.Tools.closed_section_analysis(geometry.x, geometry.y, fill(t, length(geometry.x)), lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs)
    eig = 1
    Pcrℓ = minimum(CUFSM.Tools.get_load_factor(local_buckling_P, eig))

    #local buckling, Mxx_pos 
    P = 0.0
    Mxx = 1.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.25*H, 1.0*H, 7)
    local_buckling_Mxx_pos = CUFSM.Tools.closed_section_analysis(geometry.x, geometry.y, fill(t, length(geometry.x)), lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs)
    Mcrℓ_xx_pos  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Mxx_pos, eig))

    #local buckling, Mxx_neg 
    P = 0.0
    Mxx = -1.0
    Myy = 0.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.25*H, 1.0*H, 7)
    local_buckling_Mxx_neg = CUFSM.Tools.closed_section_analysis(geometry.x, geometry.y, fill(t, length(geometry.x)), lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs)
    Mcrℓ_xx_neg  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Mxx_neg, eig))
    

    #local buckling, Myy_neg 
    P = 0.0
    Mxx = 0.0
    Myy = 1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.5*H, 1.0*H, 7)
    local_buckling_Myy_neg = CUFSM.Tools.closed_section_analysis(geometry.x, geometry.y, fill(t, length(geometry.x)), lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs)
    Mcrℓ_yy_neg  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Myy_neg, eig))

    #local buckling, Myy_pos
    P = 0.0
    Mxx = 0.0
    Myy = -1.0
    M11 = 0.0
    M22 = 0.0
    constraints = []
    springs = []

    lengths = range(0.5*(H-L), 1.0*(H-L), 7)
    local_buckling_Myy_pos = CUFSM.Tools.closed_section_analysis(geometry.x, geometry.y, fill(t, length(geometry.x)), lengths, E, ν, P, Mxx, Myy, M11, M22, constraints, springs)
    Mcrℓ_yy_pos  = minimum(CUFSM.Tools.get_load_factor(local_buckling_Myy_pos, eig))

    #gather up everything 
    properties = AngledStepBeam(input, geometry, gross_section_properties, local_buckling_P, Pcrℓ, local_buckling_Mxx_pos, Mcrℓ_xx_pos, local_buckling_Mxx_neg, Mcrℓ_xx_neg, local_buckling_Myy_neg, Mcrℓ_yy_neg, local_buckling_Myy_pos, Mcrℓ_yy_pos)

    return properties 

end

end #module