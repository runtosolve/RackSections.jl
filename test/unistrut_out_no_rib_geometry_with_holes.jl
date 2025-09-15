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



geometry = RackSections.Columns.unistrut_out_no_rib_geometry(H, D, L1, L2, R, t, dh_H, dh_D, de_H, de_D)

