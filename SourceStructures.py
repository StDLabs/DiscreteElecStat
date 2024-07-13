from MathTools.PointGenerators.PeriodStruct.periodic_structure_cartesian_points_3d import periodic_structure_cartesian_points_3d
from MathTools.PlotVisualize.SpaceMapping.plot_dots_planes_3d import plot_dots_planes_3d
from MathTools.PointGenerators.SpatialFields.field_section import field_section
from MathTools.PlotVisualize.PlaneMapping.plot_scalar_field_contour_map import plot_scalar_field_contour_map
from MathTools.PlotVisualize.SpaceMapping.plot_vector_field_section import plot_vector_field_section
from MathTools.PointGenerators.SpatialFields.field_3d import field_3d
from MathTools.PlotVisualize.SpaceMapping.plot_scalar_filed_contour_3d import plot_scalar_filed_contour_3d
from MathTools.PlotVisualize.SpaceMapping.plot_field_lines_3D import plot_field_lines_3D
import numpy as np

# structure parameters for dipole and crystal lattice:

H1 = [2, 1, 1]
D1 = [0, 0, 0, 0, 0, 0]
T1 = [5, 5, 5]
Sh1 = [0, 0, 0]

H2 = [10, 4, 4]
D2 = [0.5, 0.5, 0, 0, 0, 0]
T2 = [2, 5, 7.5]
Sh2 = [0, 0, 0]

# areas for calculations and plane sections
P1 = [5, 5, 5]
P2 = [10, 10, 10]
PL1 = [[0, 0, 1, -1]]
PL2 = [[0, 0, 1, -3.0]]

# getting arrays of point coordinates with three-dimensional point generators
G1 = periodic_structure_cartesian_points_3d(H1, D1, T1, True, Sh1)
G2 = periodic_structure_cartesian_points_3d(H2, D2, T2, True, Sh2)

plot_dots_planes_3d(G1, PL1, P1, input_type=0, planes_key=True, key_save=False, key_show=True)
plot_dots_planes_3d(G2, PL2, P2, input_type=0, planes_key=True, key_save=False, key_show=True)

# plane grid resolutions
Npl1 = [50, 50]
Npl2 = [50, 50]

# charge values
Q1 = [1, -1]
Q2 = [1 for i in range(0, len(G2))]

# function parameters
# A = [k, n, Rm, Fm]. (k = 1)
# Fi(r) = k*Qi/(r^n) if r > Rm; Fi(r) = Fm if r <= Rm
A_phi = [1, 1, 0.5, 0]  # (k = 1, n = 1)
A_E = [1, 2, 0.5, 0]  # (k = 1, n = 2)
function = 'hyperbola'
option_phi = 'scalar'
option_E = 'vector'

# calculation of scalar potential field values on grids
M1x_phi, M1y_phi, F1_phi = field_section(PL1[0], P1, Npl1, G1, Q1, A_phi, function, option_phi)
M2x_phi, M2y_phi, F2_phi = field_section(PL2[0], P2, Npl2, G2, Q2, A_phi, function, option_phi)

# calculation of vector field intensity values on grids
M1x_E, M1y_E, F1_E = field_section(PL1[0], P1, Npl1, G1, Q1, A_E, function, option_E)
M2x_E, M2y_E, F2_E = field_section(PL2[0], P2, Npl2, G2, Q2, A_E, function, option_E)
F1_E = np.transpose(np.array([F1_E]), (3, 2, 1, 0)).tolist()
F2_E = np.transpose(np.array([F2_E]), (3, 2, 1, 0)).tolist()

# plot contour maps
LN11 = 50
LN12 = 20
LN21 = 50
LN22 = 20
Title11 = 'Contour map, dipole scalar field'
Title12 = 'Contour map, dipole scalar field'
Title21 = 'Contour map, lattice scalar field'
Title22 = 'Contour map, lattice scalar field'
layout1 = 'DarkBlue'
layout2 = 'White'
plot_scalar_field_contour_map(M1x_phi, M1y_phi, F1_phi, LN11, Title11, PL1[0], P1, False, layout1)
plot_scalar_field_contour_map(M1x_phi, M1y_phi, F1_phi, LN12, Title12, PL1[0], P1, False, layout2)
plot_scalar_field_contour_map(M2x_phi, M2y_phi, F2_phi, LN21, Title21, PL2[0], P2, False, layout1)
plot_scalar_field_contour_map(M2x_phi, M2y_phi, F2_phi, LN22, Title22, PL2[0], P2, False, layout2)

# plot vector field sections
title_main1 = 'dipole vector field section'
title_bar1 = ' '
title_main2 = 'lattice vector field section'
title_bar2 = ' '
plot_vector_field_section(P1, PL1[0], F1_E, title_main1, title_bar1)
plot_vector_field_section(P2, PL2[0], F2_E, title_main2, title_bar2)

# field generators on 3D grids
N3d1 = [30, 30, 30]
N3d2 = [30, 30, 30]
M11, MQ11 = field_3d(P1, N3d1, G1, Q1, A_phi, function, option_phi)
M21, MQ21 = field_3d(P2, N3d2, G2, Q2, A_phi, function, option_phi)
M12, MQ12 = field_3d(P1, N3d1, G1, Q1, A_E, function, option_E)
M22, MQ22 = field_3d(P2, N3d2, G2, Q2, A_E, function, option_E)

# plot isosurfaces
LN1 = 20
LN2 = 20
plot_scalar_filed_contour_3d(P1, MQ11, LN1, title_main1, title_bar1)
plot_scalar_filed_contour_3d(P2, MQ21, LN2, title_main2, title_bar2)

# plot field lines
seed_resolution = 15
seed_visible = True
MQ12 = np.transpose(np.array(MQ12), (3, 0, 1, 2)).tolist()
MQ22 = np.transpose(np.array(MQ22), (3, 0, 1, 2)).tolist()
plot_field_lines_3D(P1, MQ12, seed_resolution, seed_visible, title_main1, title_bar1)
plot_field_lines_3D(P2, MQ22, seed_resolution, seed_visible, title_main1, title_bar1)
