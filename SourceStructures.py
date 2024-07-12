from MathTools.PointGenerators.PeriodStruct.periodic_structure_cartesian_points_3d import periodic_structure_cartesian_points_3d
from MathTools.PlotVisualize.SpaceMapping.plot_dots_planes_3d import plot_dots_planes_3d


H1 = [2, 1, 1]
D1 = [0, 0, 0, 0, 0, 0]
T1 = [5, 5, 5]
Sh1 = [0, 0, 0]
P1 = [2, 2, 2]
PL1 = [[0, 0, 1, -1]]

H2 = [10, 4, 4]
D2 = [0.5, 0.5, 0, 0, 0, 0]
T2 = [2, 5, 7.5]
Sh2 = [0, 0, 0]
P2 = [10, 10, 10]
PL2 = [[0, 0, 1, -1]]

G1 = periodic_structure_cartesian_points_3d(H1, D1, T1, True, Sh1)
G2 = periodic_structure_cartesian_points_3d(H2, D2, T2, True, Sh2)

plot_dots_planes_3d(G1, PL1, P1, input_type=0, planes_key=True, key_save=False, key_show=True)
plot_dots_planes_3d(G2, PL2, P2, input_type=0, planes_key=True, key_save=False, key_show=True)
