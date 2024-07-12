# Electrostatics of point charges

One of the fundamental problems in **computational electrostatics** involves computing and visualizing **electrostatic
fields** generated by **arbitrary configurations** of **point charges** in space. Due to a variety and complexity of
input parameters and requirements, computations can be highly encumbered.

This repository offers a **quick exposition** of **key abilities** of my
main [**MathTools**](https://github.com/StDLabs/MathTools) library to effectively generate and
visualize described electrostatic fields due to discrete configurations of point charges.

## Spatial structures and crystal lattices

There are several obstacles to overcome, but the first will involve the complexity of the point charges
**configuration**, which can be a **multipole**, **periodic spatial lattice** or a more complex
**functional rule-based configuration**.

I would like to show two structures that differ in complexity: a simple **electric dipole** and an interplanar area of
a **perfect crystalline lattice** with **triclinic symmetry** (parallelepipedal cells). The triclinic crystal system is
the most general crystal system, characterized by three axes of unequal lengths that intersect at oblique angles,
resulting in a unit cell shaped like a parallelepiped.

<p align="center">
    <img width="45%" src="https://github.com/StDLabs/DiscreteElecStat/blob/main/Content/Dipole2.png" alt="dipole"/>
&nbsp;
    <img width="45%" src="https://github.com/StDLabs/DiscreteElecStat/blob/main/Content/Lattice1.gif" alt="lattice"/>
</p>

For every case I used a [**three-dimensional point generator**](https://github.com/StDLabs/MathTools/blob/main/PointGenerators/PeriodStruct/periodic_structure_cartesian_points_3d.py) 
to get an array of point coordinates, and a [**universal three-dimensional plotter**](https://github.com/StDLabs/MathTools/blob/main/PlotVisualize/SpaceMapping/plot_dots_planes_3d.py) 
for points. You could see precise input data in a simple [**script file**](https://github.com/StDLabs/DiscreteElecStat/blob/main/SourceStructures.py).

## Field calculations, grid settings, input parameters

The complexity of possible input structures requires **universal tools** that can generate corresponding scalar and
vector fields using a variety of **spatial grids** and **sections**. 

First, we need a [function](https://github.com/StDLabs/MathTools/blob/main/PointGenerators/SpatialFields/Central/point_charges_field_calc.py)
that can calculate both scalar and vector fields for a specific observation point $\vec{r}$ according to **Coulomb’s law**
and **superposition principle**. The function will take our arbitrary point configuration $G=\{\vec{r}_n\}$, corresponding
set of charges $Q=\{q_n\}$ and will compute corresponding field values in every nodal points:

$$\varphi (\vec{r})=k\sum_{n} \frac{q_n}{|\vec{r}-\vec{r}_n|}$$
$$\vec{E} (\vec{r})=k\sum_{n} \frac{q_n(\vec{r}-\vec{r}_n)}{|\vec{r}-\vec{r}_n|^3}$$
