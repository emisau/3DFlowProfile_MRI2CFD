# 3DFlowProfile_MRI2CFD
Python tool to extract flow profile from phase contrast MRI images and convert it to inlet condition for fluid domain in a computer simulation.

The main idea of the code is to transfer the flow profile information from the voxel grid to the inlet surface triangulation by overlaying them in cylindrical coordinates. A cylindrical coordinates transformation is applied on both grids. Subsequently, a closest node search algorithm is used to achieve the correspondance between the grids. 
