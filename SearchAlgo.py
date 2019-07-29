import numpy as np
import math

import Gmsh

# ********************************************************************************************************************************  

def closest_value(NNodes, Time, R_inlet_max, R_voxel_max, r_inlet, Phi_inlet, RowNb, ColumnNb, Radial_Coordinate, BoundaryNodes, Angular_Coordinate, VelocityMat):

  #################################################
  ### Computation of closest value of an inlet point
  ### to a voxel-grid point
  #################################################


  temp_dR = 0
  i_idx = 0
  j_idx = 0
  Velocity_inlet = np.zeros((NNodes,len(Time)))


  for k in range(len(Time)):
    
    print('k = {}/{}'.format(k,len(Time)))

    for i_node in range (NNodes): # pick one node of the inlet
      X_loc = r_inlet[i_node] / R_inlet_max * np.cos(Phi_inlet[i_node])
      Y_loc = r_inlet[i_node] / R_inlet_max * np.sin(Phi_inlet[i_node])

      
      MinDist= 10000.0 * R_inlet_max
      
      for i_row in range(RowNb): 
        for j_col in range(ColumnNb):

          #X_MRI = (VectXcoord[j_col] - Xc[k])
          #Y_MRI = (VectYcoord[i_row] - Yc[k])
          
          X_MRI = Radial_Coordinate[i_row,j_col,k]/ R_voxel_max[k] * np.cos(Angular_Coordinate[i_row,j_col,k])
          Y_MRI = Radial_Coordinate[i_row,j_col,k]/ R_voxel_max[k] * np.sin(Angular_Coordinate[i_row,j_col,k])



          #temp_dR = (VectXcoord[j_col] - X_local[i_node]) * (VectXcoord[j_col] - X_local[i_node]) + (VectYcoord[i_row] - Y_local[i_node]) * (VectYcoord[i_row] - Y_local[i_node])
          temp_dR = (X_MRI - X_loc) * (X_MRI - X_loc) + (Y_MRI - Y_loc) * (Y_MRI - Y_loc)

          if (temp_dR < MinDist):
            i_idx = i_row
            j_idx = j_col
            MinDist = temp_dR
          
      #print('MinDist = {}'.format(MinDist))
      Velocity_inlet[i_node,k] = VelocityMat[i_idx,j_idx,k]
      if (BoundaryNodes[i_node] == True):
        Velocity_inlet[i_node,k] = 0

  '''
    print("VelocityMat:")
    print(VelocityMat[:,:,k])
    print("Velocity_inlet:")
    print(Velocity_inlet[:,k])
    print("**********************************************************************************************")
  '''


  return Velocity_inlet
