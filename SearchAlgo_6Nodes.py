##############################################################################################################
##############################################################################################################
########                                                                                              ########
######## Copyright (C) Emilie Sauvage 2017 - e.sauvage@ucl.ac.uk / sauvage.emilie@gmail.com           ######## 
######## All rights reserved.                                                                         ########
########                                                                                              ########
######## If you are using this software, please let me know in case your work is made public or       ########
######## leads to a publication.                                                                      ########
########                                                                                              ########
##############################################################################################################
##############################################################################################################

import numpy as np
import math

import Gmsh

# ********************************************************************************************************************************  

def closest_value_6Nodes(NNodes, Time, R_inlet_max, R_voxel_max, r_inlet, Phi_inlet, RowNb, ColumnNb, Radial_Coordinate, BoundaryNodes, Angular_Coordinate, VelocityMat):

  #################################################
  ### Computation of closest value of an inlet point
  ### to a voxel-grid point
  #################################################
  
  Velocity_inlet = np.zeros((NNodes,len(Time)))


  for k in range(len(Time)):
    
    print('k = {}/{}'.format(k,len(Time)))

    for i_node in range (NNodes): # pick one node of the inlet
      X_loc = r_inlet[i_node] / R_inlet_max * np.cos(Phi_inlet[i_node])
      Y_loc = r_inlet[i_node] / R_inlet_max * np.sin(Phi_inlet[i_node])

      temp_dR = 0
      i_idx = 0
      j_idx = 0
      i_idx_1st = 0
      j_idx_1st = 0
      i_idx_2nd = 0
      j_idx_2nd = 0
      i_idx_3rd = 0
      j_idx_3rd = 0
      i_idx_4th = 0
      j_idx_4th = 0
      i_idx_5th = 0
      j_idx_5th = 0
      i_idx_6th = 0
      j_idx_6th = 0
      
      MinDist_1st= 10000.0 * R_inlet_max
      MinDist_2nd= 10000.0 * R_inlet_max
      MinDist_3rd= 10000.0 * R_inlet_max
      MinDist_4th= 10000.0 * R_inlet_max
      MinDist_5th= 10000.0 * R_inlet_max
      MinDist_6th= 10000.0 * R_inlet_max

      for i_row in range(RowNb): 
        for j_col in range(ColumnNb):

          #X_MRI = (VectXcoord[j_col] - Xc[k])
          #Y_MRI = (VectYcoord[i_row] - Yc[k])
          
          X_MRI = Radial_Coordinate[i_row,j_col,k]/ R_voxel_max[k] * np.cos(Angular_Coordinate[i_row,j_col,k])
          Y_MRI = Radial_Coordinate[i_row,j_col,k]/ R_voxel_max[k] * np.sin(Angular_Coordinate[i_row,j_col,k])

          temp_dR = (X_MRI - X_loc) * (X_MRI - X_loc) + (Y_MRI - Y_loc) * (Y_MRI - Y_loc)

          if (temp_dR < MinDist_1st):
            i_idx_1st = i_row
            j_idx_1st = j_col
            MinDist_1st = temp_dR

      for i_row in range(RowNb): 
        for j_col in range(ColumnNb):

          X_MRI = Radial_Coordinate[i_row,j_col,k]/ R_voxel_max[k] * np.cos(Angular_Coordinate[i_row,j_col,k])
          Y_MRI = Radial_Coordinate[i_row,j_col,k]/ R_voxel_max[k] * np.sin(Angular_Coordinate[i_row,j_col,k])

          temp_dR = (X_MRI - X_loc) * (X_MRI - X_loc) + (Y_MRI - Y_loc) * (Y_MRI - Y_loc)

          if (i_row != i_idx_1st) and (j_col !=j_idx_1st) and  (temp_dR < MinDist_2nd): #Look for the 2nd nearest point
            i_idx_2nd = i_row
            j_idx_2nd = j_col
            MinDist_2nd = temp_dR
            
      for i_row in range(RowNb): 
        for j_col in range(ColumnNb):

          X_MRI = Radial_Coordinate[i_row,j_col,k]/ R_voxel_max[k] * np.cos(Angular_Coordinate[i_row,j_col,k])
          Y_MRI = Radial_Coordinate[i_row,j_col,k]/ R_voxel_max[k] * np.sin(Angular_Coordinate[i_row,j_col,k])

          temp_dR = (X_MRI - X_loc) * (X_MRI - X_loc) + (Y_MRI - Y_loc) * (Y_MRI - Y_loc)

          if (i_row != i_idx_1st) and (i_row != i_idx_2nd) and (j_col !=j_idx_1st) and (j_col !=j_idx_2nd) and  (temp_dR < MinDist_3rd): #Look for the 3rd nearest point
            i_idx_3rd = i_row
            j_idx_3rd = j_col
            MinDist_3rd = temp_dR
            
      for i_row in range(RowNb): 
        for j_col in range(ColumnNb):

          X_MRI = Radial_Coordinate[i_row,j_col,k]/ R_voxel_max[k] * np.cos(Angular_Coordinate[i_row,j_col,k])
          Y_MRI = Radial_Coordinate[i_row,j_col,k]/ R_voxel_max[k] * np.sin(Angular_Coordinate[i_row,j_col,k])

          temp_dR = (X_MRI - X_loc) * (X_MRI - X_loc) + (Y_MRI - Y_loc) * (Y_MRI - Y_loc)

          if (i_row != i_idx_1st) and (j_col!= j_idx_1st) and (i_row !=i_idx_2nd) and (j_col !=j_idx_2nd) and (i_row !=i_idx_3rd) and (j_col !=j_idx_3rd) and  (temp_dR < MinDist_4th): #Look for the 4th nearest point
            i_idx_4th = i_row
            j_idx_4th = j_col
            MinDist_4th = temp_dR
            
      for i_row in range(RowNb): 
        for j_col in range(ColumnNb):

          X_MRI = Radial_Coordinate[i_row,j_col,k]/ R_voxel_max[k] * np.cos(Angular_Coordinate[i_row,j_col,k])
          Y_MRI = Radial_Coordinate[i_row,j_col,k]/ R_voxel_max[k] * np.sin(Angular_Coordinate[i_row,j_col,k])

          temp_dR = (X_MRI - X_loc) * (X_MRI - X_loc) + (Y_MRI - Y_loc) * (Y_MRI - Y_loc)

          if (i_row != i_idx_1st) and (j_col!= j_idx_1st) and (i_row !=i_idx_2nd) and (j_col !=j_idx_2nd) and (i_row !=i_idx_3rd) and (j_col !=j_idx_3rd) and (i_row !=i_idx_4th) and (j_col !=j_idx_4th) and  (temp_dR < MinDist_5th): #Look for the 5th nearest point
            i_idx_5th = i_row
            j_idx_5th = j_col
            MinDist_5th = temp_dR
            
      for i_row in range(RowNb): 
        for j_col in range(ColumnNb):

          X_MRI = Radial_Coordinate[i_row,j_col,k]/ R_voxel_max[k] * np.cos(Angular_Coordinate[i_row,j_col,k])
          Y_MRI = Radial_Coordinate[i_row,j_col,k]/ R_voxel_max[k] * np.sin(Angular_Coordinate[i_row,j_col,k])

          temp_dR = (X_MRI - X_loc) * (X_MRI - X_loc) + (Y_MRI - Y_loc) * (Y_MRI - Y_loc)

          if (i_row != i_idx_1st) and (j_col!= j_idx_1st) and (i_row !=i_idx_2nd) and (j_col !=j_idx_2nd) and (i_row !=i_idx_3rd) and (j_col !=j_idx_3rd) and (i_row !=i_idx_4th) and (j_col !=j_idx_4th) and (i_row !=i_idx_5th) and (j_col !=j_idx_5th) and  (temp_dR < MinDist_6th): #Look for the 3rd nearest point
            i_idx_6th = i_row
            j_idx_6th = j_col
            MinDist_6th = temp_dR
            

      w1 = 1.0 / (MinDist_1st + 1.e-4) # The constant 1.e-4 is here so that we don't divide by zero
      w2 = 1.0 / (MinDist_2nd + 1.e-4)
      w3 = 1.0 / (MinDist_3rd + 1.e-4)
      w4 = 1.0 / (MinDist_4th + 1.e-4) # The constant 1.e-4 is here so that we don't divide by zero
      w5 = 1.0 / (MinDist_5th + 1.e-4)
      w6 = 1.0 / (MinDist_6th + 1.e-4)
      
      SumW = w1 + w2 + w3 + w4 + w5 + w6
      Velocity_inlet[i_node,k] = ( VelocityMat[i_idx_1st,j_idx_1st,k] * w1 + \
                                   VelocityMat[i_idx_2nd,j_idx_2nd,k] * w2 + \
                                   VelocityMat[i_idx_3rd,j_idx_3rd,k] * w3 + \
                                   VelocityMat[i_idx_4th,j_idx_4th,k] * w4 + \
                                   VelocityMat[i_idx_5th,j_idx_5th,k] * w5 + \
                                   VelocityMat[i_idx_6th,j_idx_6th,k] * w6 ) / SumW
      
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
