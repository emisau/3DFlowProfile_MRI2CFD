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


def adimention_MRI_voxel_grid(RowNb, ColumnNb, Time, Xc, Yc, VectXcoord, VectYcoord, R_voxel_max, Radial_Coordinate, VelocityMat, VelocityMat_Parabolic):

  Angular_Coordinate = np.zeros((RowNb,ColumnNb,len(Time)))
  Recal_Xcoord = np.zeros((RowNb,ColumnNb,len(Time)))
  Recal_Ycoord = np.zeros((RowNb,ColumnNb,len(Time)))

  # See coordinates transformation https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations
  for k in range(len(Time)):
    for i in range(RowNb):
        dY = (VectYcoord[i] - Yc[k])
        for j in range(ColumnNb):
          dX = (VectXcoord[j] - Xc[k])
          r = np.hypot(dX,dY)  ### Equivalent to:  r = np.sqrt(X*X + Y*Y)
          if (dX >=0.0 ) and (dY >= 0.0):
            Angular_Coordinate[i,j,k] = np.arctan(dY/dX)
          elif (dX < 0.0) and (dY >= 0.0):
            Angular_Coordinate[i,j,k] = math.pi + np.arctan(dY/dX)
          elif (dX < 0.0) and (dY < 0.0):
            Angular_Coordinate[i,j,k] = math.pi + np.arctan(dY/dX)
          elif (dX >=0.0) and (dY < 0.0):
            Angular_Coordinate[i,j,k] = 2*math.pi + np.arctan(dY/dX)
            
          Recal_Xcoord[i,j,k] = Radial_Coordinate[i,j,k]/ R_voxel_max[k] * np.cos(Angular_Coordinate[i,j,k])
          Recal_Ycoord[i,j,k] = Radial_Coordinate[i,j,k]/ R_voxel_max[k] * np.sin(Angular_Coordinate[i,j,k]) 


  # Find global (over all times) maximum and minimum of velocity
  vel_min = VelocityMat[0,0,0]
  vel_max = VelocityMat[0,0,0]

  for k in range(len(Time)):
    for i in range(RowNb):
      for j in range(ColumnNb):
          vel_min = min(vel_min, VelocityMat[i,j,k])
          vel_max = max(vel_max, VelocityMat[i,j,k])
          vel_min = min(vel_min, VelocityMat_Parabolic[i,j,k])
          vel_max = max(vel_max, VelocityMat_Parabolic[i,j,k])

  return Angular_Coordinate, Recal_Xcoord, Recal_Ycoord, vel_min, vel_max
  
