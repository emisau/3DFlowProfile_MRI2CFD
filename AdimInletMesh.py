##############################################################################################################
##############################################################################################################
########                                                                                              ########
######## Copyright (C) Emilie Sauvage 2017 - e.sauvage@ucl.ac.uk / sauvage.emilie@gmail.com           ######## 
######## All rights reserved.                                                                         ########
########                                                                                              ########
######## If you are using this software, please let me know in case your work leads to a publication. ########
########                                                                                              ########
##############################################################################################################
##############################################################################################################
import numpy as np
import math

import Gmsh

# ********************************************************************************************************************************  

def adimention_Inlet_grid(meshfile, RowNb, ColumnNb, Time, VectXcoord, VectYcoord, VelocityMat, NormalVector):

  #################################################
  ### coordinates transformation of the points
  ### at the inlet into cylindrical coordinates
  #################################################
  R = 0
  Xc_inlet = 0
  Yc_inlet = 0
  Zc_inlet = 0
  Sum_x = 0
  Sum_y = 0
  Sum_z = 0
  Sum_R = 0
  dX = 0
  dY = 0
  dZ = 0

  x,y,z, lines, triags, tetras = Gmsh.read_gmsh_file(meshfile)

  ### Compute the center of the inlet surface 
  NNodes = len(x)
  Xc_inlet, Yc_inlet, Zc_inlet = Gmsh.surface_gravity_center(triags, x, y, z)
  print('The new gravity center = {}, {}, {}'.format(Xc_inlet, Yc_inlet, Zc_inlet))


  #################################################
  ### Define the local coordinate system
  ### on the surface of the inlet
  #################################################

  ##Vector of the base -  taken from VMTK (I think...)
  X_vect = [0.011,0,0]
  Y_vect = [0,0.011,0]
  N_vect = [0.6843449, -0.5662259, 0.4594129]


  ## Change of coordinates:
  r_inlet = np.zeros(NNodes)
  Phi_inlet = np.zeros(NNodes)
  X_local = np.zeros(NNodes)
  Y_local = np.zeros(NNodes)

  # Transformation from Cartesian to cylindrical coordinates (x,y) --> (r, phi)
  # See coordinates transformation https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations, https://en.wikipedia.org/wiki/Cylindrical_coordinate_system

# Removed the Z-component

  R_inlet_max = 0.0
  for i_node in range (NNodes):
    dX = x[i_node] - Xc_inlet
    dY = y[i_node] - Yc_inlet
    dZ = z[i_node] - Zc_inlet
    r_inlet[i_node] = np.sqrt(dX*dX + dY*dY) #This is the magnitude of "r" - the norm of "r"
    if (r_inlet[i_node] > R_inlet_max):
      R_inlet_max = r_inlet[i_node]
    r_vect_comp = [dX,dY] #May be the computation of r_vect_comp is useless
    dot_product_rX = (dX*X_vect[0] + dY*X_vect[1])
    dot_product_rY = (dX*Y_vect[0] + dY*Y_vect[1])
    ##### TO CHECK:  Mara has removed the Z-components from these expressions above
    
    
    arccos_argument = dot_product_rX / (r_inlet[i_node] * np.linalg.norm(X_vect))
    
    if ( dot_product_rY  >= 0.0):
      Phi_inlet[i_node] = np.arccos(arccos_argument)
    else:
      Phi_inlet[i_node] = 2.0 * math.pi - np.arccos(arccos_argument)

  #print ('BEFORE ADIM: the r-vector at the inlet looks like: {}'.format(r_inlet))

  print ('the averaged-radius of the inlet is: {}'.format(R_inlet_max))

  for i_node in range (NNodes):
    X_local[i_node] = r_inlet[i_node] * np.cos(Phi_inlet[i_node])   #Parametric coord system --> u
    Y_local[i_node] = r_inlet[i_node] * np.sin(Phi_inlet[i_node])   #Parametric coord system --> v

  return R_inlet_max, NNodes, X_local, Y_local, r_inlet, Phi_inlet
