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
import csv
import os

# My own module called 'Gmsh' (functions in file Gmsh.py)
import Gmsh
import VoxelGrid
import Parabole
import AdimVoxGrid
import AdimInletMesh
import SearchAlgo
import SearchAlgo_3Nodes
import SearchAlgo_6Nodes
import Plot
import WriteFiles
import WritePoints


from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.tri as tri

from subprocess import call

def plot_to_file_2profiles(X,Y,Vel1,Vel2,z_min, z_max, filename):

  fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True, subplot_kw={'projection':'3d'}, 
                                 figsize=plt.figaspect(0.5)*1.5)
  fig.tight_layout()


  Xtmp, Ytmp = np.meshgrid(X, Y)
  surf = ax1.plot_surface(Xtmp, Ytmp, Vel1, rstride=1, cstride=1, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)
  surf = ax2.plot_surface(Xtmp, Ytmp, Vel2, rstride=1, cstride=1, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)
  ax1.set_zlim(z_min, z_max)
  ax2.set_zlim(z_min, z_max)

  ax1.zaxis.set_major_locator(LinearLocator(10))
  ax1.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
  ax2.zaxis.set_major_locator(LinearLocator(10))
  ax2.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
  
  fig.colorbar(surf, shrink=0.5, aspect=5)



  #plt.show()
  fig.savefig(filename)
  plt.close()
  
# ****************************************************************  

def plot_to_file_3profiles(X1,Y1,X2,Y2,Vel1,Vel2,Vel3,z_min, z_max, filename):

  fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=False, sharey=False, subplot_kw={'projection':'3d'}, 
                                 figsize=plt.figaspect(0.5)*1.5)
  fig.tight_layout()


  Xtmp1, Ytmp1 = np.meshgrid(X1, Y1)
  surf = ax1.plot_surface(Xtmp1, Ytmp1, Vel1, rstride=1, cstride=1, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)
  surf = ax2.plot_surface(Xtmp1, Ytmp1, Vel2, rstride=1, cstride=1, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)
  ax1.set_zlim(z_min, z_max)
  ax2.set_zlim(z_min, z_max)


  triang = tri.Triangulation(X2, Y2)
  
  #-----------------------------------------------------------------------------
  # Refine data
  #-----------------------------------------------------------------------------
  refiner = tri.UniformTriRefiner(triang)
  tri_refi, Vel3_refi = refiner.refine_field(Vel3, subdiv=3)
  

  #print(X2.shape)
  #print(Y2.shape)
  #print(Vel3.shape)
  
  #-----------------------------------------------------------------------------
  # Plot the triangulation and the high-res iso-contours
  #-----------------------------------------------------------------------------
  #plt.figure()
  #plt.gca().set_aspect('equal')
  ax3.triplot(triang, lw=0.5, color='black')
  ax3.set_zlim(z_min, z_max)

  levels = np.arange(z_min, z_max, 0.025)
  #cmap = cm.get_cmap(name='terrain', lut=None)
  ax3.plot_trisurf(X2, Y2, Vel3, cmap=cm.coolwarm, linewidth=0.2)
  ax3.set_zlim(z_min, z_max)

  levels = np.arange(z_min, z_max, 0.025)
  cmap = cm.get_cmap(name='terrain', lut=None)
  #plt.tricontourf(tri_refi, Vel3_refi.flatten(), levels=levels, cmap=cmap)
  ax3.plot_trisurf(X2, Y2, Vel3, cmap=cm.jet, linewidth=0.2)
  #print(Vel3)
  #plt.tricontourf(tri_refi, Vel3_refi, levels=levels, cmap=cmap)
  #plt.tricontour(tri_refi, z_test_refi, levels=levels,
                #colors=['0.25', '0.5', '0.5', '0.5', '0.5'],
                #linewidths=[1.0, 0.5, 0.5, 0.5, 0.5])
  

  ax1.zaxis.set_major_locator(LinearLocator(10))
  ax1.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
  ax2.zaxis.set_major_locator(LinearLocator(10))
  ax2.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
  ax3.zaxis.set_major_locator(LinearLocator(10))
  ax3.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
  
  fig.colorbar(surf, shrink=0.5, aspect=5)



  #plt.show()
  fig.savefig(filename)
  plt.close()
# ****************************************************************  


if __name__ == '__main__':


  
  MRI_Inputfile =  "TOF3_3DProfileData.csv"
  NormalVector = [0.6843449, -0.5662259, 0.4594129]  ##This is the normal vector to the inlet surface obtained from VMTK
  BaseFolder_path = "/home/emilie/Work/MyLittleCodes/Python/3DFlowProfile_MRI2CFD_MariaBoumpouli/TOF3_m"



  ResolX, ResolY, RowNb, ColumnNb, VectXcoord, VectYcoord, Time, VelocityMat = VoxelGrid.read_MRI_voxel_grid(MRI_Inputfile)

  print("Extracting ellipse areas ...")
  XR, XL, YT, YB, axis1, axis2, elipse_area, R_voxel_max = Parabole.extract_ellipse_area(ResolX, ResolY, RowNb, ColumnNb, Time, VectXcoord, VectYcoord, VelocityMat)

  print("Computing equivalent parabolic profile ...")
  Xc, Yc, MeanVel_perPhase, VelocityMat_Parabolic, Radial_Coordinate = Parabole.compute_equivalent_parabolic(ResolX, ResolY, RowNb, ColumnNb, Time, VectXcoord, VectYcoord, XR, XL, YT, YB, axis1, axis2, elipse_area, VelocityMat)
  
  print("Adimensionalisation of voxel grid ...")
  Angular_Coordinate, Recal_Xcoord, Recal_Ycoord, vel_min, vel_max = AdimVoxGrid.adimention_MRI_voxel_grid(RowNb, ColumnNb, Time, Xc, Yc,VectXcoord, VectYcoord, R_voxel_max, Radial_Coordinate, VelocityMat, VelocityMat_Parabolic)

  meshfile = "inlet_TOF3.msh"
  x,y,z, lines, triags, tetras = Gmsh.read_gmsh_file(meshfile)
  FinalEdgeList, BoundaryNodes = Gmsh.find_boundary_nodes_SurfMesh(triags)
  
  print("Adimensionalisation of inlet grid ...")
  R_inlet_max, NNodes, X_local, Y_local, r_inlet, Phi_inlet = AdimInletMesh.adimention_Inlet_grid(meshfile, RowNb, ColumnNb, Time, VectXcoord, VectYcoord, VelocityMat, NormalVector)

  print('R-max is: {}'.format(R_inlet_max))
  print("Searching for closest values ....")
  #Velocity_inlet = SearchAlgo.closest_value(NNodes, Time, R_inlet_max, R_voxel_max, r_inlet, \
                                            #Phi_inlet, RowNb, ColumnNb, Radial_Coordinate, BoundaryNodes, Angular_Coordinate, VelocityMat)
  Velocity_inlet = SearchAlgo_3Nodes.closest_value_3Nodes(NNodes, Time, R_inlet_max, R_voxel_max, r_inlet, \
                                            Phi_inlet, RowNb, ColumnNb, Radial_Coordinate, BoundaryNodes, Angular_Coordinate, VelocityMat)
  #Velocity_inlet = SearchAlgo_6Nodes.closest_value_6Nodes(NNodes, Time, R_inlet_max, R_voxel_max, r_inlet, \
                                            #Phi_inlet, RowNb, ColumnNb, Radial_Coordinate, BoundaryNodes, Angular_Coordinate, VelocityMat)
  
# ********************************************************************************************************************************  


##################################################"
#### Plot the solution                         ###"
##################################################"

#----------------------------------------------------------------------------------------------------
NameBase = "3Profiles_"
for k in range(len(Time)):
     output_filename = NameBase + str(k).zfill(4) + "_" + str(int(round(Time[k]))) + ".png"
     print('Saving file {} ...'.format(output_filename))
     Plot.plot_3profiles(VectXcoord, VectYcoord, X_local, Y_local, VelocityMat[:,:,k], VelocityMat_Parabolic[:,:,k], Velocity_inlet[:,k], vel_min, vel_max, output_filename)


# mencoder mf://*.png -mf w=800:h=600:fps=4:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi
print("Encoding images into movie ...")
cmnd = "mencoder mf://" + NameBase + "*.png -mf w=800:h=600:fps=4:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output_3Profiles.avi"
call(cmnd,shell=True)
#----------------------------------------------------------------------------------------------------

NameBase = "1Profiles_onMesh_"
for k in range(len(Time)):
     output_filename = NameBase + str(k).zfill(4) + "_" + str(int(round(Time[k]))) + ".png"
     print('Saving file {} ...'.format(output_filename))
     Plot.plot_1profile_onMesh(X_local, Y_local, Velocity_inlet[:,k], vel_min, vel_max, output_filename)
     
print("Encoding images into movie ...")
#cmnd = "mencoder mf://*.png -mf w=800:h=600:fps=4:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi"
#cmnd = "mencoder mf://output_filename -mf w=800:h=600:fps=4:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi"
cmnd = "mencoder mf://" + NameBase + "*.png -mf w=800:h=600:fps=4:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output_1Profile.avi"

call(cmnd,shell=True)

# ********************************************************************************************************************************  

##################################################################################"
  ### Write Velocity_inlet and points_inlet at a csv file  ###"
  
##################################################################################"



#NormalVector = [0.6843449, -0.5662259, 0.4594129]  ##This is the normal vector to the inlet surface obtained from VMTK
#BaseFolder_path = "/home/maria/Desktop/SRV/3DVelocityProfile/new_code/TOF3_m/"

NameOfInletSurface = "inflow"
OpenFOAMFolder_path = os.path.join(BaseFolder_path, "TOF3_tv", "constant", "boundaryData", NameOfInletSurface)
if not os.path.exists(OpenFOAMFolder_path):
  os.makedirs(OpenFOAMFolder_path)
os.chdir(OpenFOAMFolder_path)
    
    
n = 3 ## Number of decimal point after rounding -- This is used for the name of the folder
k = 1e-3 ## This a multiplying factor to convert the time values from milisecond to second

for t in range(len(Time)):
  Time[t] = k * Time[t]
  TimeFolderName = '{0:.{1}f}'.format(Time[t], n)
  SubFolder_path = os.path.join(OpenFOAMFolder_path,TimeFolderName)
  if not os.path.exists(SubFolder_path):
    os.makedirs(SubFolder_path)
  os.chdir(SubFolder_path)
  WriteFiles.VelocityForOpenFOAM(NNodes, t, Velocity_inlet, NormalVector)
  os.chdir(OpenFOAMFolder_path)
  
WriteFiles.NodesCoordForOpenFOAM(NNodes, x, y, z)


