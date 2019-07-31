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

# My own module called 'Gmsh' (functions in file Gmsh.py)
import Gmsh

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.tri as tri

from subprocess import call

def plot_to_file_2profiles(X,Y,Vel1,Vel2,z_min, z_max, filename):
  #ax1 = fig.gca(projection='3d')
  #ax2 = fig.gca(projection='3d')
  #fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, projection='3d')

  #fig = plt.figure()
  # figaspect(0.5) makes the figure twice as wide as it is tall. 
  # Then the *1.5 increases the size of the figure. 
  # The labels etc won't increase so this is a way to make the graph look less cluttered by the labels.
  #fig = plt.figure(figsize=plt.figaspect(0.5)*1.5)
  #ax1 = fig.add_subplot(121, projection='3d')
  #ax2 = fig.add_subplot(122, projection='3d')

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
  #ax1 = fig.gca(projection='3d')
  #ax2 = fig.gca(projection='3d')
  #fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, projection='3d')

  #fig = plt.figure()
  # figaspect(0.5) makes the figure twice as wide as it is tall. 
  # Then the *1.5 increases the size of the figure. 
  # The labels etc won't increase so this is a way to make the graph look less cluttered by the labels.
  #fig = plt.figure(figsize=plt.figaspect(0.5)*1.5)
  #ax1 = fig.add_subplot(121, projection='3d')
  #ax2 = fig.add_subplot(122, projection='3d')

  fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=True, sharey=True, subplot_kw={'projection':'3d'}, 
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


###################################################
### Detect the size of the field and count the number of paheses (time steps)  ###  
### 1st pass reading of the input file                                         ###
###################################################
  data_file = open("3DVelocityShapeProfiles.csv")
  count = 0

  unit_scale_factor = 0.1

  for line in data_file:
    if "TopLeftCorner" in line:
      #print(line)
      temp = line.split(",")
      Xmin = float(temp[1]) * unit_scale_factor
      Ymax = float(temp[2]) * unit_scale_factor
      print('Xmin = {}, Ymax = {}'.format(Xmin,Ymax))
    if "BottomRightCorner" in line:
      #print(line)
      temp = line.split(",")
      Xmax = float(temp[1]) * unit_scale_factor
      Ymin = float(temp[2]) * unit_scale_factor
      print('Xmax = {}, Ymin = {}'.format(Xmax,Ymin))
    if "ResolutionRows" in line:
      #print(line)
      temp = line.split(",")
      ResolY = float(temp[1]) * unit_scale_factor
      print('ResolY = {}'.format(ResolY))
    if "ResolutionCols" in line:
      #print(line)
      temp = line.split(",")
      ResolX = float(temp[1]) * unit_scale_factor
      print('ResolX = {}'.format(ResolX))
    if "Time [ms]" in line:
      count = count+1

  SpanX = (Xmax - Xmin)
  SpanY = (Ymax - Ymin)
  RowNb = int(round(math.fabs(SpanY / ResolY))) + 1
  ColumnNb = int(round(math.fabs(SpanX / ResolX))) + 1
  print('Nb of Rows = {0:2d}'.format(RowNb))
  print('Nb of Columns = {0:2d}'.format(ColumnNb))

  VectXcoord = np.zeros(ColumnNb)
  VectYcoord = np.zeros(RowNb)

  SignX = 1.0
  if Xmin > Xmax:
    SignX = -1.0;
  for i in range (ColumnNb):
    VectXcoord[i] = (Xmin + SignX*i*ResolX)
    
  SignY = 1.0
  if Ymin > Ymax:
    SignY = -1.0;
  for j in range (RowNb):
    VectYcoord[j] = (Ymin + SignY * j * ResolY)

  print('Vector of X-coordinates = {}'.format(VectXcoord))
  print('Vector of Y-coordinates = {}'.format(VectYcoord))
  print('Nb of phases = {}'.format(count))

  data_file.close()

#################################################################
### Store mean velocity matrices for each time for all voxels ###
### 2nd pass reading of the input file                        ###
#################################################################

  Time = np.zeros(count)
  VelocityMat = np.zeros((RowNb,ColumnNb,len(Time)))
  count = 0

  data_file = open("3DVelocityShapeProfiles.csv")
  for line in data_file:
      if "Time [ms]" in line:
        temp = line.split(",")
        Time[count] = float(temp[1])
        
        for i in range(RowNb):
          line = data_file.readline()
          temp = line.split(",")
          for j in range(ColumnNb):
            VelocityMat[i,j,count] = float(temp[j])
          #print(line)

        count = count +1
        
  print('Vector Time = {}'.format(Time))

  data_file.close()


  ## *********************************************************

  #for k in range(len(Time)):
    #output_filename = "output_" + str(k).zfill(4) + "_" + str(int(round(Time[k]))) + ".png"
    #print('Saving file {} ...'.format(output_filename))
    #plot_to_file(VectXcoord, VectYcoord, VelocityMat[:,:,k], output_filename)


  # *********************************************************

########################
### Computation of Elipse area  ###
########################


XL = np.zeros(len(Time))
XR = np.zeros(len(Time))
YT = np.zeros(len(Time))
YB = np.zeros(len(Time))

            
for k in range(len(Time)):
    for j in range (ColumnNb): #Select a "vertical" starting from the left to the right
      norm = np.linalg.norm(VelocityMat[:,j,k], np.inf) #Compute the norm of the j-th column
      if norm != 0: 
        XL[k] = VectXcoord[j] 
        break
    for j in range (ColumnNb-j-1): #Select a "vertical" starting from the right to the left
      norm = np.linalg.norm(VelocityMat[:,ColumnNb-j-1,k], np.inf)
      if norm != 0:
        XR[k] = VectXcoord[ColumnNb -j-1]
        break
    for i in range (RowNb): #Select a "horizontal" starting from the top to the bottom
      norm = np.linalg.norm(VelocityMat[i,:,k], np.inf)
      if norm != 0:
        YT[k] = VectYcoord[i]
        break
    for i in range (RowNb-i-1): #Select a "horizontal" starting from the bottom to the top
      norm = np.linalg.norm(VelocityMat[RowNb-i-1,:,k], np.inf)
      if norm != 0:
        YB[k] = VectYcoord[ColumnNb -j-1]
        break

axis1 = np.zeros(len(Time))
axis2 = np.zeros(len(Time))
elipse_area = np.zeros(len(Time))
R_voxel_max = np.zeros(len(Time))

for k in range (len(Time)):
  axis1[k] = 0.5 * abs(XL[k] - XR[k]) # This is first ellipse half-axis
  axis2[k] = 0.5 * abs(YB[k] - YT[k]) # This is second ellipse half-axis
  elipse_area = math.pi * axis1 * axis2
  R_voxel_max[k] = max(axis1[k], axis2[k])

#####################################
### Computation of velocity integral for each section  ###
#####################################


vel_integral  = np.zeros(len(Time))
for k in range(len(Time)):
  temp = 0
  for i in range(RowNb):
      for j in range(ColumnNb):
          temp = temp + VelocityMat[i,j,k]
  #vel_integral[k] = temp #I verified this expression with the Excel sheet (sum of all entries at t = 157 ms]
  vel_integral[k] = temp * (ResolX)*(ResolY)


print('Vector of velocity integrals = {}'.format(vel_integral))

########################################
### Computation of "global" mean velocity over the section  ###
########################################

MeanVel_perPhase = np.zeros(len(Time))
for k in  range(len(Time)):
  MeanVel_perPhase[k] = vel_integral[k] / elipse_area[k]

print('Vector of mean velocity per section = {}'.format(MeanVel_perPhase))

#################################################
### Computation of of the corresponding parabolic profile for each phase   ###
#################################################

VelocityMat_Parabolic = np.zeros((RowNb,ColumnNb,len(Time)))
Radial_Coordinate = np.zeros((RowNb,ColumnNb,len(Time)))
Angular_Coordinate = np.zeros((RowNb,ColumnNb,len(Time)))
Recal_Xcoord = np.zeros((RowNb,ColumnNb,len(Time)))
Recal_Ycoord = np.zeros((RowNb,ColumnNb,len(Time)))

R = 0
Xc = np.zeros(len(Time))
Yc = np.zeros(len(Time))
dX = 0
dY = 0

for k in range(len(Time)):
   r = 0
   R = np.sqrt(axis1[k]*axis2[k]) #This is in mm
   Xc[k] = 0.5 * (XL[k] + XR[k])
   Yc[k] = 0.5 * (YB[k] + YT[k])
   for i in range(RowNb):
      dY = (VectYcoord[i] - Yc[k])
      for j in range(ColumnNb):
         dX = (VectXcoord[j] - Xc[k])
         r = np.hypot(dX,dY)  ### Equivalent to:  r = np.sqrt(X*X + Y*Y)
         Radial_Coordinate[i,j,k] = r

         if r < R:
           VelocityMat_Parabolic[i,j,k] = 2*MeanVel_perPhase[k]*(1 - (r*r)/(R*R))
         else:
           VelocityMat_Parabolic[i,j,k] = 0
   # Add the other if-conditions

# For control: compute the largest difference between original and computed coordinates
MaxDiffX = 0.0
MaxDiffY = 0.0

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
           
         #Recal_Xcoord[i,j,k] = Xc[k] + Radial_Coordinate[i,j,k] * np.cos(Angular_Coordinate[i,j,k])
         #Recal_Ycoord[i,j,k] = Yc[k] + Radial_Coordinate[i,j,k] * np.sin(Angular_Coordinate[i,j,k]) 
         
         
         Recal_Xcoord[i,j,k] = Radial_Coordinate[i,j,k]/ R_voxel_max[k] * np.cos(Angular_Coordinate[i,j,k])
         Recal_Ycoord[i,j,k] = Radial_Coordinate[i,j,k]/ R_voxel_max[k] * np.sin(Angular_Coordinate[i,j,k]) 
   #print ('AFTER ADIM: the Radial_Coord adiemtionalized at the inlet looks like: {}'.format(Radial_Coordinate[:,:,k]) )
   #print ('Recal_Xcoord and Recal_Ycoord are: {}'.format(Recal_Xcoord[:,:,k], Recal_Ycoord[:,:,k]) )
         
         #print('Original = [{0:.13f},{1:.13f}]'.format(VectXcoord[j], VectYcoord[i]))
         #print('Computed = [{0:.13f},{1:.13f}]\n'.format(Recal_Xcoord[i,j,k], Recal_Ycoord[i,j,k]))

         #MaxDiffX = max(MaxDiffX,abs(Recal_Xcoord[i,j,k]-VectXcoord[j]))
         #MaxDiffY = max(MaxDiffY,abs(Recal_Ycoord[i,j,k]-VectYcoord[i]))

#print('--------------------------------------------------------------------')
#print('Maximum difference between original and recomputed coordinates:')
#print('dx= {0:.18f}, dy = {1:.18f}'.format(MaxDiffX,MaxDiffY))
#print('--------------------------------------------------------------------\n')



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

targetmesh = "PlanarDisk_2D.msh"
x,y,z, lines, triags, tetras = Gmsh.read_gmsh_file(targetmesh)

### Compute the center of the inlet surface 
NNodes = len(x)

#for i in range (NNodes):
  #Sum_x = Sum_x + x[i]
  #Sum_y = Sum_y + y[i]
  #Sum_z = Sum_z + z[i]

#Xc_inlet = Sum_x / NNodes
#Yc_inlet = Sum_y / NNodes
#Zc_inlet = Sum_z / NNodes
#print ('The coordinates of the inlet center are: Xc = {}, Yc = {}, Zc = {}'.format(Xc_inlet, Yc_inlet, Zc_inlet))



#NEW VERSION OF COMPUTING THE GRAVITY CENTER
Xc_inlet, Yc_inlet, Zc_inlet = Gmsh.surface_gravity_center(triags, x, y, z)
print('The new gravity center = {}, {}, {}'.format(Xc_inlet, Yc_inlet, Zc_inlet))



#################################################
### Define the local coordinate system
### on the surface of the inlet
#################################################

##Vector of the base:
X_vect = [1.4,0,0]
Y_vect = [0,1.4,0]
N_vect = [0,0,1]

## Change of coordinates:
r_inlet = np.zeros(NNodes)
Phi_inlet = np.zeros(NNodes)
X_local = np.zeros(NNodes)
Y_local = np.zeros(NNodes)

# Transformation from Cartesian to cylindrical coordinates (x,y) --> (r, phi)
# See coordinates transformation https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations
R_max = 0.0
for i_node in range (NNodes):
  dX = x[i_node] - Xc_inlet
  dY = y[i_node] - Yc_inlet
  dZ = z[i_node] - Zc_inlet
  r_inlet[i_node] = np.sqrt(dX*dX + dY*dY + dZ*dZ) #This is the magnitude of "r" - the norm of "r"
  if (r_inlet[i_node] > R_max):
    R_max = r_inlet[i_node]
  r_vect_comp = [dX,dY,dZ] #May be the computation of r_vect_comp is useless
  dot_product_rX = (dX*X_vect[0] + dY*X_vect[1] + dZ*X_vect[2])
  dot_product_rY = (dX*Y_vect[0] + dY*Y_vect[1] + dZ*Y_vect[2])
  
  arccos_argument = dot_product_rX / (r_inlet[i_node] * np.linalg.norm(X_vect))
  
  if ( dot_product_rY  >= 0.0):
    Phi_inlet[i_node] = np.arccos(arccos_argument)
  else:
    Phi_inlet[i_node] = 2.0 * math.pi - np.arccos(arccos_argument)

#print ('BEFORE ADIM: the r-vector at the inlet looks like: {}'.format(r_inlet))

print ('the averaged-radius of the inlet is: {}'.format(R_max))

for i_node in range (NNodes):
  X_local[i_node] = r_inlet[i_node] * np.cos(Phi_inlet[i_node])   #Parametric coord system --> u
  Y_local[i_node] = r_inlet[i_node] * np.sin(Phi_inlet[i_node])   #Parametric coord system --> v


#################################################
### Computation of closest value of an inlet point
### to a voxel-grid point
#################################################


temp_dR = 0
i_idx = 0
j_idx = 0
Velocity_inlet = np.zeros((NNodes,len(Time)))


for k in range(len(Time)):

  for i_node in range (NNodes): # pick one node of the inlet
     X_loc = r_inlet[i_node] / R_max * np.cos(Phi_inlet[i_node])
     Y_loc = r_inlet[i_node] / R_max * np.sin(Phi_inlet[i_node])

    
     MinDist= 10000.0 * R_max
    
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

'''
  print("VelocityMat:")
  print(VelocityMat[:,:,k])
  print("Velocity_inlet:")
  print(Velocity_inlet[:,k])
  print("**********************************************************************************************")
'''


##################################################
#### Plot the solution   ###
##################################################
for k in range(len(Time)):
     output_filename = "output_" + str(k).zfill(4) + "_" + str(int(round(Time[k]))) + ".png"
     print('Saving file {} ...'.format(output_filename))
     #plot_to_file(VectXcoord, VectYcoord,VelocityMat[:,:,k],VelocityMat_Parabolic[:,:,k], vel_min, vel_max, output_filename)
     #plot_to_file_3profiles(VectXcoord, VectYcoord, VectXcoord, VectYcoord, VelocityMat[:,:,k], VelocityMat_Parabolic[:,:,k], VelocityMat[:,:,k], vel_min, vel_max, output_filename)
     plot_to_file_3profiles(VectXcoord, VectYcoord, X_local, Y_local, VelocityMat[:,:,k], VelocityMat_Parabolic[:,:,k], Velocity_inlet[:,k], vel_min, vel_max, output_filename)


# mencoder mf://*.png -mf w=800:h=600:fps=4 :type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi
print("Encoding images into movie ...")
cmnd = "mencoder mf://*.png -mf w=800:h=600:fps=4 :type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi"
call(cmnd,shell=True)



#1) given Xvoxel, Yvoxel 
#2) Loop over the matrix with all Phi and with all r and select such pair (Phi, r) that:
 #the distance of [Xcenter + r * cos(Phi), Ycenter + r*shi(Phi)] - [Xvoxel, Yvoxel] is minimal
 
#you would not need Phi if your velocity profile were perfectly symmetric



