import numpy as np
import math

# ********************************************************************************************************************************  

def extract_ellipse_area(ResolX, ResolY, RowNb, ColumnNb, Time, VectXcoord, VectYcoord, VelocityMat):
  
  ###################################"
  ### Computation of Elipse area  ###"
  ###################################"


  XL = np.zeros(len(Time))
  XR = np.zeros(len(Time))
  YT = np.zeros(len(Time))
  YB = np.zeros(len(Time))


  for k in range(len(Time)):
    for j in range (ColumnNb): #Select a "vertical" starting from the left to the right
      norm = np.linalg.norm(VelocityMat[:,j,k], np.inf) #Compute the norm of the j-th column
      if norm != 0: 
        XL[k] = VectXcoord[j]
        #XL[k] = XL[k] - 2*ResolX
        break
    for j in range (ColumnNb-j-1): #Select a "vertical" starting from the right to the left
      norm = np.linalg.norm(VelocityMat[:,ColumnNb-j-1,k], np.inf)
      if norm != 0:
        XR[k] = VectXcoord[ColumnNb -j-1]
        #XR[k] = XR[k] + 2*ResolX
        break
    for i in range (RowNb): #Select a "horizontal" starting from the top to the bottom
      norm = np.linalg.norm(VelocityMat[i,:,k], np.inf)
      if norm != 0:
        YT[k] = VectYcoord[i]
        #YT[k] = YT[k] - 2*ResolY
        break
    for i in range (RowNb-i-1): #Select a "horizontal" starting from the bottom to the top
      norm = np.linalg.norm(VelocityMat[RowNb-i-1,:,k], np.inf)
      if norm != 0:
        YB[k] = VectYcoord[ColumnNb -i-1]
        #YB[k] = YB[k] + 2*ResolY
        break

  axis1 = np.zeros(len(Time))
  axis2 = np.zeros(len(Time))
  elipse_area = np.zeros(len(Time))
  R_voxel_max = np.zeros(len(Time))

  for k in range (len(Time)):
    axis1[k] = 0.5 * abs(XL[k] - XR[k]) # This is first ellipse half-axis
    axis2[k] = 0.5 * abs(YB[k] - YT[k]) # This is second ellipse half-axis
    elipse_area[k] = math.pi * axis1[k] * axis2[k]
    R_voxel_max[k] = max(axis1[k], axis2[k])
  
  return XR,XL,YT,YB, axis1, axis2, elipse_area, R_voxel_max





# ********************************************************************************************************************************  


def compute_equivalent_parabolic(ResolX, ResolY, RowNb, ColumnNb, Time, VectXcoord, VectYcoord, XR, XL, YT, YB, axis1, axis2, elipse_area, VelocityMat):

  ##########################################################"
  ### Computation of velocity integral for each section  ###"
  ##########################################################"


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

  return Xc, Yc, MeanVel_perPhase, VelocityMat_Parabolic, Radial_Coordinate
