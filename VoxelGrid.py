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


def read_MRI_voxel_grid(filename):

  ##################################################################################"
  ### Detect the size of the field and count the number of phases (time steps)  ###"
  ### 1st pass reading of the input file                                         ###"
  ##################################################################################"
  data_file = open(filename, 'r')

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

  #################################################################"
  ### Store mean velocity matrices for each time for all voxels ###"
  ### 2nd pass reading of the input file                        ###"
  #################################################################"

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

  return ResolX, ResolY, RowNb, ColumnNb, VectXcoord, VectYcoord, Time, VelocityMat
