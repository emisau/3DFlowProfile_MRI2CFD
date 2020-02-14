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


def VelocityForOpenFOAM(NNodes, t, Velocity_inlet, NormalVector):
  
  n = 6 ## n represents the number of decimal points desired after rounding -- this is used to print the value of velocity components in the text files
  k = 1e-2 ## This a multiplying factor to convert the velocity values from cm/s to m/s

  OutputFileName = "U.c"
  outfile = open(OutputFileName,"w")
  outfile.write("// Data on points"+"\n")
  outfile.write(str(NNodes)+"\n")
  outfile.write("\n") ## Skip a line
  outfile.write("("+"\n")
  
  for i_node in range(NNodes):
    NormalVelocity = np.zeros(3) ##This is a dummy vector
    for j in range(3):
      NormalVelocity[j] = k * NormalVector[j] * Velocity_inlet[i_node,t]
    outfile.write("("+'{0:.{1}f}'.format(NormalVelocity[0], n)+ " " +'{0:.{1}f}'.format(NormalVelocity[1], n)+ " " +'{0:.{1}f}'.format(NormalVelocity[2], n)+")"+"\n")
    
  outfile.write("\n")  ## Skip a line 
  outfile.write(")"+"\n")
  outfile.write("\n") ## Skip a line
  outfile.write("// ************************************************************************* //")
  outfile.close()



def NodesCoordForOpenFOAM(NNodes, x, y, z):
  
  n = 10 ## n represents the number of decimal points desired after rounding -- this is used to print the value of velocity components in the text files
  
  OutputFileName = "points.c"
  outfile = open(OutputFileName,"w")
  outfile.write("// Data on points"+"\n")
  outfile.write(str(NNodes)+"\n")
  outfile.write("\n") ## Skip a line
  outfile.write("("+"\n")
  
  for i_node in range(NNodes):
    outfile.write("("+'{0:.{1}f}'.format(x[i_node], n)+ " " +'{0:.{1}f}'.format(y[i_node], n)+ " " +'{0:.{1}f}'.format(z[i_node], n)+")"+"\n")
    
  outfile.write("\n")  ## Skip a line 
  outfile.write(")"+"\n")
  outfile.write("\n") ## Skip a line
  outfile.write("// ************************************************************************* //")
  outfile.close()














