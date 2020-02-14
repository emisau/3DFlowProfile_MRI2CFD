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


from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.tri as tri



def plot_3profiles(X1,Y1,X2,Y2,Vel1,Vel2,Vel3,z_min, z_max, filename):

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

  
  #-----------------------------------------------------------------------------
  # Plot the triangulation and the high-res iso-contours
  #-----------------------------------------------------------------------------

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
  
# ********************************************************************************************************************************  
  
def plot_1profile_onMesh(X1, Y1, Vel1, z_min, z_max, filename):

  #fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=True, sharey=True, subplot_kw={'projection':'3d'}, 
                                 #figsize=plt.figaspect(0.5)*1.5)
  #fig.tight_layout()
  fig = plt.figure()
  ax = fig.gca(projection='3d')
  triang = tri.Triangulation(X1, Y1)

 #-----------------------------------------------------------------------------
  # Refine data
  #-----------------------------------------------------------------------------
  refiner = tri.UniformTriRefiner(triang)
  tri_refi, Vel1_refi = refiner.refine_field(Vel1, subdiv=3)
  

  #-----------------------------------------------------------------------------
  # Plot the triangulation and the high-res iso-contours
  #-----------------------------------------------------------------------------

  ax.triplot(triang, lw=0.5, color='black')
  ax.set_zlim(z_min, z_max)

  levels = np.arange(z_min, z_max, 0.025)
  cmap = cm.get_cmap(name='terrain', lut=None)
  #plt.tricontourf(tri_refi, Vel3_refi.flatten(), levels=levels, cmap=cmap)
  ax.plot_trisurf(X1, Y1, Vel1, cmap=cm.jet, linewidth=0.2)
 
  ax.zaxis.set_major_locator(LinearLocator(10))
  ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
  
  #fig.colorbar(surf, shrink=0.5, aspect=5)


  #plt.show()
  fig.savefig(filename)
  plt.close()
