import numpy as np

#########################################
### Read a GMSH file
### Store connectivity and node coord
#########################################

def read_gmsh_file(filename):

  # Array which remembers the number of nodes for each element type defined in gmsh:
  gmsh_elem_nodes = np.array([2,3,4,4,8,6,5,3,6,9,10,27,18,14,1,8,20,15,13,9,10,12,15,15,21,4,5,6,20,35,56])

  print ('Starting to read the input file ...')
  mesh_file = open(filename, 'r')
  ##mesh_file = open("PlanarDisk_2D.msh")

  buffer = "";

  ## START READING THE FILE

  while buffer.find('Nodes') != 1 :
    buffer = mesh_file.readline();

  NNodes = int(mesh_file.readline())
  print('Nb of nodes = {0:2d}'.format(NNodes))


  x = np.zeros(NNodes,float)
  y = np.zeros(NNodes,float)
  z = np.zeros(NNodes,float)

  for inode in range(NNodes):
    nodeinfo = (mesh_file.readline()).split() # Line containing all coordinates
    x[inode] = float(nodeinfo[1])
    y[inode] = float(nodeinfo[2])
    z[inode] = float(nodeinfo[3])

  print ('\t...done')

  #for inode in range(NNodes):
  #  print "x[%d] = %f" % (inode,x[inode])

  buffer = "";

  medit_idx = np.zeros(NNodes,int)

  while buffer.find('Elements') != 1 :
    buffer = mesh_file.readline()

  NElements = int(mesh_file.readline())

  print('Nb of elements = {0:2d}'.format(NElements))


  ############## FIRST READING - WE USE IT TO COUNT THE ELEMENT TYPES ############## 
  # (HOW MANY LINES, TRIANGLES ETC.)

  NrP1Line = 0
  NrP1Tri = 0
  NrP1Tet = 0

  for ielem in range(NElements):
    eleminfo = (mesh_file.readline()).split()
    elem_nr = int(eleminfo[0])
    elem_type = int(eleminfo[1])

    if elem_type == 1: # Simple line segment (P1 line)
      NrP1Line = NrP1Line + 1
  
    elif elem_type == 2: # 3-node triangle
      NrP1Tri = NrP1Tri + 1

    elif elem_type == 4:
      NrP1Tet = NrP1Tet+1

  mesh_file.close() # The first pass finished, start the second:

  print('Nb of elements = {0:2d}'.format(NElements))

  mesh_file = open(filename, 'r')

  buffer = "";

  ## START READING THE FILE

  while buffer.find('Elements') != 1 :
    buffer = mesh_file.readline();

  buffer = mesh_file.readline(); #This reads the number of elements again
                              #we have this information from the first reading

  # Create matrices which hold connectivity:
  LineConnectivity = np.zeros( (NrP1Line,3), int )
  TriagConnectivity = np.zeros( (NrP1Tri,4), int )
  TetraConnectivity = np.zeros( (NrP1Tet,5), int )

  iline = 0
  itriag = 0
  itetra = 0

  for ielem in range(NElements):
    eleminfo = (mesh_file.readline()).split()
    elem_nr = int(eleminfo[0])
    elem_type = int(eleminfo[1])
    nr_tags = int(eleminfo[2])
    phys_tag = int(eleminfo[3])


    if elem_type == 1: # 2-node line
      LineConnectivity[iline,0] = int(eleminfo[2+nr_tags+1]) - 1 # 2 nodes of the line
      LineConnectivity[iline,1] = int(eleminfo[2+nr_tags+2]) - 1
      LineConnectivity[iline,2] = phys_tag # physical tag of triag needed for medit file format
      iline = iline + 1

    elif elem_type == 2: # 3-node triangle
      # We save the information about the nodes of the triangle - they will have
      # the same physical entity index as the triangle
      TriagConnectivity[itriag,0] = int(eleminfo[2+nr_tags+1]) - 1 # 3 nodes of the triag
      TriagConnectivity[itriag,1] = int(eleminfo[2+nr_tags+2]) - 1
      TriagConnectivity[itriag,2] = int(eleminfo[2+nr_tags+3]) - 1
      TriagConnectivity[itriag,3] = phys_tag # physical tag of triag needed for medit file format
      itriag = itriag + 1

    elif elem_type == 4: # 4-node tetrahedron
      TetraConnectivity[itetra,0] = int(eleminfo[2+nr_tags+1]) - 1 # 4 nodes of the tetra
      TetraConnectivity[itetra,1] = int(eleminfo[2+nr_tags+2]) - 1
      TetraConnectivity[itetra,2] = int(eleminfo[2+nr_tags+3]) - 1
      TetraConnectivity[itetra,3] = int(eleminfo[2+nr_tags+4]) - 1
      TetraConnectivity[itetra,4] = phys_tag # physical tag of tetra needed for medit file format
      itetra = itetra + 1

  mesh_file.close() # Second pass finished
  
  return x,y,z, LineConnectivity, TriagConnectivity, TetraConnectivity

# =============================================================================

def surface_gravity_center(triags,xcoord,ycoord,zcoord):

  tot_surface = 0.0

  grav_center = np.array([0.0, 0.0, 0.0])

  for i in range(triags.shape[0]):
     # Get coordinates of vertices of each triangle
     x0 = xcoord[triags[i,0]]
     x1 = xcoord[triags[i,1]]
     x2 = xcoord[triags[i,2]]

     y0 = ycoord[triags[i,0]]
     y1 = ycoord[triags[i,1]]
     y2 = ycoord[triags[i,2]]
     
     z0 = zcoord[triags[i,0]]
     z1 = zcoord[triags[i,1]]
     z2 = zcoord[triags[i,2]]

     # To compute the surface, evaluate
     # the vector product of two edges and divide by 2

     edge0 = [x1-x0, y1-y0, z1-z0]
     edge1 = [x2-x1, y2-y1, z2-z1]

     cross = np.cross(edge0, edge1)
     tri_area = 0.5 * np.linalg.norm(cross)
     tri_barycenter = 1.0/3.0 * np.array([x0 + x1 + x2, y0 + y1 + y2, z0 + z1 + z2])

     tot_surface = tot_surface + tri_area
     grav_center = grav_center + tri_area * tri_barycenter

  # Finally, divide by total surface area
  grav_center = 1./tot_surface * grav_center

  return grav_center[0], grav_center[1], grav_center[2]

# =============================================================================

def find_boundary_nodes_SurfMesh(triags):
  
  NTriag = triags.shape[0]
  NEdges = 3 * NTriag
  
  # TmpEdgeList is a matrix which has as many rows as number of edges in the mesh
  # and has 4 columns: 
  # Begin Node: beginning node of the edge, 
  # End Node: end node of the edge, 
  # Left Elem: element ID on the left of the edge, 
  # Right Elem: elment ID on the right of the edge
  TmpEdgeList = np.zeros((NEdges, 4), int) 
  for i_edge in range (NEdges):
    TmpEdgeList[i_edge,:] = -1.0
    #TmpEdgeList[i_edge,3] = -1.0
              
  lines_filled = 0
  for i_triag in range (NTriag):
    for i_node in range(3):
      node_0 = triags[i_triag,i_node]
      node_1 = triags[i_triag,(i_node + 1) % 3]
      sibling_found = False
      for i_edge in range (NEdges):
        if (TmpEdgeList[i_edge,0] == node_1) and (TmpEdgeList[i_edge,1] == node_0) :
          TmpEdgeList[i_edge,3] = i_triag
          sibling_found = True
          break
      if (sibling_found == False):
        TmpEdgeList[lines_filled,0] = node_0
        TmpEdgeList[lines_filled,1] = node_1
        TmpEdgeList[lines_filled,2] = i_triag                  
        lines_filled = lines_filled + 1
         
  FinalEdgeList = np.zeros((lines_filled, 4), int)
  # Store only that part of TmpEdgeList which was actually filled
  # with edge data
  FinalEdgeList = TmpEdgeList[:lines_filled,:]

  # The number of nodes must be equal to the maximum node number
  # in triangle connectivity
  NNodes = FinalEdgeList[:,:2].max() + 1

  BoundaryNodes = np.zeros(NNodes, bool)
  BoundaryNodes[:] = False
  for i_edge in range (lines_filled):
    if(FinalEdgeList[i_edge,3] == -1):
      BoundaryNodes[FinalEdgeList[i_edge,0]] = True
      BoundaryNodes[FinalEdgeList[i_edge,1]] = True

  #print(BoundaryNodes)
  #print(FinalEdgeList)
  #print("Number of edges = {}".format(lines_filled))
  #print("Number of nodes = {}".format(NNodes))


         
  return FinalEdgeList, BoundaryNodes
    
