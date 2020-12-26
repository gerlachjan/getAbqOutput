# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 12:44:19 2020

@author: Jan Gerlach
"""

from abaqus import *
from odbAccess import*
from abaqusConstants import*
from odbSection import *

import __main__
import visualization

import numpy as np
import math


def unVoigt3d(a):
    
    if a.shape == (6,):
        
        sig = np.zeros(shape=(3,3), dtype=float)
        vorfak = [1, 1, 1, 1, 1, 1]
        indizes = [(0,0), (1,1), (2,2), (0,1), (0,2), (1,2)]
        
        for k in range(6):
            n,m = indizes[k]
            sig[n,m] = sig[m,n] = a[k]/vorfak[k]
            
        return sig
    
    else:
        print('Error in tensor transformation: voigt(,6) to tensor(3,3)')



def Voigt3d(a):
    
    if a.shape == (3,3):
        
        sig = np.zeros(shape=(6,))
        vorfak = [1, 1, 1, 1, 1, 1]
        indizes = [(0,0), (1,1), (2,2), (0,1), (0,2), (1,2)]
        
        for k in range(6):
            n,m = indizes[k]
            sig[k] = vorfak[k]*a[n,m]
        
        return sig
    
    
    else:
        
        print('Error in tensor transformation: tensor(3,3) to voigt(,6)')
        
        
def tensor_swap(init_tensor,swaptype):
    swaped_tensor = np.zeros((3,3))
    
    if swaptype == 'xy':
        swaped_tensor[0,0] = init_tensor[1,1] #xx to yy
        swaped_tensor[1,1] = init_tensor[0,0] #yy to xx
        swaped_tensor[2,2] = init_tensor[2,2] #zz to zz
        swaped_tensor[0,1] = init_tensor[1,0] #xy to yx
        swaped_tensor[0,2] = init_tensor[1,2] #xz to yz
        swaped_tensor[1,2] = init_tensor[0,2] #yz to xz
    elif swaptype =='xz':
        swaped_tensor[0,0] = init_tensor[2,2] #xx to zz
        swaped_tensor[1,1] = init_tensor[1,1] #yy to yy
        swaped_tensor[2,2] = init_tensor[0,0] #zz to xx
        swaped_tensor[0,1] = init_tensor[2,1] #xy to zy
        swaped_tensor[0,2] = init_tensor[2,0] #xz to zx
        swaped_tensor[1,2] = init_tensor[1,0] #yz to yx
    elif swaptype =='yz':
        swaped_tensor[0,0] = init_tensor[0,0] #xx to xx
        swaped_tensor[1,1] = init_tensor[2,2] #yy to zz
        swaped_tensor[2,2] = init_tensor[1,1] #zz to yy
        swaped_tensor[0,1] = init_tensor[0,2] #xy to xz
        swaped_tensor[0,2] = init_tensor[0,1] #xz to xy
        swaped_tensor[1,2] = init_tensor[2,1] #yz to zy
        
    else:
        print('Wrong swaptype')
    
    #symmetry of tensor
    swaped_tensor[1,0] = swaped_tensor[0,1] #yx = xy
    swaped_tensor[2,0] = swaped_tensor[0,2] #zx = xz
    swaped_tensor[2,1] = swaped_tensor[1,2] #zy = yz
        
    return swaped_tensor

def mirrorData(init_tensor_list):
    tensor_xy_list = []
    tensor_xz_list = []
    tensor_yz_list = []
    
    for tensor in init_tensor_list:
        swaped_tensor_xy = tensor_swap(tensor,'xy')
        swaped_tensor_xz = tensor_swap(tensor,'xz')
        swaped_tensor_yz = tensor_swap(tensor,'yz')
        
        tensor_xy_list.append(swaped_tensor_xy)
        tensor_xz_list.append(swaped_tensor_xz)
        tensor_yz_list.append(swaped_tensor_yz)
        
    voigt_xy_list   = []
    voigt_xz_list   = []
    voigt_yz_list   = []
    
    for tensor in tensor_xy_list:    
        voigt_xy_list.append(Voigt3d(tensor))
    for tensor in tensor_xz_list:    
        voigt_xz_list.append(Voigt3d(tensor)) 
    for tensor in tensor_yz_list:    
        voigt_yz_list.append(Voigt3d(tensor)) 
        
    return voigt_xy_list,voigt_xz_list,voigt_yz_list
######################################################################################
######################################################################################
#main program

#validation 0 = false, 1 = true
validation = 0

#initialization of elemtype, path, set- and instance names
elem_type='C3D8R'
#elem_type = 'CAX4R'
#elem_type = 'CPS4R'

odbname='tensile_r5_t2p5_imp'


import os

try:
    dirname = os.path.dirname(os.path.abspath(__file__))
except NameError:  # We are the main py2exe script, not a module
    import sys
    dirname = os.path.dirname(os.path.abspath(sys.argv[0]))

#dirname = os.path.dirname(__file__)
#workdir = os.path.join(dirname,'output')
#os.chdir(workdir)

#odbpath = os.path.join(dirname,'odb_files',odbname + '.' + 'odb')

odb_path = 'odb_files/tensile_r5_t2p5_imp.odb'
instance_name = 'PART-1-1'

#output dimension
dim = 3
#dim = 2

#data extraction at following sets
elem_setname = 'OUTPUT_ELEM'
node_setname = 'OUTPUT_NODE'


odb = openOdb(path=odb_path)
scratchOdb = session.ScratchOdb(odb)
currentSession = session.odbs[odb_path]

#extract data at the elements integration/node points at Part Level
rootElem=odb.rootAssembly.instances[instance_name].elementSets[elem_setname]  #sets are defined at Part-Level
rootNode=odb.rootAssembly.instances[instance_name].nodeSets[node_setname]     #sets are defined at Part-Level

#output at lastframe
lastFrame = odb.steps['Step-1'].frames[-1]



stressfield_IP      = lastFrame.fieldOutputs['S'].getSubset(region=rootElem, position=INTEGRATION_POINT, elementType=elem_type).values
strainfield_IP      = lastFrame.fieldOutputs['LE'].getSubset(region=rootElem, position=INTEGRATION_POINT, elementType=elem_type).values
peeqstrainfield_IP  = lastFrame.fieldOutputs['PEEQ'].getSubset(region=rootElem, position=INTEGRATION_POINT, elementType=elem_type).values
coordfield          = lastFrame.fieldOutputs['COORD'].getSubset(region=rootNode,position=NODAL).values

stress_list         = []
maxPrincipal_list   = []
minPrincipal_list   = []
midPrincipal_list   = []
maxPrincipal_strain_list   = []
minPrincipal_strain_list   = []
midPrincipal_strain_list   = []
abs_maxPrincipal_list = []
mises_list          = []
elemLabel_list      = []
strain_list         = []
peeq_list           = []
coordinate_list     = []
nLabel_list         = []
connectivity_list   = []


    
for elem in rootElem.elements:
    connectivity_list.append(elem.connectivity)
    
connectivity_list   = np.array(connectivity_list)  


#create local coordinate system
local_nodes = []
systeme = []
elementnumber = 0
coordSys_matrix_list = []
for connectivity in connectivity_list:
    node_xy_plane = []
    elementnumber = elementnumber +1
    session.viewports['Viewport: 1'].setValues(displayedObject=currentSession)
    session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
        DEFORMED, ))
    node1 = odb.rootAssembly.instances[instance_name].nodes[connectivity[0]-1]
    z_initial = node1.coordinates[2]
    for element in connectivity:
        node_temp = odb.rootAssembly.instances[instance_name].nodes[element-1]
        z_coord_temp = node_temp.coordinates[2]
        dz = abs(z_coord_temp - z_initial)
        if element == connectivity[0]:
            continue
        elif dz < 0.001:
            node_xy_plane.append(node_temp)
            
        
    node2 = node_xy_plane[0]
    node3 = node_xy_plane[1]            
                
    #scratchOdb.rootAssembly.DatumCsysByThreeNodes(name='CSYS-Elem'+str(elementnumber), coordSysType=CARTESIAN, origin=node1, point1=node2, point2=node3)   
    coordSys = odb.rootAssembly.DatumCsysByThreeNodes(name='CSYS-Elem'+str(elementnumber), coordSysType=CARTESIAN, origin=node1, point1=node2, point2=node3)
    coordSys_matrix = np.array([[coordSys.xAxis],
                                [coordSys.yAxis],
                                [coordSys.zAxis]])
    coordSys_matrix = coordSys_matrix.reshape(3,3)
    coordSys_matrix = np.transpose(coordSys_matrix)
    coordSys_matrix_list.append(coordSys_matrix)
    
coordSys_matrix_list = np.array(coordSys_matrix_list)
    
     


for elem in stressfield_IP:
    stress              = elem.data
    stress_list.append(stress)
    
    maxPrincipal_list.append(elem.maxPrincipal)
    minPrincipal_list.append(elem.minPrincipal)
    midPrincipal_list.append(elem.midPrincipal)
    
    if abs(elem.maxPrincipal) > abs(elem.minPrincipal):
        abs_maxPrincipal_list.append(elem.maxPrincipal)
    else:
        abs_maxPrincipal_list.append(0)
    
    stress_mises        = elem.mises
    mises_list.append(stress_mises)
    
    elemLabel_list.append(elem.elementLabel)

  
for elem in strainfield_IP:
    strain = elem.data
    strain_list.append(strain)
    
    maxPrincipal_strain_list.append(elem.maxPrincipal)
    minPrincipal_strain_list.append(elem.minPrincipal)
    midPrincipal_strain_list.append(elem.midPrincipal)

for elem in peeqstrainfield_IP:
    peeq = elem.data
    peeq_list.append(peeq)

for node in coordfield:
    if node.nodeLabel in nLabel_list:
        continue
    else:
        nLabel_list.append(node.nodeLabel)
        coordinates = node.data
        coordinate_list.append(coordinates)
       
stress_list         = np.array(stress_list)
maxPrincipal_list   = np.array(maxPrincipal_list)
mises_list          = np.array(mises_list)
elemLabel_list      = np.array(elemLabel_list)
strain_list         = np.array(strain_list)    
peeq_list           = np.array(peeq_list)
coordinate_list     = np.array(coordinate_list)   
nLabel_list         = np.array(nLabel_list)         


#add additional columns for 13 and 23 components
if dim == 2:
    zero_row = np.zeros(len(nLabel_list))
    two_zero_row = np.zeros((len(elemLabel_list),2))
    
    coordinate_list = np.column_stack((coordinate_list,zero_row))
    stress_list = np.column_stack((stress_list,two_zero_row))
    strain_list = np.column_stack((strain_list,two_zero_row))


IP_coordinate_list = []
for connectivity in connectivity_list:
    elem_coordinates =[]
    for node in connectivity:
        idx = np.where(nLabel_list == node)
        elem_coordinates.append(coordinate_list[idx[0][0]])
    
    elem_coordinates = np.array(elem_coordinates)
    IP_coordinate_list.append(np.mean(elem_coordinates,axis=0))

IP_coordinate_list = np.array(IP_coordinate_list)  

stress_list_unVoigt = []
stress_hyd_list_unVoigt = []
stress_hyd_scalar_list = []
stress_dev_list_unVoigt = []         
for stress_elem in stress_list:
    sig = unVoigt3d(stress_elem)
    sig_h = 1.0/3.0*np.eye(3)*np.trace(sig)
    sig_h_scalar = 1.0/3.0*np.trace(sig)
    sig_dev = sig - sig_h
    stress_list_unVoigt.append(sig)
    stress_hyd_list_unVoigt.append(sig_h)
    stress_hyd_scalar_list.append(sig_h_scalar)
    stress_dev_list_unVoigt.append(sig_dev)

stress_hyd_list_unVoigt = np.array(stress_hyd_list_unVoigt)
stress_hyd_scalar_list = np.array(stress_hyd_scalar_list)
stress_list_unVoigt = np.array(stress_list_unVoigt)
stress_dev_list_unVoigt = np.array(stress_dev_list_unVoigt)


strain_tensor_list = []
for strain_elem in strain_list:
    eps = unVoigt3d(strain_elem)
    strain_tensor_list.append(eps)

strain_tensor_list = np.array(strain_tensor_list)

#%%###################################################################################
###################################################################################
#calculation of invariants
inv1_list = []
inv1_dev_list = []
inv2_dev_list = []
inv3_dev_list = []


for tensor in stress_dev_list_unVoigt:
    inv1 = np.trace(tensor)
    inv1_list.append(inv1)

for tensor in stress_dev_list_unVoigt:
    inv1_dev = np.trace(tensor)
    inv2_dev = 1.0/6.0*((tensor[0,0]-tensor[1,1])**2 + (tensor[0,0]-tensor[2,2])**2 + (tensor[1,1]-tensor[2,2])**2) + tensor[0,1]**2 + tensor[0,2]**2 + tensor[1,2]**2
    inv3_dev = np.linalg.det(tensor)

    inv1_dev_list.append(inv1_dev)
    inv2_dev_list.append(inv2_dev)
    inv3_dev_list.append(inv3_dev)
inv1_dev_list = np.array(inv1_dev_list)
inv2_dev_list = np.array(inv2_dev_list)
inv3_dev_list = np.array(inv3_dev_list)



stress_tensor_local_list=[]
strain_tensor_local_list=[]
for idx in range(len(coordSys_matrix_list)):
    stress_tensor_local = np.tensordot(np.tensordot(np.transpose(coordSys_matrix_list[idx]),stress_list_unVoigt[idx],axes=1),coordSys_matrix_list[idx],axes=1)
    stress_tensor_local_list.append(stress_tensor_local)
    
    strain_tensor_local = np.tensordot(np.tensordot(np.transpose(coordSys_matrix_list[idx]),strain_tensor_list[idx],axes=1),coordSys_matrix_list[idx],axes=1)
    strain_tensor_local_list.append(strain_tensor_local)

stress_tensor_local_list = np.array(stress_tensor_local_list)
strain_tensor_local_list = np.array(strain_tensor_local_list)


###################################################################
###################################################################

stress_voigt_xy_local_list, stress_voigt_xz_local_list,stress_voigt_yz_local_list = mirrorData(stress_tensor_local_list)
strain_voigt_xy_local_list, strain_voigt_xz_local_list,strain_voigt_yz_local_list = mirrorData(strain_tensor_local_list)
strain_voigt_xy_list, strain_voigt_xz_list,strain_voigt_yz_list = mirrorData(strain_tensor_list)
stress_voigt_xy_list, stress_voigt_xz_list,stress_voigt_yz_list = mirrorData(stress_list_unVoigt)


stress_voigt_local_list =[]
for tensor in stress_tensor_local_list:    
    stress_voigt_local_list.append(Voigt3d(tensor)) 

strain_voigt_local_list =[]
for tensor in strain_tensor_local_list:    
    strain_voigt_local_list.append(Voigt3d(tensor)) 

    
stress_voigt_local_list = np.array(stress_voigt_local_list)
strain_voigt_local_list = np.array(strain_voigt_local_list)

#calculation of stress triaxiality + lode parameter
eta = inv1_list/(3.0*np.sqrt(3.0*inv2_dev_list))
eta = stress_hyd_scalar_list/(mises_list)  #triaxiality
xi = 3.0*np.sqrt(3.0)/2.0*inv3_dev_list/inv2_dev_list**(3.0/2.0) #Lode parameter
#theta = np.arccos(xi)/3.0 #lode angle

####################################################################################
####################################################################################
                   

#%%calculation of damage parameter(normalized cockcroft-latham criterion)
damage = (abs_maxPrincipal_list/mises_list)*peeq_list
#damage = np.zeros(len(elemLabel_list))


#########################################################################
#########################################################################
if validation == True:
    maxPrincipal_list = np.array(maxPrincipal_list)
    midPrincipal_list = np.array(midPrincipal_list)
    minPrincipal_list = np.array(minPrincipal_list)
    mises_strain_principal = math.sqrt(0.5*((maxPrincipal_list-midPrincipal_list)**2 + (maxPrincipal_list-minPrincipal_list)**2 + (midPrincipal_list-minPrincipal_list)**2))
    stress_tensor_local_list
    mises_stress_local = math.sqrt(0.5*((stress_voigt_local_list[0][0] - stress_voigt_local_list[0][1])**2 + (stress_voigt_local_list[0][0] - stress_voigt_local_list[0][2])**2 + (stress_voigt_local_list[0][1] - stress_voigt_local_list[0][2])**2 + 6*(stress_voigt_local_list[0][3]**2 + stress_voigt_local_list[0][4]**2 + stress_voigt_local_list[0][5]**2)))                          
    mises_stress_global = math.sqrt(0.5*((stress_list[0][0] - stress_list[0][1])**2 + (stress_list[0][0] - stress_list[0][2])**2 + (stress_list[0][1] - stress_list[0][2])**2 + 6*(stress_list[0][3]**2 + stress_list[0][4]**2 + stress_list[0][5]**2)))
    print(['v.Mises Principal: ',mises_strain_principal])
    print(['v.Mises Local: ',mises_stress_local])
    print(['v.Mises Global: ',mises_stress_global])
    print(['v.Mises Abaqus: ',mises_stress_global])
    print(['local stress components: ',stress_voigt_local_list])
    
#########################################################################
#########################################################################

         
#%%generate output



merged_data_strain_stress = np.column_stack((elemLabel_list,IP_coordinate_list,strain_list,stress_list))

merged_data_local_strain_stress = np.column_stack((elemLabel_list,IP_coordinate_list,strain_voigt_local_list,stress_voigt_local_list))
merged_data_local_xy_swap_strain_stress = np.column_stack((elemLabel_list,IP_coordinate_list,strain_voigt_xy_local_list,stress_voigt_xy_local_list))
merged_data_local_xz_swap_strain_stress = np.column_stack((elemLabel_list,IP_coordinate_list,strain_voigt_xz_local_list,stress_voigt_xz_local_list))
merged_data_local_yz_swap_strain_stress = np.column_stack((elemLabel_list,IP_coordinate_list,strain_voigt_yz_local_list,stress_voigt_yz_local_list))


merged_data_local = np.column_stack((elemLabel_list,IP_coordinate_list,strain_voigt_local_list,damage))
merged_data_local_xy_swap = np.column_stack((elemLabel_list,IP_coordinate_list,strain_voigt_xy_local_list,damage))
merged_data_local_xz_swap = np.column_stack((elemLabel_list,IP_coordinate_list,strain_voigt_xz_local_list,damage))
merged_data_local_yz_swap = np.column_stack((elemLabel_list,IP_coordinate_list,strain_voigt_yz_local_list,damage))



merged_data             = np.column_stack((elemLabel_list,IP_coordinate_list,strain_list,damage))
merged_data_xy_swap = np.column_stack((elemLabel_list,IP_coordinate_list,strain_voigt_xy_list,damage))
merged_data_xz_swap = np.column_stack((elemLabel_list,IP_coordinate_list,strain_voigt_xz_list,damage))
merged_data_yz_swap = np.column_stack((elemLabel_list,IP_coordinate_list,strain_voigt_yz_list,damage))


merged_data_principal   = np.column_stack((elemLabel_list,IP_coordinate_list,maxPrincipal_strain_list,midPrincipal_strain_list,minPrincipal_strain_list,damage))
merged_data_strain_stress_principal   = np.column_stack((elemLabel_list,IP_coordinate_list,maxPrincipal_strain_list,midPrincipal_strain_list,minPrincipal_strain_list,maxPrincipal_list,midPrincipal_list,minPrincipal_list))
stress_char             = np.column_stack((eta,xi))



np.savetxt("output/strain-stress_principals.csv",merged_data_strain_stress_principal,delimiter=';')
np.savetxt("output/global_strain-stress_tensor.csv",merged_data_strain_stress,delimiter=';')
np.savetxt("output/local_strain-stress_tensor.csv",merged_data_local_strain_stress,delimiter=';')
#np.savetxt("local_strain-stress_tensor_xy-swap.csv",merged_data_local_xy_swap_strain_stress,delimiter=';')
#np.savetxt("local_strain-stress_tensor_xz-swap.csv",merged_data_local_xz_swap_strain_stress,delimiter=';')
#np.savetxt("local_strain-stress_tensor_yz-swap.csv",merged_data_local_yz_swap_strain_stress,delimiter=';')

# np.savetxt("local_strain-tensor_comp_0p1.csv",merged_data_local,delimiter=';')
# np.savetxt("local_strain-tensor_xy-swap_comp_0p1.csv",merged_data_local_xy_swap,delimiter=';')
# np.savetxt("local_strain-tensor_xz-swap_comp_0p1.csv",merged_data_local_xz_swap,delimiter=';')
# np.savetxt("local_strain-tensor_yz-swap_comp_0p1.csv",merged_data_local_yz_swap,delimiter=';')

#np.savetxt("strain-tensor.csv",merged_data,delimiter=';')
#np.savetxt("strain-tensor_xy-swap.csv",merged_data_xy_swap,delimiter=';')
#np.savetxt("strain-tensor_xz-swap.csv",merged_data_xz_swap,delimiter=';')
#np.savetxt("strain-tensor_yz-swap.csv",merged_data_yz_swap,delimiter=';')


#np.savetxt("data_extracted_strain-tensor.csv",merged_data,delimiter=';')
# np.savetxt("data_extracted_strain-principals.csv",merged_data_principal,delimiter=';')
#np.savetxt("triaxility-lode.csv",stress_char,delimiter=';')

print('Data extraction finished')
