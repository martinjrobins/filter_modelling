import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import glob

def get_particle_data(filename):
     # load a vtk file as input
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()

    # Get the coordinates of nodes in the mesh
    nodes_vtk_array= reader.GetOutput().GetPoints().GetData()
    #The "Temperature" field is the third scalar in my vtk file
    angle_vtk_array = reader.GetOutput().GetPointData().GetArray("angle")
    angle_numpy_array = vtk_to_numpy(angle_vtk_array)
    count_vtk_array = reader.GetOutput().GetPointData().GetArray("count")
    count_numpy_array = vtk_to_numpy(count_vtk_array)

#Get the coordinates of the nodes and their temperatures
    nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)
    x,y,z= nodes_nummpy_array[:,0] , nodes_nummpy_array[:,1] , nodes_nummpy_array[:,2]

    return x,y,z,angle_numpy_array,count_numpy_array

for simtype in ['regular','hexagon','random']:
    particles_files = glob.glob('./%s_particles*.vtu'%simtype)
    fibres_files = glob.glob('./%s_fibres*.vtu'%simtype)
    dead_particles_files = glob.glob('./%s_dead_particles*.vtu'%simtype)
    for particles_file,fibres_file,dead_particles_file in zip(sorted(particles_files),sorted(fibres_files),sorted(dead_particles_files)):
        print particles_file
        x,y,z,angle,count = get_particle_data(particles_file)
        print fibres_file
        x,y,z,angle,count = get_particle_data(fibres_file)
        print dead_particles_file
        x,y,z,angle,count = get_particle_data(dead_particles_file)

