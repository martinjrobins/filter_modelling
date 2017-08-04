import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import griddata
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import glob
from collections import namedtuple
from math import pi

ParticlesData = namedtuple('ParticlesData', ['x', 'y', 'angle', 'count'])

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

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

    return ParticlesData(x,y,angle_numpy_array,count_numpy_array)

def plot_main(particles,fibres,dparticles,filename):
    print 'plotting',filename
    fig = plt.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(2, 4)
    vis_ax = plt.subplot(gs[:, 0:2])
    angle_ax = plt.subplot(gs[0, 2:4])
    xcoord_ax = plt.subplot(gs[1, 2])
    ycoord_ax = plt.subplot(gs[1, 3])

    #vis_ax.tick_params(axis='both',which='both',bottom='off',top='off',labelbottom='off')
    vis_ax.axis('off')

    vis_ax.scatter(particles.x,particles.y,2,'k',lw = 0)
    vis_ax.scatter(fibres.x,fibres.y,300,fibres.count,lw = 0)
    vis_ax.scatter(dparticles.x,dparticles.y,5,'r',lw = 0)
    vis_ax.arrow(-0.1, -0.2, 0, 1, head_width=0.1, head_length=0.2, fc='k', ec='k')
    vis_ax.arrow(-0.1, -0.2, 1, 0, head_width=0.1, head_length=0.2, fc='k', ec='k')
    vis_ax.text(-0.4,0.9,'y')
    vis_ax.text(1,-0.6,'x')
    vis_ax.set_xlim([0,10])
    vis_ax.set_ylim([-1,11])
    vis_ax.set_aspect('equal', 'datalim')


    angles = dparticles.angle-pi/2
    angles[angles<-pi] = 2*pi+angles[angles<-pi]
    angle_ax.hist(np.abs(angles),10,range=(0,pi),normed=1,color='r')
    angle_ax.set_xlabel('angle to fibre centre')
    angle_ax.set_ylabel('normed count')
    angle_ax.set_ylim([0,0.7])
    angle_ax.set_xlim([0,pi])
    angle_ax.set_xticks([0,pi/2,pi])
    angle_ax.set_xticklabels(['N','E/W','S'])

    xcoord_ax.hist(dparticles.x,10,range=(0,10),normed=1,color='r')
    xcoord_ax.set_xlabel('xcoordinate')
    xcoord_ax.set_ylabel('normed count')
    xcoord_ax.set_ylim([0,0.2])
    xcoord_ax.set_xlim([0,10])

    ycoord_ax.hist(dparticles.y,10,range=(0,10),normed=1,color='r')
    ycoord_ax.set_xlabel('ycoordinate')
    ycoord_ax.set_ylabel('normed count')
    ycoord_ax.set_ylim([0,0.5])
    ycoord_ax.set_xlim([0,10])

    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

def plot_compare(hexagon_dparticles,regular_dparticles,random_dparticles,filename):
    print 'plotting',filename

    linewidth = 3
    fig = plt.figure(figsize=(8, 8))
    gs = gridspec.GridSpec(2, 2)
    angle_ax = plt.subplot(gs[0, 0:2])
    xcoord_ax = plt.subplot(gs[1, 0])
    ycoord_ax = plt.subplot(gs[1, 1])

    angles = hexagon_dparticles.angle-pi/2
    angles[angles<-pi] = 2*pi+angles[angles<-pi]
    angle_ax.hist(np.abs(angles),10,histtype='step',range=(0,pi),normed=1,label='hexagon',lw=linewidth)
    angles = regular_dparticles.angle-pi/2
    angles[angles<-pi] = 2*pi+angles[angles<-pi]
    angle_ax.hist(np.abs(angles),10,histtype='step',range=(0,pi),normed=1,label='regular',lw=linewidth)
    angles = random_dparticles.angle-pi/2
    angles[angles<-pi] = 2*pi+angles[angles<-pi]
    angle_ax.hist(np.abs(angles),10,histtype='step',range=(0,pi),normed=1,label='random',lw=linewidth)
    angle_ax.set_xlabel('angle to fibre centre')
    angle_ax.set_ylabel('normed count')
    angle_ax.set_ylim([0,0.7])
    angle_ax.set_xlim([0,pi])
    angle_ax.set_xticks([0,pi/2,pi])
    angle_ax.set_xticklabels(['N','E/W','S'])

    xcoord_ax.hist(hexagon_dparticles.x,10,range=(0,10),histtype='step',normed=1,label='hexagon',lw=linewidth)
    xcoord_ax.hist(regular_dparticles.x,10,range=(0,10),histtype='step',normed=1,label='regular',lw=linewidth)
    xcoord_ax.hist(random_dparticles.x,10,range=(0,10),histtype='step',normed=1,label='random',lw=linewidth)
    xcoord_ax.set_xlabel('xcoordinate')
    xcoord_ax.set_ylabel('normed count')
    xcoord_ax.set_ylim([0,0.2])
    xcoord_ax.set_xlim([0,10])

    ycoord_ax.hist(hexagon_dparticles.y,10,range=(0,10),histtype='step',normed=1,label='hexagon',lw=linewidth)
    ycoord_ax.hist(regular_dparticles.y,10,range=(0,10),histtype='step',normed=1,label='regular',lw=linewidth)
    ycoord_ax.hist(random_dparticles.y,10,range=(0,10),histtype='step',normed=1,label='random',lw=linewidth)
    ycoord_ax.set_xlabel('ycoordinate')
    ycoord_ax.set_ylabel('normed count')
    ycoord_ax.set_ylim([0,0.5])
    ycoord_ax.set_xlim([0,10])
    ycoord_ax.legend(loc='upper left')

    plt.tight_layout()
    plt.savefig(filename)
    plt.close()


hexagon_dead_particles_files = glob.glob('./hexagon_dead_particles*.vtu')
regular_dead_particles_files = glob.glob('./regular_dead_particles*.vtu')
random_dead_particles_files = glob.glob('./random_dead_particles*.vtu')
for hexagon_file,regular_file,random_file,i in zip(sorted(hexagon_dead_particles_files),sorted(regular_dead_particles_files),sorted(random_dead_particles_files),range(len(hexagon_dead_particles_files))):
    print hexagon_file
    hexagon = get_particle_data(hexagon_file)
    print regular_file
    regular = get_particle_data(regular_file)
    print random_file
    random = get_particle_data(random_file)
    #plot_compare(hexagon,regular,random,'compare%05d.png'%(i))


for simtype in ['hexagon','regular','random']:
    particles_files = glob.glob('./%s_particles*.vtu'%simtype)
    fibres_files = glob.glob('./%s_fibres*.vtu'%simtype)
    dead_particles_files = glob.glob('./%s_dead_particles*.vtu'%simtype)
    for particles_file,fibres_file,dead_particles_file,i in zip(sorted(particles_files),sorted(fibres_files),sorted(dead_particles_files),range(len(particles_files))):
        print particles_file
        particles = get_particle_data(particles_file)
        print fibres_file
        fibres = get_particle_data(fibres_file)
        print dead_particles_file
        dparticles = get_particle_data(dead_particles_file)
        plot_main(particles,fibres,dparticles,'%s_main%05d.png'%(simtype,i))

