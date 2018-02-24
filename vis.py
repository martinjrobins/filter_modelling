import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from scipy.interpolate import griddata
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import glob
from collections import namedtuple
from math import pi
import xml.etree.ElementTree
import itertools

ParticlesData = namedtuple('ParticlesData', ['x', 'y', 'angle', 'count','stokes_vel_x','stokes_vel_y','electro_vel_x','electro_vel_y','charge'])

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

def get_simulation_data(filename,particle_types=['particles','fibres','dead_particles']):
     # load a xml file as input
    root = xml.etree.ElementTree.parse(filename).getroot()
    ret = []
    for particle_type in particle_types:
        particles = root.find(particle_type)
        position = particles.find('position')
        n_particles = int(position.find('count').text)
        print 'found',n_particles,'particles'
        x = np.zeros(n_particles)
        y = np.zeros(n_particles)

        for i,xy in enumerate(position.findall('item')):
            x[i] = float(xy[0][1].text)
            y[i] = float(xy[0][2].text)

        angle_xml = particles.find('angle')
        angle = np.zeros(n_particles)
        for i,a in enumerate(angle_xml.iter('item')):
            angle[i] = float(a.text);

        count_xml = particles.find('count')
        count = np.zeros(n_particles)
        for i,c in enumerate(angle_xml.iter('item')):
            count[i] = float(c.text);
        ret.append(ParticlesData(x,y,angle,count))

    return ret

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
    stokes_vtk_array = reader.GetOutput().GetPointData().GetArray("stokes_velocity")
    stokes_numpy_array = vtk_to_numpy(stokes_vtk_array)
    electro_vtk_array = reader.GetOutput().GetPointData().GetArray("electro_velocity")
    electro_numpy_array = vtk_to_numpy(electro_vtk_array)
    charge_vtk_array = reader.GetOutput().GetPointData().GetArray("charge")
    charge_numpy_array = vtk_to_numpy(charge_vtk_array)

    #Get the coordinates of the nodes and their temperatures
    nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)
    x,y,z= nodes_nummpy_array[:,0] , nodes_nummpy_array[:,1] , nodes_nummpy_array[:,2]

    return ParticlesData(x,y,angle_numpy_array,count_numpy_array,stokes_numpy_array[:,0],stokes_numpy_array[:,1],electro_numpy_array[:,0],electro_numpy_array[:,1],charge_numpy_array)


def circles(ax, x, y, s, c='b', vmin=None, vmax=None, **kwargs):
    """
    Make a scatter plot of circles.
    Similar to plt.scatter, but the size of circles are in data scale.
    Parameters
    ----------
    x, y : scalar or array_like, shape (n, )
        Input data
    s : scalar or array_like, shape (n, )
        Radius of circles.
    c : color or sequence of color, optional, default : 'b'
        `c` can be a single color format string, or a sequence of color
        specifications of length `N`, or a sequence of `N` numbers to be
        mapped to colors using the `cmap` and `norm` specified via kwargs.
        Note that `c` should not be a single numeric RGB or RGBA sequence
        because that is indistinguishable from an array of values
        to be colormapped. (If you insist, use `color` instead.)
        `c` can be a 2-D array in which the rows are RGB or RGBA, however.
    vmin, vmax : scalar, optional, default: None
        `vmin` and `vmax` are used in conjunction with `norm` to normalize
        luminance data.  If either are `None`, the min and max of the
        color array is used.
    kwargs : `~matplotlib.collections.Collection` properties
        Eg. alpha, edgecolor(ec), facecolor(fc), linewidth(lw), linestyle(ls),
        norm, cmap, transform, etc.
    Returns
    -------
    paths : `~matplotlib.collections.PathCollection`
    Examples
    --------
    a = np.arange(11)
    circles(a, a, s=a*0.2, c=a, alpha=0.5, ec='none')
    plt.colorbar()
    License
    --------
    This code is under [The BSD 3-Clause License]
    (http://opensource.org/licenses/BSD-3-Clause)
    """

    if np.isscalar(c):
        kwargs.setdefault('color', c)
        c = None

    if 'fc' in kwargs:
        kwargs.setdefault('facecolor', kwargs.pop('fc'))
    if 'ec' in kwargs:
        kwargs.setdefault('edgecolor', kwargs.pop('ec'))
    if 'ls' in kwargs:
        kwargs.setdefault('linestyle', kwargs.pop('ls'))
    if 'lw' in kwargs:
        kwargs.setdefault('linewidth', kwargs.pop('lw'))
    # You can set `facecolor` with an array for each patch,
    # while you can only set `facecolors` with a value for all.

    zipped = np.broadcast(x, y, s)
    pats = [patches.Circle((x_, y_), s_)
               for x_, y_, s_ in zipped]
    collection = PatchCollection(pats, **kwargs)
    if c is not None:
        c = np.broadcast_to(c, zipped.shape).ravel()
        collection.set_array(c)
        collection.set_clim(vmin, vmax)

    ax.add_collection(collection)
    ax.autoscale_view()
    plt.draw_if_interactive()
    #if c is not None:
        #plt.sci(collection)
    return collection

def plot_microstructure(fibres,filename):
    print 'plotting',filename
    fig = plt.figure(figsize=(6, 8))
    vis_ax = plt.gca()

    #vis_ax.tick_params(axis='both',which='both',bottom='off',top='off',labelbottom='off')
    vis_ax.axis('off')

    circles(vis_ax,fibres.x,fibres.y,0.3,'k',lw = 0)
    vis_ax.add_patch(patches.Rectangle((0,-1),10,12,fill=False))
    vis_ax.set_aspect('equal', 'datalim')
    vis_ax.set_xlim(-1,11)
    vis_ax.set_ylim(-2,12)

    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

def plot_charge(fibres,filename):
    print 'plotting',filename
    fig = plt.figure(figsize=(6, 8))
    vis_ax = plt.gca()

    #vis_ax.tick_params(axis='both',which='both',bottom='off',top='off',labelbottom='off')
    vis_ax.axis('off')

    circles(vis_ax,fibres.x,fibres.y,0.3,c=fibres.charge,lw = 0)

    vis_ax.add_patch(patches.Rectangle((0,-1),10,12,fill=False))
    vis_ax.set_aspect('equal', 'datalim')
    vis_ax.set_xlim(-1,11)
    vis_ax.set_ylim(-2,12)

    plt.tight_layout()
    plt.savefig(filename)
    plt.close()


def plot_velocity_fields(fibres,particles,filename,is_stokes):
    print 'plotting',filename
    fig = plt.figure(figsize=(6, 8))
    vis_ax = plt.gca()

    #vis_ax.tick_params(axis='both',which='both',bottom='off',top='off',labelbottom='off')
    vis_ax.axis('off')

    if is_stokes:
        circles(vis_ax,fibres.x,fibres.y,0.3,c='k',lw = 0)
    else:
        circles(vis_ax,fibres.x,fibres.y,0.3,c=fibres.charge,lw = 0)
    vis_ax.add_patch(patches.Rectangle((0,-1),10,12,fill=False))
    vis_ax.set_aspect('equal', 'datalim')
    vis_ax.set_xlim(-1,11)
    vis_ax.set_ylim(-2,12)

    X = particles.x
    Y = particles.y
    if is_stokes:
        U = particles.stokes_vel_x
        V = particles.stokes_vel_y
    else:
        U = particles.electro_vel_x
        V = particles.electro_vel_y

    C = np.sqrt(U**2 + V**2)
    vis_ax.quiver(X,Y,U,V,C)

    plt.tight_layout()
    plt.savefig(filename)
    plt.close()



def plot_main(particles,fibres,dparticles,filename,just_lhs=False):
    print 'plotting',filename
    if just_lhs:
        fig = plt.figure(figsize=(6, 8))
        vis_ax = plt.gca()
    else:
        fig = plt.figure(figsize=(12, 8))
        gs = gridspec.GridSpec(2, 4)
        vis_ax = plt.subplot(gs[:, 0:2])
        angle_ax = plt.subplot(gs[0, 2:4])
        xcoord_ax = plt.subplot(gs[1, 2])
        ycoord_ax = plt.subplot(gs[1, 3])


    #vis_ax.tick_params(axis='both',which='both',bottom='off',top='off',labelbottom='off')
    vis_ax.axis('off')

    vis_ax.scatter(particles.x,particles.y,7,'k',lw = 0)
    #vis_ax.scatter(fibres.x,fibres.y,300+70*(10-fibres.y),fibres.count,lw = 0)
    #vis_ax.scatter(fibres.x,fibres.y,300,fibres.count,lw = 0)
    circles(vis_ax,fibres.x,fibres.y,s=0.3,c=fibres.count,lw = 0)
    #vis_ax.scatter(dparticles.x,dparticles.y,7,'r',lw = 0)
    vis_ax.arrow(0.1, -1.5, 0, 1, head_width=0.1, head_length=0.2, fc='k', ec='k')
    vis_ax.arrow(0.1, -1.5, 1, 0, head_width=0.1, head_length=0.2, fc='k', ec='k')
    vis_ax.text(-0.3,-0.6,'y')
    vis_ax.text(1,-1.9,'x')
    vis_ax.set_xlim([0,10])
    vis_ax.set_ylim([-1,11])
    vis_ax.set_aspect('equal', 'datalim')

    if not just_lhs:
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

def plot_compare_elec(all_random_dparticles,values,filename):
    print 'plotting',filename

    linewidth = 3
    fig = plt.figure(figsize=(8, 8))
    gs = gridspec.GridSpec(2, 2)
    angle_ax = plt.subplot(gs[0, 0:2])
    xcoord_ax = plt.subplot(gs[1, 0])
    ycoord_ax = plt.subplot(gs[1, 1])

    for random_dparticles,value in zip(all_random_dparticles,values):
        angles = random_dparticles.angle-pi/2
        angles[angles<-pi] = 2*pi+angles[angles<-pi]
        angle_ax.hist(np.abs(angles),10,histtype='step',range=(0,pi),normed=1,label='$\phi_0=%f$'%value,lw=linewidth)

    angle_ax.set_xlabel('angle to fibre centre')
    angle_ax.set_ylabel('normed count')
    angle_ax.set_ylim([0,0.7])
    angle_ax.set_xlim([0,pi])
    angle_ax.set_xticks([0,pi/2,pi])
    angle_ax.set_xticklabels(['N','E/W','S'])

    for random_dparticles,value in zip(all_random_dparticles,values):
        xcoord_ax.hist(random_dparticles.x,10,range=(0,10),histtype='step',normed=1,label='$\phi_0=%f$'%value,lw=linewidth)

    xcoord_ax.set_xlabel('xcoordinate')
    xcoord_ax.set_ylabel('normed count')
    xcoord_ax.set_ylim([0,0.2])
    xcoord_ax.set_xlim([0,10])

    for random_dparticles,value in zip(all_random_dparticles,values):
        ycoord_ax.hist(random_dparticles.y,10,range=(0,10),histtype='step',normed=1,label='$\phi_0=%f$'%value,lw=linewidth)

    ycoord_ax.set_xlabel('ycoordinate')
    ycoord_ax.set_ylabel('normed count')
    ycoord_ax.set_ylim([0,0.5])
    ycoord_ax.set_xlim([0,10])
    ycoord_ax.legend(loc='upper left')

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

random_sims = ['%02drandom'%i for i in range(0,21)]
hexagon_dead_particles_files = sorted(glob.glob('./hexagon_dead_particles*.vtu'))
regular_dead_particles_files = sorted(glob.glob('./regular_dead_particles*.vtu'))
hexagon_fibre_files = sorted(glob.glob('./hexagon_fibre*.vtu'))
regular_fibre_files = sorted(glob.glob('./regular_fibre*.vtu'))
random_fibre_files = sorted(glob.glob('./00random_fibre*.vtu'))

regular_simulation_files = sorted(glob.glob('./regular_simulation*.xml'))
hexagon_simulation_files = sorted(glob.glob('./hexagon_simulation*.xml'))

random_simulation_files_tmp = []
random_dead_particles_files_tmp = []
for i in random_sims:
    random_simulation_files_tmp.append(sorted(glob.glob('./%s_simulation*.xml'%i)))
    random_dead_particles_files_tmp.append(sorted(glob.glob('./%s_dead_particles*.vtu'%i)))

random_simulation_files = []
random_dead_particles_files = []
for tpl in itertools.izip(*random_simulation_files_tmp):
    random_simulation_files.append(tpl)
for tpl in itertools.izip(*random_dead_particles_files_tmp):
    random_dead_particles_files.append(tpl)

print hexagon_dead_particles_files
print random_dead_particles_files


if len(hexagon_fibre_files) > 0:
    plot_microstructure(get_particle_data(hexagon_fibre_files[0]),'microstructureH.pdf')
    plot_charge(get_particle_data(hexagon_fibre_files[0]),'microstructure_charge_H.pdf')
if len(regular_fibre_files) > 0:
    plot_microstructure(get_particle_data(regular_fibre_files[0]),'microstructureS.pdf')
    plot_charge(get_particle_data(regular_fibre_files[0]),'microstructure_charge_S.pdf')
if len(random_fibre_files) > 0:
    plot_microstructure(get_particle_data(random_fibre_files[0]),'microstructureR.pdf')
    plot_charge(get_particle_data(random_fibre_files[0]),'microstructure_charge_R.pdf')

elec_dirs = ['electro_0_0.0', 'electro_1_0.04','electro_1_0.08','electro_1_0.1']
elec_value = [0.0,0.04,0.08,0.1]
elec_random = []
elec_hex = []
elec_reg = []
#for elec_dir in elec_dirs:
#    print elec_dir
#    hexagon = get_particle_data('./%s/hexagon_dead_particles19800.vtu'%elec_dir)
#    elec_hex.append(hexagon)
#    regular = get_particle_data('./%s/regular_dead_particles19800.vtu'%elec_dir)
#    elec_reg.append(regular)
#    random_dead_particles_files_tmp = []
#    random_files = sorted(glob.glob('./%s/*random_dead_particles19800.vtu'%elec_dir))
#    random = ParticlesData(np.zeros(0),np.zeros(0),np.zeros(0),np.zeros(0),np.zeros(0),np.zeros(0),np.zeros(0),np.zeros(0),np.zeros(0))
#    for random_file in random_files:
#        print random_file
#        tmp = get_particle_data(random_file)
#        random = ParticlesData(np.concatenate([random.x,tmp.x]),
#                                 np.concatenate([random.y,tmp.y]),
#                                 np.concatenate([random.angle,tmp.angle]),
#                                 np.concatenate([random.count,tmp.count]),
#                                 np.concatenate([random.stokes_vel_x,tmp.stokes_vel_x]),
#                                 np.concatenate([random.stokes_vel_y,tmp.stokes_vel_y]),
#                                 np.concatenate([random.electro_vel_x,tmp.electro_vel_x]),
#                                 np.concatenate([random.electro_vel_y,tmp.electro_vel_y]),
#                                 np.concatenate([random.charge,tmp.charge])
#                                 )
#    elec_random.append(random)
#
#plot_compare_elec(elec_random,elec_value,'compare_elec_random.pdf')
#plot_compare_elec(elec_hex,elec_value,'compare_elec_hexagon.pdf')
#plot_compare_elec(elec_reg,elec_value,'compare_elec_regular.pdf')



#for hexagon_file,regular_file,random_files,i in zip(hexagon_dead_particles_files,regular_dead_particles_files,random_dead_particles_files,range(len(hexagon_dead_particles_files))):
#    print hexagon_file
#    hexagon = get_particle_data(hexagon_file)
#    print regular_file
#    regular = get_particle_data(regular_file)
#    random = ParticlesData(np.zeros(0),np.zeros(0),np.zeros(0),np.zeros(0),np.zeros(0),np.zeros(0),np.zeros(0),np.zeros(0),np.zeros(0))
#    for random_file in random_files:
#        print random_file
#        tmp = get_particle_data(random_file)
#        random = ParticlesData(np.concatenate([random.x,tmp.x]),
#                                 np.concatenate([random.y,tmp.y]),
#                                 np.concatenate([random.angle,tmp.angle]),
#                                 np.concatenate([random.count,tmp.count]),
#                                 np.concatenate([random.stokes_vel_x,tmp.stokes_vel_x]),
#                                 np.concatenate([random.stokes_vel_y,tmp.stokes_vel_y]),
#                                 np.concatenate([random.electro_vel_x,tmp.electro_vel_x]),
#                                 np.concatenate([random.electro_vel_y,tmp.electro_vel_y]),
#                                 np.concatenate([random.charge,tmp.charge])
#                                 )
#    plot_compare(hexagon,regular,random,'compare%05d.pdf'%(i))

#for hexagon_file,regular_file,random_files,i in zip(hexagon_simulation_files,regular_simulation_files,random_simulation_files,range(len(hexagon_simulation_files))):
#    print hexagon_file
#    hexagon, = get_simulation_data(hexagon_file,['dead_particles'])
#    print regular_file
#    regular, = get_simulation_data(regular_file,['dead_particles'])
#    print random_files
#    random = ParticlesData(np.zeros(0),np.zeros(0),np.zeros(0),np.zeros(0))
#    for random_file in random_files:
#        tmp, = get_simulation_data(random_file,['dead_particles'])
#        random = ParticlesData(np.concatenate([random.x,tmp.x]),
#                                 np.concatenate([random.y,tmp.y]),
#                                 np.concatenate([random.angle,tmp.angle]),
#                                 np.concatenate([random.count,tmp.count]))
#    plot_compare(hexagon,regular,random,'compare%05d.pdf'%(i))



for simtype in ['hexagon','regular']+random_sims:
#for simtype in ['hexagon']+random_sims:
#for simtype in ['regular']:
    particles_files = sorted(glob.glob('./%s_particles19800.vtu'%simtype))
    fibres_files = sorted(glob.glob('./%s_fibres19800.vtu'%simtype))
    dead_particles_files = sorted(glob.glob('./%s_dead_particles19800.vtu'%simtype))
    for particles_file,fibres_file,dead_particles_file,i in zip(particles_files,fibres_files,dead_particles_files,range(len(particles_files))):
        print particles_file
        particles = get_particle_data(particles_file)
        print fibres_file
        fibres = get_particle_data(fibres_file)
        print dead_particles_file
        dparticles = get_particle_data(dead_particles_file)
        plot_velocity_fields(fibres,particles,'%s_velocity_stokes%05d.pdf'%(simtype,i),True)
        plot_velocity_fields(fibres,particles,'%s_velocity_electro%05d.pdf'%(simtype,i),False)
        plot_main(particles,fibres,dparticles,'%s_main%05d.pdf'%(simtype,i))
        plot_main(particles,fibres,dparticles,'%s_main_just_lhs%05d.pdf'%(simtype,i),True)

#for simtype in ['hexagon','regular']+random_sims:
#for simtype in ['hexagon']+random_sims:
#    sim_files = sorted(glob.glob('./%s_simulation*.xml'%simtype))
#    for sim_file,i in zip(sim_files,range(len(sim_files))):
#        print sim_file
#        particles,fibres,dparticles = get_simulation_data(sim_file)
#        plot_main(particles,fibres,dparticles,'%s_main%05d.pdf'%(simtype,i))
#        plot_main(particles,fibres,dparticles,'%s_main_just_lhs%05d.pdf'%(simtype,i),True)

