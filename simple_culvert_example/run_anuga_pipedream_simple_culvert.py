#------------------------------------------------------------------------------
# IMPORT NECESSARY MODULES
#------------------------------------------------------------------------------
print (' ABOUT to Start Simulation:- Importing Modules')


import anuga
import numpy as np

#------------------------------------------------------------------------------
# FILENAMES, MODEL DOMAIN and VARIABLES
#------------------------------------------------------------------------------

basename = 'simple_culvert'
outname =  'anuga_pipedream_simple_culvert'

rf = 2  # refinement factor for domain

dt = 0.01     # yield step
out_dt = 2.0 # output step
ft = 400     # final timestep

verbose = False

#------------------------------------------------------------------------------
# CREATING MESH
#------------------------------------------------------------------------------

domain = anuga.rectangular_cross_domain(60*rf,20*rf, len1=60, len2=20)

#------------------------------------------------------------------------------
# SETUP COMPUTATIONAL DOMAIN
#------------------------------------------------------------------------------
domain.set_minimum_storable_height(0.0001) 
domain.set_name(outname) 
print (domain.statistics())


#------------------------------------------------------------------------------
# APPLY ELEVATION
#------------------------------------------------------------------------------
def topography(x,y):

    z = 5*np.ones_like(x)

    channel = np.logical_and(y>5,y<15)

    z = np.where(np.logical_and(channel,x<10), x/300, z)
    z = np.where(np.logical_and(channel,x>20), x/300, z)

    return z

domain.set_quantity('elevation', topography, location='centroids')

#------------------------------------------------------------------------------
# APPLY MANNING'S ROUGHNESSES
#------------------------------------------------------------------------------

domain.set_quantity('friction', 0.035)

#------------------------------------------------------------------------------
# SETUP BOUNDARY CONDITIONS
#------------------------------------------------------------------------------

print ('Available boundary tags', domain.get_boundary_tags())

Br = anuga.Reflective_boundary(domain)
Bd = anuga.Dirichlet_boundary([-1.0,0,0])

domain.set_boundary({'left': Bd, 'bottom': Br, 'top': Br, 'right': Br})


inlet1_anuga_region = anuga.Region(domain, polygon=[[20.0,5.0], [22.0, 5.0], [22.0, 15.0], [20.0, 15.0]])
outlet_anuga_region = anuga.Region(domain, polygon=[[8.0,5.0], [10.0, 5.0], [10.0, 15.0], [8.0, 15.0]])

anuga_Lweirs = np.array([20.0, 20.0])
anuga_Amanholes = np.array([20.0, 20.0])

inlet1_anuga_inlet_op = anuga.Inlet_operator(domain, inlet1_anuga_region, Q=0.0, zero_velocity=True)
outlet_anuga_inlet_op = anuga.Inlet_operator(domain, outlet_anuga_region, Q=0.0, zero_velocity=False)

anuga_beds = np.array([inlet1_anuga_inlet_op.inlet.get_average_elevation(),
                       outlet_anuga_inlet_op.inlet.get_average_elevation()])

print(anuga_beds)


input_Q = 1.0
line=[[59.0, 5.0],[59.0, 15.0]]
anuga.Inlet_operator(domain, line, input_Q)

#------------------------------------------------------------------------------
# PIPEDREAM
#------------------------------------------------------------------------------

print('Setup pipedream structures')
from pipedream_solver.hydraulics import SuperLink
import matplotlib.pyplot as plt
import pandas as pd


# Details availabe from https://mattbartos.com/pipedream/geometry-reference.html and
# https://github.com/mdbartos/pipedream/blob/master/pipedream_solver/geometry.py

superjunctions = pd.DataFrame({'name': [0, 1],
                               'id': [0, 1],
                               'z_inv': [0.04, 0.00],
                               'h_0': 2*[0],
                               'bc': 2*[False],
                               'storage': 2*['functional'],
                               'a': 2*[0.],
                               'b': 2*[1.],
                               'c': 2*[10.],
                               'max_depth': 2*[np.inf],
                               'map_x': 2*[0],
                               'map_y': 2*[0]})

superlinks = pd.DataFrame({'name': [0],
                           'id': [0],
                           'sj_0': [0],
                           'sj_1': [1],
                           'in_offset': 1*[0.],
                           'out_offset': 1*[0.],
                           'dx': [10],
                           'n': 1*[0.013],
                           'shape': 1*['rect_closed'],
                           'g1': 1*[1.0],
                           'g2': 1*[10.0],
                           'g3': 1*[0.1],
                           'g4': 1*[0.],
                           'Q_0': 1*[0.],
                           'h_0': 1*[1e-5],
                           'ctrl': 1*[False],
                           'A_s': 1*[0.],
                           'A_c': 1*[0.],
                           'C': 1*[0.]})

superlink = SuperLink(superlinks, superjunctions, internal_links=6)


#--------------------------------------------------------------------------
# Setup storage for output
#--------------------------------------------------------------------------
H_js = []
losses = []

Q_iks =[]
Q_uks =[]
Q_dks =[]
time_series = []
anuga_ws = []
Q_ins = []



def Calculate_Q(head1D, depth2D, bed2D, Lweir, Amanhole, cw=0.67, co=0.67):
    """
    Reference:
    A methodology for linking 2D overland flow models with the sewer network model SWMM 5.1 
    based on dynamic link libraries
    Leandro, Jorge, Martins, Ricardo
    Water Science and Technology
    2016, 73, 3017-3026

    cw is the weir discharge coefficient, 
    w is the weir crest width [m], 
    Amh is the manhole area [m2] 
    co is the orifice discharge coefficient.

    """

    import numpy as np
    from anuga import g

    Q = np.zeros_like(head1D)

    # if head1D < bed2D use Weir Equation (Reference Equation (10)):
    Q = np.where(head1D<bed2D, cw*Lweir*depth2D*np.sqrt(2*g*depth2D), Q)

    # If head1D > bed2D and  head1D < (depth2D + bed2d) use orifice equation (Equation (11))
    Q = np.where(np.logical_and(bed2D<=head1D, head1D<depth2D+bed2D) , co*Amanhole*np.sqrt(2*g*(depth2D+bed2D-head1D)), Q)

    # Otherwise if h1d >= h2d + Z2d use orifice equation (Equation (11)) surcharge
    Q = np.where(head1D>=depth2D+bed2D,  -co*Amanhole*np.sqrt(2*g*(head1D-depth2D-bed2D)), Q)

    return Q

#---------------------------------------------------------------------------
print('Start Evolve')
#---------------------------------------------------------------------------

# slow the response of the coupling calculation
time_average = 10 # sec
Q_in_old = np.array([0.0, 0.0])

for t in domain.evolve(yieldstep=dt, outputstep=out_dt, finaltime=ft):
    #print('\n')
    if domain.yieldstep_counter%domain.output_frequency == 0:
        domain.print_timestepping_statistics()

    anuga_depths = np.array([inlet1_anuga_inlet_op.inlet.get_average_depth(),
                             outlet_anuga_inlet_op.inlet.get_average_depth()])
    
    anuga_stages = np.array([inlet1_anuga_inlet_op.inlet.get_average_stage(),
                             outlet_anuga_inlet_op.inlet.get_average_stage()])


    # Compute the water volumes
    link_volume = ((superlink._A_ik * superlink._dx_ik).sum() +
                   (superlink._A_SIk * superlink._h_Ik).sum())
    node_volume = (superlink._A_sj * (superlink.H_j - superlink._z_inv_j)).sum()
    sewer_volume = link_volume + node_volume

    boundary_flux = domain.get_boundary_flux_integral()
    total_volume_correct = t*input_Q + boundary_flux 
    
    total_volume_real = domain.get_water_volume() + sewer_volume
    loss = total_volume_real - total_volume_correct

    if domain.yieldstep_counter%domain.output_frequency == 0:
        print('    Loss         ', loss)
        print('    TV correct   ', total_volume_correct)
        print('    domain volume', domain.get_water_volume())
        print('    node_volume  ', node_volume)
        print('    sewer_volume ', sewer_volume)
        print('    anuga_depths ', anuga_depths)
        print('    anuga_beds   ', anuga_beds)
        print('    MOInvert     ', superlink._z_inv_j)
        print('    Head         ', superlink.H_j)
        print('    anuga_stages ', anuga_stages)

    # Append data
    time_series.append(t)
    losses.append(loss)
    H_js.append(superlink.H_j.copy())
    anuga_ws.append(anuga_stages.copy())
    

    # record flow time series in each pipe
    Q_iks.append(superlink.Q_ik.copy())
    Q_uks.append(superlink.Q_uk.copy())
    Q_dks.append(superlink.Q_dk.copy())

        
    # Calculate discharge at inlets and smooth its response
    Q_in = Calculate_Q(superlink.H_j, anuga_depths, anuga_beds, anuga_Lweirs, anuga_Amanholes)

    Q_in = ((time_average - dt)*Q_in_old + dt*Q_in)/time_average
    Q_in_old = Q_in

    Q_ins.append(Q_in.copy())

    if domain.yieldstep_counter%domain.output_frequency == 0:
        print('    Q            ', Q_in)
    
    # Simulate sewer with flow input (or outflow)
    superlink.step(Q_in=Q_in, dt=dt)
    superlink.reposition_junctions()

    # And consequently set anuga sim with flow output (or inflow)
    inlet1_anuga_inlet_op.set_Q(-Q_in[0])
    outlet_anuga_inlet_op.set_Q(-Q_in[1])


H_j = np.vstack(H_js)
anuga_j = np.vstack(anuga_ws)
Q_ins = np.vstack(Q_ins)

plt.ion()

plt.figure(1)
plt.plot(time_series, H_j[:,0], label='Pipe Inlet 0')
plt.plot(time_series, H_j[:,1], label='Pipe Inlet 1')
plt.plot(time_series, anuga_j[:,0], label='Anuga Inlet 0')
plt.plot(time_series, anuga_j[:,1], label='Anuga Inlet 1')
plt.legend()
plt.title('Head at junctions')
plt.xlabel('Time (s)')
plt.ylabel('Head (m)')
plt.savefig('Figure1.png')
plt.show()

plt.figure(2)
plt.clf()
plt.plot(time_series, losses)
plt.title('Losses')
plt.savefig('Figure2.png')
plt.show()

plt.figure(3)
plt.clf()
plt.plot(time_series, Q_dks)
plt.title('Q_dks')
plt.savefig('Figure3.png')
plt.show()

plt.figure(4)
plt.clf()
plt.plot(time_series, Q_ins[:,0], label='Inlet 0')
plt.plot(time_series, Q_ins[:,1], label='Inlet 1')
plt.legend()
plt.title('Q_in')
plt.savefig('Figure4.png')
plt.show()

input('Enter key ...')
