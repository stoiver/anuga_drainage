#------------------------------------------------------------------------------
print('IMPORT NECESSARY MODULES')
#------------------------------------------------------------------------------

import anuga
import numpy as np

#------------------------------------------------------------------------------
print('FILENAMES, MODEL DOMAIN and VARIABLES')
#------------------------------------------------------------------------------

basename = 'simple_culvert'
outname =  'anuga_pipedream_simple_culvert'


rf = 20  # refinement factor for domain, if too coarse the inlets will overlap the wall

dt = 0.05     # yield step
out_dt = 1.0 # output step
ft = 400     # final timestep

verbose = False


#------------------------------------------------------------------------------
print('SETUP COMPUTATIONAL DOMAIN')
#------------------------------------------------------------------------------
domain = anuga.rectangular_cross_domain(3*rf, rf, len1=60, len2=20)
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
print('SETUP BOUNDARY CONDITIONS')
#------------------------------------------------------------------------------

print ('Available boundary tags', domain.get_boundary_tags())

Br = anuga.Reflective_boundary(domain)
Bd = anuga.Dirichlet_boundary([-1.0,0,0])

domain.set_boundary({'left': Bd, 'bottom': Br, 'top': Br, 'right': Br})


#------------------------------------------------------------------------------
print('SETUP ANUGA INFLOW INLET')
#------------------------------------------------------------------------------

input_Q = 1.0
line=[[59.0, 5.0],[59.0, 15.0]]
anuga.Inlet_operator(domain, line, input_Q)

#------------------------------------------------------------------------------
print('SETUP ANUGA INLETS FOR COUPLING')
#------------------------------------------------------------------------------

print('Setup anuga inlets for coupling')
inlet1_anuga_region = anuga.Region(domain, polygon=[[20.0,5.0], [22.0, 5.0], [22.0, 15.], [20.0, 15.0]])
outlet_anuga_region = anuga.Region(domain, polygon=[[8.0,5.0], [10.0, 5.0], [10.0, 15.0], [8.0, 15.0]])

anuga_length_weirs = np.array([20.0, 20.0])
anuga_area_manholes = np.array([20.0, 20.0])

inlet1_anuga_inlet_op = anuga.Inlet_operator(domain, inlet1_anuga_region, Q=0.0, zero_velocity=True)
outlet_anuga_inlet_op = anuga.Inlet_operator(domain, outlet_anuga_region, Q=0.0, zero_velocity=False)

anuga_beds = np.array([inlet1_anuga_inlet_op.inlet.get_average_elevation(),
                       outlet_anuga_inlet_op.inlet.get_average_elevation()])

print('anuga beds', anuga_beds)

#------------------------------------------------------------------------------
print('Setup PIPEDREAM')
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
print('Setup storage for output')
#--------------------------------------------------------------------------
H_js = []
losses = []

Q_iks =[]
Q_uks =[]
Q_dks =[]
time_series = []
anuga_ws = []
Q_ins = []


from coupling import calculate_Q

#---------------------------------------------------------------------------
print('Set time averaging of Q')
#---------------------------------------------------------------------------

# slow the response of the coupling calculation
time_average = 1 # sec
Q_in_old = np.array([0.0, 0.0])

#---------------------------------------------------------------------------
print('Start Evolve')
#---------------------------------------------------------------------------
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
    Q_in = calculate_Q(superlink.H_j, anuga_depths, anuga_beds, anuga_length_weirs, anuga_area_manholes)

    Q_in = ((time_average - dt)*Q_in_old + dt*Q_in)/time_average
    Q_in_old = Q_in

    Q_ins.append(Q_in.copy())

    if domain.yieldstep_counter%domain.output_frequency == 0:
        print('    Q            ', Q_in)
    
    # Simulate sewer with flow input (or outflow)
    superlink.step(Q_in=Q_in, dt=dt)
    #superlink.reposition_junctions()

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
