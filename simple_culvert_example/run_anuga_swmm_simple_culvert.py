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
outname =  'anuga_swmm_simple_culvert'

rf = 20  # refinement factor for domain, if too coarse the inlets will overlap the wall

dt = 1.0     # yield step
out_dt = 2.0 # output step
ft = 10      # final timestep

verbose   = False
visualise = False

#------------------------------------------------------------------------------
# CREATING MESH
#------------------------------------------------------------------------------

domain = anuga.rectangular_cross_domain(3*rf, rf, len1=60, len2=20)

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
# SETUP ANUGA INLETS FOR COUPLING
#------------------------------------------------------------------------------

print('Setup anuga inflw inlet')

input_Q = 1.0
line=[[59.0, 5.0],[59.0, 15.0]]
anuga.Inlet_operator(domain, line, input_Q)

#------------------------------------------------------------------------------
# SETUP BOUNDARY CONDITIONS
#------------------------------------------------------------------------------

print ('Available boundary tags', domain.get_boundary_tags())

Br = anuga.Reflective_boundary(domain)
Bd = anuga.Dirichlet_boundary([-1.0,0,0])

domain.set_boundary({'left': Bd, 'bottom': Br, 'top': Br, 'right': Br})


inlet1_anuga_region = anuga.Region(domain, polygon=[[20.0,5.0], [22.0, 5.0], [22.0, 15.0], [20.0, 15.0]])
outlet_anuga_region = anuga.Region(domain, polygon=[[8.0,5.0], [10.0, 5.0], [10.0, 15.0], [8.0, 15.0]])

anuga_length_weirs = np.array([20.0, 20.0])
anuga_area_manholes = np.array([20.0, 20.0])

inlet1_anuga_inlet_op = anuga.Inlet_operator(domain, inlet1_anuga_region, Q=0.0, zero_velocity=True)
outlet_anuga_inlet_op = anuga.Inlet_operator(domain, outlet_anuga_region, Q=0.0, zero_velocity=False)

anuga_beds = np.array([inlet1_anuga_inlet_op.inlet.get_average_elevation(),
                       outlet_anuga_inlet_op.inlet.get_average_elevation()])

print(anuga_beds)


#------------------------------------------------------------------------------
# SWMM
#------------------------------------------------------------------------------

print('Setup swmm simulation')
from pyswmm import Simulation, Nodes, Links
import matplotlib.pyplot as plt
import pandas as pd

sim = Simulation('./swmm_input.inp')
sim.start()

swmm_inlet = Nodes(sim)['Inlet']
swmm_outlet = Nodes(sim)['Outlet']
swmm_outfall = Nodes(sim)['Outfall']
swmm_culvert = Links(sim)['Culvert']
swmm_outpipe = Links(sim)['Outpipe']

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


from coupling import calculate_Q

#---------------------------------------------------------------------------
print('Start Evolve')
#---------------------------------------------------------------------------

# slow the response of the coupling calculation
time_average = dt # sec
Q_in_old = np.array([0.0, 0.0])

for t in domain.evolve(yieldstep=dt, outputstep=dt, finaltime=ft):
    #print('\n')
    if domain.yieldstep_counter%domain.output_frequency == 0:
        domain.print_timestepping_statistics()

    anuga_depths = np.array([inlet1_anuga_inlet_op.inlet.get_average_depth(),
                             outlet_anuga_inlet_op.inlet.get_average_depth()])
    
    anuga_stages = np.array([inlet1_anuga_inlet_op.inlet.get_average_stage(),
                             outlet_anuga_inlet_op.inlet.get_average_stage()])


    # FIXME: Compute the water volumes in the swmm model
    link_volume = 0.0
    node_volume = 0.0
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
        print('    anuga_stages ', anuga_stages)

    # Append data
    time_series.append(t)
    losses.append(loss)
    anuga_ws.append(anuga_stages.copy())
    

    # FIXME: record flow time series in each pipe
    inlet_head  = swmm_inlet.head
    outlet_head = swmm_outlet.head
    outfall_head = swmm_outfall.head

    inlet_invert = swmm_inlet.invert_elevation
    outlet_invert = swmm_outlet.invert_elevation
    outfall_invert = swmm_outfall.invert_elevation

    if domain.yieldstep_counter%domain.output_frequency == 0:
        print('    Inlet Head   ', inlet_head)
        print('    Outlet Head  ', outlet_head)
        print('    Outfall Head ', outfall_head)
        print('    Inlet invert   ', inlet_invert)
        print('    Outlet invert  ', outlet_invert)
        print('    Outfall invert ', outfall_invert)

    node_heads = np.array([inlet_head, outlet_head])

        
    # Calculate discharge at inlets and smooth its response
    Q_in = calculate_Q(node_heads, anuga_depths, anuga_beds, anuga_length_weirs, anuga_area_manholes)

    Q_in = ((time_average - dt)*Q_in_old + dt*Q_in)/time_average
    Q_in_old = Q_in

    Q_ins.append(Q_in.copy())

    if domain.yieldstep_counter%domain.output_frequency == 0:
        print('    Q            ', Q_in)
    
    # Simulate sewer with flow input (or outflow)
    swmm_inlet.generated_inflow(Q_in[0])
    swmm_outlet.generated_inflow(Q_in[1])
    sim.step_advance(dt)
    #print(sim.current_time)
    sim.next()

    # determine how much actually flowed into 1D model

    if domain.yieldstep_counter%domain.output_frequency == 0:    
        print('    Inlet flooding ', swmm_inlet.flooding)
        print('    Inlet depth    ', swmm_inlet.depth)
        print('    Inlet volume   ', swmm_inlet.volume)
        print('    Inlet surcharge', swmm_inlet.surcharge_depth)

        print('    Inlet lat inflow', swmm_inlet.lateral_inflow)
        print('    Inlet gen inflow', Q_in[0])
        print('    Inlet tot inflow', swmm_inlet.total_inflow)
        print('    Inlet tot outflow', swmm_inlet.total_outflow)

    # And consequently set anuga sim with flow output (or inflow)
    inlet1_anuga_inlet_op.set_Q(-Q_in[0])
    outlet_anuga_inlet_op.set_Q(-Q_in[1])


sim.close()

if visualise:
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
