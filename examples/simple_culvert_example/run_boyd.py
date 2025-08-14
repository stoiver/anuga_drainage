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
outname =  'domain_boyd'

rf = 20  # refinement factor for domain, if too coarse the inlets will overlap the wall

dt = 20.0     # yield step
out_dt = dt   # output step
ft = 400      # final timestep

verbose = False

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
# SETUP BOUNDARY CONDITIONS
#------------------------------------------------------------------------------

print ('Available boundary tags', domain.get_boundary_tags())

Br = anuga.Reflective_boundary(domain)
Bd = anuga.Dirichlet_boundary([-1.0,0,0])

domain.set_boundary({'left': Bd, 'bottom': Br, 'top': Br, 'right': Br})


#------------------------------------------------------------------------------
# SETUP ANUGA INLETS FOR COUPLING
#------------------------------------------------------------------------------

print('Setup anuga inflow inlet')

input_Q = 10.0
line=[[59.0, 7.0],[59.0, 13.0]]
polygon = [[59.0, 7.0],[59.0, 13.0], [50.0, 13.0], [50.0, 7.0]]
inlet_region = anuga.Region(domain, polygon= polygon)
anuga.Inlet_operator(domain, inlet_region, input_Q)

#------------------------------------------------------------------------------
# BOYD PIPE CULVERT
#------------------------------------------------------------------------------

print('Setup boyd culvert')

losses = {'inlet':0.0, 'outlet':0.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
ep0 = np.array([21, 10.0]) # upstream
ep1 = np.array([9,  10.0]) # downstream


el0 = [ [21.0, 5.5], [21.0, 14.5]]
el1 = [ [9.0, 5.5], [9.0, 14.5]]

eq0 = [23.0,10.0]
eq1 = [7.0, 10.0]

ep0_tid = domain.get_triangle_containing_point(ep0)
ep1_tid = domain.get_triangle_containing_point(ep1)

eq0_tid = domain.get_triangle_containing_point(eq0)
eq1_tid = domain.get_triangle_containing_point(eq1)

stage_c = domain.get_quantity('stage').centroid_values
elev_c  = domain.get_quantity('elevation').centroid_values
xmom_c  = domain.get_quantity('xmomentum').centroid_values

    
invert_elevations=[0.07, 0.03]  

culvert = anuga.Boyd_box_operator(domain,
    losses=losses,
    width=10.0,
    height=1.0,
    #end_points=[ep0, ep1],
    exchange_lines=[el0,el1],
    enquiry_points=[eq0,eq1],
    invert_elevations=invert_elevations,
    use_momentum_jet=True,
    use_velocity_head=True,
    manning=0.01,
    logging=False,
    label='boyd_box', 
    verbose=False)
    
print(culvert.statistics())
    
accumulated_flow = culvert.accumulated_flow

#------------------------------------------------------------------------------
# EVOLVE SYSTEM THROUGH TIME
#------------------------------------------------------------------------------

import time
t0 = time.time()

for t in domain.evolve(yieldstep=dt, outputstep=out_dt, finaltime=ft):
    #print('\n')
    if domain.yieldstep_counter%domain.output_frequency == 0:
        domain.print_timestepping_statistics()

    anuga_stages = np.array([stage_c[eq1_tid], stage_c[eq0_tid]])                 
    anuga_depths = np.array([stage_c[eq1_tid]-elev_c[eq1_tid], stage_c[eq0_tid]-elev_c[eq0_tid]])    
    anuga_beds   = np.array([elev_c[eq1_tid], elev_c[eq0_tid]]) 
    anuga_xmom   = np.array([xmom_c[eq1_tid], xmom_c[eq0_tid]]) 

    sewer_volume = 0.0

    boundary_flux = domain.get_boundary_flux_integral()
    total_volume_correct = t*input_Q + boundary_flux 
    
    total_volume_real = domain.get_water_volume() + sewer_volume
    loss = total_volume_real - total_volume_correct

    old_accumulated_flow = accumulated_flow
    accumulated_flow = culvert.accumulated_flow

    yield_step_flow = (accumulated_flow - old_accumulated_flow)/dt

    if domain.yieldstep_counter%domain.output_frequency == 0:
        print('    Loss (m^3)            ', loss)
        print('    TV correct (m^3)      ', total_volume_correct)
        print('    domain volume (m^3)   ', domain.get_water_volume())
        print('    sewer_volume (m^3)    ', sewer_volume)
        print('    anuga_beds (m)        ', anuga_beds)
        print('    anuga_stages (m)      ', anuga_stages)
        print('    anuga_depths (m)      ', anuga_depths)
        print('    anuga_xmom (m^2/s)    ', anuga_xmom) 
        print('    culvert vel (m/s)     ', culvert.velocity)
        print('    culvert depth (m)     ', culvert.outlet_depth)
        print('    acc flow rate (m^3/s) ', yield_step_flow)

        #print('time,ins inst Q,aver Q,vel,driving E,delta_total_E')
        
 

print ('Finished')
