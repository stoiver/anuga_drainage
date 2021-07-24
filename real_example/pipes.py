#------------------------------------------------------------------------------
# IMPORT NECESSARY MODULES
#------------------------------------------------------------------------------
print (' ABOUT to Start Simulation:- Importing Modules')

import anuga, numpy, time, os, glob
from anuga import file_function, Polygon_function, read_polygon, create_mesh_from_regions, Domain, Inlet_operator
import anuga.utilities.spatialInputUtil as su

from anuga import distribute, myid, numprocs, finalize, barrier
from anuga import Inlet_operator, Boyd_box_operator, Boyd_pipe_operator
from anuga import Rate_operator
from anuga import Region

import numpy as np
import math

#------------------------------------------------------------------------------
# FILENAMES, MODEL DOMAIN and VARIABLES
#------------------------------------------------------------------------------

basename = 'model/terrain'
outname = 'pipes'
meshname = 'model/terrain.tsh'

#------------------------------------------------------------------------------
# CREATING MESH
#------------------------------------------------------------------------------
CatchmentDictionary = {'model/kerb/kerb1.csv':0.01, 'model/kerb/kerb2.csv':0.01}
    
bounding_polygon = anuga.read_polygon('model/domain.csv')
interior_regions = anuga.read_polygon_dir(CatchmentDictionary, 'model/kerb')

create_mesh_from_regions(bounding_polygon,
    boundary_tags={'south': [0], 'east': [1], 'north': [2], 'west': [3]},
    maximum_triangle_area=0.1,
    interior_regions=interior_regions,
    filename=meshname,
    use_cache=False,
    verbose=True)

#------------------------------------------------------------------------------
# SETUP COMPUTATIONAL DOMAIN
#------------------------------------------------------------------------------

domain = anuga.Domain(meshname, use_cache=False, verbose=True)
domain.set_minimum_storable_height(0.025)
domain.set_name(outname) 

print (domain.statistics())

#------------------------------------------------------------------------------
# APPLY MANNING'S ROUGHNESSES
#------------------------------------------------------------------------------

domain.set_quantity('friction', 0.03)

# Set a Initial Water Level over the Domain
domain.set_quantity('stage', 0)

domain.set_quantity('elevation', filename=basename+'.csv', use_cache=False, verbose=True, alpha=0.99)

#------------------------------------------------------------------------------
# SETUP BOUNDARY CONDITIONS
#------------------------------------------------------------------------------

print ('Available boundary tags', domain.get_boundary_tags())

Br = anuga.Reflective_boundary(domain)  
Bd = anuga.Dirichlet_boundary([0,0,0])

domain.set_boundary({'interior': Br, 'exterior': Br, 'west': Br, 'south': Br, 'north': Br, 'east': Br})
 
# ------------------------------------------------------------------------------
# Setup inject water
# ------------------------------------------------------------------------------
input_rate = 0.102
input1_anuga_region = Region(domain, radius=1.0, center=(305694.91,6188013.94))
input1_anuga_inlet_op = Inlet_operator(domain, input1_anuga_region, Q=input_rate) # i made flow exactly the same as in DRAINS example

# ------------------------------------------------------------------------------
# Setup pipedream inlets
# ------------------------------------------------------------------------------
radius=.25

inlet1_anuga_region = Region(domain, radius=radius, center=(305698.51,6188004.63))
inlet2_anuga_region = Region(domain, radius=radius, center=(305703.39,6187999.00))
inlet3_anuga_region = Region(domain, radius=radius, center=(305713.18,6188002.02))
inlet4_anuga_region = Region(domain, radius=radius, center=(305727.24,6188004.61))
outlet_anuga_region = Region(domain, radius=radius, center=(305736.68,6188026.65))

inlet1_anuga_inlet_op = Inlet_operator(domain, inlet1_anuga_region, Q=0.0, zero_velocity=True)
inlet2_anuga_inlet_op = Inlet_operator(domain, inlet2_anuga_region, Q=0.0, zero_velocity=True)
inlet3_anuga_inlet_op = Inlet_operator(domain, inlet3_anuga_region, Q=0.0, zero_velocity=True)
inlet4_anuga_inlet_op = Inlet_operator(domain, inlet4_anuga_region, Q=0.0, zero_velocity=True)
outlet_anuga_inlet_op = Inlet_operator(domain, outlet_anuga_region, Q=0.0, zero_velocity=False)



x = domain.centroid_coordinates[:, 0]
y = domain.centroid_coordinates[:, 1]
indices = np.where(x < 10)

from pipedream_solver.hydraulics import SuperLink
import matplotlib.pyplot as plt
import pandas as pd

superjunctions = pd.DataFrame({'name' : [0, 1, 2, 3, 4], 'id' : [0, 1, 2, 3, 4], 'z_inv' : [37.5, 36.4, 34.5, 33.4, 31.0], 'h_0' : 5*[1e-5], 'bc' : 5*[False], 'storage' : 5*['functional'], 'a' : 5*[0.], 'b' : 5*[0.], 'c' : 5*[1.], 'max_depth' : 5*[np.inf], 'map_x' : 5*[0], 'map_y' : 5*[0]})

superlinks = pd.DataFrame({'name' : [0, 1, 2, 3], 'id' : [0, 1, 2, 3], 'sj_0' : [0, 1, 2, 3], 'sj_1' : [1, 2, 3, 4], 'in_offset' : 4*[0.], 'out_offset' : 4*[0.], 'dx' : [7.443, 10.251, 14.295, 24.0], 'n' : 4*[0.013], 'shape' : 4*['circular'], 'g1' : [0.375, 0.375, 0.375, 0.45], 'g2' : 4*[0.], 'g3' : 4*[0.], 'g4' : 4*[0.], 'Q_0' : 4*[0.], 'h_0' : 4*[1e-5], 'ctrl' : 4*[False], 'A_s' : 4*[0.25], 'A_c' : 4*[0.], 'C' : 4*[0.] })

superlink = SuperLink(superlinks, superjunctions, internal_links=10)

surface_elev = np.array([38.529, 37.432, 35.531, 34.393, 31.997]) # surface elevation from terrain dem

dt = 0.5    # yield step
ft = 100  # final timestep

H_js = []
losses = []

Q_iks =[]
Q_uks =[]
Q_dks =[]

for t in domain.evolve(yieldstep=dt, finaltime=ft):
    print('\n')
    domain.print_timestepping_statistics()


    # Diagnostics 
    # Compute volumes
    link_volume = ((superlink._A_ik * superlink._dx_ik).sum() +
                   (superlink._A_SIk * superlink._h_Ik).sum())
    node_volume = (superlink._A_sj * (superlink.H_j - superlink._z_inv_j)).sum()
    sewer_volume = link_volume + node_volume
    total_volume_correct = t * input_rate
    total_volume_real = domain.get_water_volume() + sewer_volume
    loss = total_volume_real - total_volume_correct

    # Append data
    losses.append(loss)
    H_js.append(superlink.H_j.copy())
    
    # record flow time series in each pipe
    Q_iks.append(superlink.Q_ik.copy())
    Q_uks.append(superlink.Q_uk.copy())
    Q_dks.append(superlink.Q_dk.copy())

    # Setup the Coupling
    anuga_depths = np.array([inlet1_anuga_inlet_op.inlet.get_average_depth(),
                             inlet2_anuga_inlet_op.inlet.get_average_depth(),
                             inlet3_anuga_inlet_op.inlet.get_average_depth(),
                             inlet4_anuga_inlet_op.inlet.get_average_depth(),
                             outlet_anuga_inlet_op.inlet.get_average_depth()])

    # Compute inflow/outflow to sewer
    C_w = 0.67
    L_w = 0.5**2 * np.pi
    Q_in = np.where(superlink.H_j <= surface_elev,
                    C_w * L_w * np.sqrt(anuga_depths) * anuga_depths,
                    C_w * L_w * np.sqrt(np.abs(anuga_depths - (superlink.H_j - surface_elev)))
                    * (anuga_depths - (superlink.H_j - surface_elev)) )

    # Simulate sewer with flow input
    superlink.step(Q_in=Q_in, dt=dt)
    superlink.reposition_junctions()

    # Add/remove flows from surface domain
    inlet1_anuga_inlet_op.set_Q(-Q_in[0])
    inlet2_anuga_inlet_op.set_Q(-Q_in[1])
    inlet3_anuga_inlet_op.set_Q(-Q_in[2])
    inlet4_anuga_inlet_op.set_Q(-Q_in[3])
    outlet_anuga_inlet_op.set_Q(-Q_in[4])


H_j = np.vstack(H_js)

plt.plot(H_j[:,0], label='Inlet 1')
plt.plot(H_j[:,1], label='Inlet 2')
plt.plot(H_j[:,2], label='Inlet 3')
plt.plot(H_j[:,3], label='Inlet 4')
plt.plot(H_j[:,4], label='Outlet')
plt.legend()
plt.title('Head at junctions')
plt.xlabel('Time (s)')
plt.ylabel('Head (m)')
plt.show()

plt.plot(losses)
plt.title('losses')
plt.show()

plt.plot(Q_uks)
plt.title('Q_uks')
plt.show()

plt.plot(Q_iks)
plt.title('Q_iks')
plt.show()

plt.plot(Q_dks)
plt.title('Q_dks')
plt.show()
