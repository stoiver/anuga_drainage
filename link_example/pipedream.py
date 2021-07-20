#------------------------------------------------------------------------------
# IMPORT NECESSARY MODULES
#------------------------------------------------------------------------------
print (' ABOUT to Start Simulation:- Importing Modules')

import anuga, anuga.parallel, numpy, time, os, glob
from anuga.operators.rate_operators import Polygonal_rate_operator
from anuga import file_function, Polygon_function, read_polygon, create_mesh_from_regions, Domain, Inlet_operator
import anuga.utilities.spatialInputUtil as su

from anuga import distribute, myid, numprocs, finalize, barrier
from anuga.parallel.parallel_operator_factory import Inlet_operator, Boyd_box_operator, Boyd_pipe_operator
from anuga import Rate_operator
from anuga import Region

import numpy as np
import math
#------------------------------------------------------------------------------
# FILENAMES, MODEL DOMAIN and VARIABLES
#------------------------------------------------------------------------------

basename = 'terrain'
outname = 'pipedream'
meshname = 'terrain.msh'


W=296600.
N=6180070.

E=296730.
S=6179960.

#------------------------------------------------------------------------------
# CREATING MESH
#------------------------------------------------------------------------------

bounding_polygon = [[W, S], [E, S], [E, N], [W, N]]


create_mesh_from_regions(bounding_polygon,
    boundary_tags={'south': [0], 'east': [1], 'north': [2], 'west': [3]},
    maximum_triangle_area=1.0,
    filename=meshname,
    use_cache=False, 
    verbose=True)

#------------------------------------------------------------------------------
# SETUP COMPUTATIONAL DOMAIN
#------------------------------------------------------------------------------

domain = Domain(meshname, use_cache=False, verbose=True)
domain.set_minimum_storable_height(0.0001) 
domain.set_name(outname) 
print (domain.statistics())

#------------------------------------------------------------------------------
# APPLY MANNING'S ROUGHNESSES
#------------------------------------------------------------------------------

domain.set_quantity('friction', 0.035)

domain.set_quantity('elevation', filename=basename+'.csv', use_cache=False, verbose=True, alpha=0.1)

#------------------------------------------------------------------------------
# SETUP BOUNDARY CONDITIONS
#------------------------------------------------------------------------------

print ('Available boundary tags', domain.get_boundary_tags())

Br = anuga.Reflective_boundary(domain)
Bd = anuga.Dirichlet_boundary([0,0,0])

domain.set_boundary({'west': Bd, 'south': Br, 'north': Bd, 'east': Bd})

inlet1_anuga_region = Region(domain, radius=0.5, center=(296660.390,6180017.186))
outlet_anuga_region = Region(domain, radius=0.5, center=(296649.976,6180038.872))

inlet1_anuga_inlet_op = Inlet_operator(domain, inlet1_anuga_region, Q=0.0, zero_velocity=True)
outlet_anuga_inlet_op = Inlet_operator(domain, outlet_anuga_region, Q=0.0, zero_velocity=False)

line=[[296669.258,6179974.191],[296677.321,6179976.449]]
Inlet_operator(domain, line, 1.0)

#------------------------------------------------------------------------------
# PIPEDREAM
#------------------------------------------------------------------------------

from pipedream_solver.hydraulics import SuperLink
import matplotlib.pyplot as plt
import pandas as pd


superjunctions = pd.DataFrame({'name' : [0, 1], 'id' : [0, 1], 'z_inv' : [12.4, 12.2], 'h_0' : 2*[1e-5], 'bc' : 2*[False], 'storage' : 2*['functional'], 'a' : 2*[0.], 'b' : 2*[1.], 'c' : 2*[10.], 'max_depth' : 2*[np.inf], 'map_x' : 2*[0], 'map_y' : 2*[0]})

superlinks = pd.DataFrame({'name' : [0], 'id' : [0], 'sj_0' : [0], 'sj_1' : [1], 'in_offset' : 1*[0.], 'out_offset' : 1*[0.], 'dx' : [24], 'n' : 1*[0.013], 'shape' : 1*['circular'], 'g1' : [1.], 'g2' : 1*[0.], 'g3' : 1*[0.], 'g4' : 1*[0.], 'Q_0' : 1*[0.], 'h_0' : 1*[1e-5], 'ctrl' : 1*[False], 'A_s' : 1*[0.], 'A_c' : 1*[0.], 'C' : 1*[0.] })

superlink = SuperLink(superlinks, superjunctions, internal_links=20)

surface_elev = np.array([12.4, 12.2]) - 0.01

input_velocity = 1
dt = 1    # yield step
ft = 4000   # final timestep

H_js = []
losses = []

Q_iks =[]
Q_uks =[]
Q_dks =[]

for t in domain.evolve(yieldstep=dt, finaltime=ft):
    print('\n')
    domain.print_timestepping_statistics()

    anuga_depths = np.array([inlet1_anuga_inlet_op.inlet.get_average_depth(),
                             outlet_anuga_inlet_op.inlet.get_average_depth()])

    # Compute inflow/outflow to sewer
    C_o = 0.67
    A_o = 2 * np.pi
    Q_in = C_o * A_o * np.sign(anuga_depths - (superlink.H_j - superlink._z_inv_j)) * np.sqrt(np.abs(anuga_depths - (superlink.H_j - superlink._z_inv_j)))

    # Simulate sewer with flow input
    superlink.step(Q_in=Q_in, dt=dt)
    superlink.reposition_junctions()

    # Add/remove flows from surface domain
    inlet1_anuga_inlet_op.set_Q(-Q_in[0])
    outlet_anuga_inlet_op.set_Q(-Q_in[1])

    # Compute volumes
    link_volume = ((superlink._A_ik * superlink._dx_ik).sum() +
                   (superlink._A_SIk * superlink._h_Ik).sum())
    node_volume = (superlink._A_sj * (superlink.H_j - superlink._z_inv_j)).sum()
    sewer_volume = link_volume + node_volume
    total_volume_correct = t * input_velocity
    total_volume_real = domain.get_water_volume() + sewer_volume
    loss = total_volume_real - total_volume_correct

    # Append data
    losses.append(loss)
    H_js.append(superlink.H_j.copy())

    # record flow time series in each pipe
    Q_iks.append(superlink.Q_ik.copy())
    Q_uks.append(superlink.Q_ik.copy())
    Q_dks.append(superlink.Q_ik.copy())

H_j = np.vstack(H_js)

plt.plot(H_j[:,0], label='Inlet 1')
plt.plot(H_j[:,1], label='Inlet 2')
plt.legend()
plt.title('Head at junctions')
plt.xlabel('Time (s)')
plt.ylabel('Head (m)')
plt.show()

plt.clf()
plt.plot(losses)
plt.show()


plt.clf()
plt.plot(Q_dks)
plt.show()
