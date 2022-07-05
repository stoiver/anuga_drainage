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

basename = 'model/terrain'
outname = 'pipes'
meshname = 'model/terrain.tsh'

#------------------------------------------------------------------------------
# CREATING MESH
#------------------------------------------------------------------------------
riverWall_csv_files = glob.glob('model/wall/*.csv') # Make a list of the csv files in BREAKLINES
(riverWalls, riverWall_parameters) = su.readListOfRiverWalls(riverWall_csv_files)

CatchmentDictionary = {'model/kerb/kerb1.csv':0.01, 'model/kerb/kerb2.csv':0.01}
    
bounding_polygon = anuga.read_polygon('model/domain.csv')
interior_regions = anuga.read_polygon_dir(CatchmentDictionary, 'model/kerb')

import numpy as num
b_polygon = num.array(bounding_polygon)
print(b_polygon)
import matplotlib.pyplot as plt

plt.figure()
plt.plot(b_polygon[:,0],b_polygon[:,1])
for i, v in enumerate(b_polygon):
    plt.annotate(str(v), xy=v, xytext=(-7,7), textcoords='offset points')
plt.pause(0.01)


create_mesh_from_regions(bounding_polygon,
    boundary_tags={'inflow': [12], 'bottom': [0,1,2,3,4,5], 'top': [7,8,9,10,11], 'outflow': [6]},
    #boundary_tags=None,
    maximum_triangle_area=0.1,
    breaklines=riverWalls.values(),
    interior_regions=interior_regions,
    filename=meshname,
    use_cache=False,
    verbose=True)

#------------------------------------------------------------------------------
# SETUP COMPUTATIONAL DOMAIN
#------------------------------------------------------------------------------

domain = anuga.Domain(meshname, use_cache=False, verbose=True)
domain.set_minimum_storable_height(0.0)
domain.riverwallData.create_riverwalls(riverWalls) 
domain.set_name(outname) 

print (domain.statistics())

#------------------------------------------------------------------------------
# APPLY MANNING'S ROUGHNESSES
#------------------------------------------------------------------------------

domain.set_quantity('friction', 0.025)

# Set a Initial Water Level over the Domain
domain.set_quantity('stage', 0)

domain.set_quantity('elevation', filename=basename+'.csv', use_cache=False, verbose=True, alpha=0.99)

#------------------------------------------------------------------------------
# SETUP BOUNDARY CONDITIONS
#------------------------------------------------------------------------------

print ('Available boundary tags', domain.get_boundary_tags())

Br = anuga.Reflective_boundary(domain)  
Bd = anuga.Dirichlet_boundary([0,0,0])
#Bt = anuga.Flather_external_stage_zero_velocity_boundary()

domain.set_boundary({'inflow': Br, 'bottom': Br, 'outflow': Bd, 'top': Br})
#domain.set_boundary({'exterior' : Bd})
 
# ------------------------------------------------------------------------------
# Setup inject water
# ------------------------------------------------------------------------------
input_rate = 0.05 #  0.102 # i made inflow exactly the same as in DRAINS example
input1_anuga_region = Region(domain, radius=1.0, center=(305694.91,6188013.94))
input1_anuga_inlet_op = Inlet_operator(domain, input1_anuga_region, Q=input_rate) 

# ------------------------------------------------------------------------------
# Setup pipedream inlets
# ------------------------------------------------------------------------------
radius=0.25

inlet1_anuga_region = Region(domain, radius=radius, center=(305698.51,6188004.63))
inlet2_anuga_region = Region(domain, radius=radius, center=(305703.39,6187999.00))
inlet3_anuga_region = Region(domain, radius=radius, center=(305713.18,6188002.02))
inlet4_anuga_region = Region(domain, radius=radius, center=(305727.24,6188004.61))
outlet_anuga_region = Region(domain, radius=radius, center=(305736.68,6188026.65))

inlet1_anuga_inlet_op = Inlet_operator(domain, inlet1_anuga_region, Q=0.0, zero_velocity=True)
inlet2_anuga_inlet_op = Inlet_operator(domain, inlet2_anuga_region, Q=0.0, zero_velocity=True)
inlet3_anuga_inlet_op = Inlet_operator(domain, inlet3_anuga_region, Q=0.0, zero_velocity=True)
inlet4_anuga_inlet_op = Inlet_operator(domain, inlet4_anuga_region, Q=0.0, zero_velocity=True)
inlet4_anuga_inlet_op = Inlet_operator(domain, outlet_anuga_region, Q=0.0, zero_velocity=True)

anuga_elevs = np.array([inlet1_anuga_inlet_op.inlet.get_average_elevation(),
                        inlet2_anuga_inlet_op.inlet.get_average_elevation(),
                        inlet3_anuga_inlet_op.inlet.get_average_elevation(),
                        inlet4_anuga_inlet_op.inlet.get_average_elevation(),
                        inlet4_anuga_inlet_op.inlet.get_average_elevation()])

x = domain.centroid_coordinates[:, 0]
y = domain.centroid_coordinates[:, 1]
indices = np.where(x < 10)

from pipedream_solver.hydraulics import SuperLink
import matplotlib.pyplot as plt
import pandas as pd

superjunctions = pd.DataFrame({'name': [0, 1, 2, 3, 4],
                               'id': [0, 1, 2, 3, 4],
                               'z_inv': [37.5, 36.4, 34.5, 32.0, 32.0],
                               'h_0': 5*[1e-5],
                               'bc': 5*[False],
                               'storage': 5*['functional'],
                               'a': 5*[0.],
                               'b': 5*[0.],
                               'c': 5*[1.],
                               'max_depth': 5*[np.inf],
                               'map_x': 5*[0],
                               'map_y': 5*[0]})

superlinks = pd.DataFrame({'name': [0, 1, 2, 3],
                           'id': [0, 1, 2, 3],
                           'sj_0': [0, 1, 2, 3],
                           'sj_1': [1, 2, 3, 4],
                           'in_offset': 4*[0.],
                           'out_offset': 4*[0.],
                           'dx': [7.4, 10.3, 14.3, 24.0],
                           'n': 4*[0.013],
                           'shape': 4*['circular'],
                           'g1': [0.375, 0.375, 0.375, 0.45],
                           'g2': 4*[0.],
                           'g3': 4*[0.],
                           'g4': 4*[0.],
                           'Q_0': 4*[0.],
                           'h_0': 4*[1e-5],
                           'ctrl': 4*[False],
                           'A_s': 4*[1.12],
                           'A_c': 4*[0.],
                           'C': 4*[0.]})  # A_s surface area of pit (1sqm) + 1.2m lintel (0.1x1.2m long = 0.12sqm)

superlink = SuperLink(superlinks, superjunctions, internal_links=20)

surface_elev = np.array([38.5489161, 37.46895542, 35.59301592, 32.00173483, 32.00173483]) # surface elevation from terrain dem

dt = 1.0     # yield step
ft = 250.0   # final timestep

H_js = []
losses = []

Q_iks =[]
Q_uks =[]
Q_dks =[]
times =[]

Q_in_old = np.zeros_like(anuga_elevs)
time_average = 10 # sec

for t in domain.evolve(yieldstep=dt, finaltime=ft):
    #print('\n')
    if domain.yieldstep_counter%domain.output_frequency == 0:
        domain.print_timestepping_statistics()

    #print(domain.volumetric_balance_statistics())

    anuga_depths = np.array([inlet1_anuga_inlet_op.inlet.get_average_depth(),
                             inlet2_anuga_inlet_op.inlet.get_average_depth(),
                             inlet3_anuga_inlet_op.inlet.get_average_depth(),
                             inlet4_anuga_inlet_op.inlet.get_average_depth(),
                             inlet4_anuga_inlet_op.inlet.get_average_depth()])

    anuga_stages = np.array([inlet1_anuga_inlet_op.inlet.get_average_stage(),
                             inlet2_anuga_inlet_op.inlet.get_average_stage(),
                             inlet3_anuga_inlet_op.inlet.get_average_stage(),
                             inlet4_anuga_inlet_op.inlet.get_average_stage(),
                             inlet4_anuga_inlet_op.inlet.get_average_stage()])



    # Compute inflow/outflow to sewer
    C_w = 0.67
    L_w = 0.25**2 * np.pi

    Q_in = np.where(superlink.H_j <= anuga_elevs,
                    C_w * L_w * np.sqrt(anuga_depths) * anuga_depths,
                    C_w * L_w * np.sqrt(superlink.H_j - anuga_elevs)
                    * (anuga_elevs- superlink.H_j ))


    # average it out
    Q_in = ((time_average - dt)*Q_in_old + dt*Q_in)/time_average
    Q_in_old = Q_in

    #Q_in = np.where(superlink.H_j <= anuga_stages,
    #               C_w * L_w * np.sqrt(anuga_depths) * anuga_depths, 0.0 )


    if domain.yieldstep_counter%domain.output_frequency == 0:    
        print(anuga_depths)
        print(anuga_elevs)
        print(anuga_stages)
        
        print(superlink.H_j)
        print(Q_in)

 
    # Simulate sewer with flow input
    superlink.step(Q_in=Q_in, dt=dt)
    superlink.reposition_junctions()

    # Add/remove flows from surface domain
    inlet1_anuga_inlet_op.set_Q(-Q_in[0])
    inlet2_anuga_inlet_op.set_Q(-Q_in[1])
    inlet3_anuga_inlet_op.set_Q(-Q_in[2])
    inlet4_anuga_inlet_op.set_Q(-Q_in[3])
    inlet4_anuga_inlet_op.set_Q(-Q_in[4])

    # Compute volumes
    link_volume = ((superlink._A_ik * superlink._dx_ik).sum() +
                   (superlink._A_SIk * superlink._h_Ik).sum())
    node_volume = (superlink._A_sj * (superlink.H_j - superlink._z_inv_j)).sum()
    sewer_volume = link_volume + node_volume

    boundary_flow = domain.get_boundary_flux_integral()
    total_volume_correct = t*input_rate + boundary_flow

    total_volume_real = domain.get_water_volume() + sewer_volume
    loss = total_volume_real - total_volume_correct
    

    if domain.yieldstep_counter%domain.output_frequency == 0:
        print('boundary flow', boundary_flow)

    # Append data
    times.append(t)
    losses.append(loss)
    H_js.append(superlink.H_j.copy())
    
    # record flow time series in each pipe
    Q_iks.append(superlink.Q_ik.copy())
    Q_uks.append(superlink.Q_uk.copy())
    Q_dks.append(superlink.Q_dk.copy())

    H_j = np.vstack(H_js)

    if domain.yieldstep_counter%domain.output_frequency == 0:
        plt.clf()
        plt.plot(times,H_j[:,0], label='Inlet 1')
        plt.plot(times,H_j[:,1], label='Inlet 2')
        plt.plot(times,H_j[:,2], label='Inlet 3')
        plt.plot(times,H_j[:,3], label='Inlet 4')
        plt.plot(times,H_j[:,4], label='Outlet')
        plt.legend()
        plt.title('Head at junctions')
        plt.xlabel('Time (s)')
        plt.ylabel('Head (m)')
        plt.pause(0.01)



H_j = np.vstack(H_js)


plt.plot(times,H_j[:,0], label='Inlet 1')
plt.plot(times,H_j[:,1], label='Inlet 2')
plt.plot(times,H_j[:,2], label='Inlet 3')
plt.plot(times,H_j[:,3], label='Inlet 4')
plt.plot(times,H_j[:,4], label='Outlet')
plt.legend()
plt.title('Head at junctions')
plt.xlabel('Time (s)')
plt.ylabel('Head (m)')
plt.savefig('Figure_01.png')
plt.show()

plt.plot(times,losses)
plt.title('losses (total_volume_real - total_volume_correct)')
plt.savefig('Figure_02.png')
plt.show()

plt.plot(times,Q_iks)
plt.title('Link flows (m^3/s)')
plt.savefig('Figure_03.png')
plt.show()

plt.plot(times,Q_uks)
plt.title('Flows into upstream ends of superlinks (m^3/s) ')
plt.savefig('Figure_04.png')
plt.show()

plt.plot(times,Q_dks)
plt.title('Flows into downstream ends of superlinks (m^3/s) ')
plt.savefig('Figure_05.png')
plt.show()
