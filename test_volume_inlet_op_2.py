"""
This is the first prototype (testing) for the ANUGA & SWMM coupling project

In this testing, we are expecting to create a one-pipe testing. Flowing water out from ANUGA to SWMM, and using the SWMM
calculate the water flow activities in the pipe, and flows back to ANUGA.

we can validate this testing by monitor the change of the water total volume. It should remains the same between flowing
to SWMM and flowing back to ANUGA.
"""


# ------------------------------------------------------------------------------
# Import necessary modules
# ------------------------------------------------------------------------------
from anuga import Dirichlet_boundary
from anuga import Domain
from anuga import Reflective_boundary
from anuga import Rate_operator
from anuga import Inlet_operator
from anuga import Region
from anuga import rectangular_cross_domain

import anuga
import numpy as num

# ------------------------------------------------------------------------------
# Setup computational domain
# ------------------------------------------------------------------------------

length = 20.
width = 4.
dx = dy = 0.2 # .1           # Resolution: Length of subdivisions on both axes

domain = rectangular_cross_domain(int(length / dx), int(width / dy),
                                               len1=length, len2=width)
domain.set_name('total_volume_testing')  # Output name based on script name. You can add timestamp=True
print(domain.statistics())


# ------------------------------------------------------------------------------
# Setup initial conditions
# ------------------------------------------------------------------------------
def topography(x, y):
    """Complex topography defined by a function of vectors x and y."""

    z = 0 * x - 5

    # higher pools
    id = x < 10
    z[id] = -3

    # wall
    id = (10 < x) & (x < 15)
    z[id] = 0

    # inflow pipe hole, located at (2, 2), r = 0.5, depth 0.1
    id = (x - 7) ** 2 + (y - 2) ** 2 < 0.3 ** 2
    z[id] -= 0.0

    # inflow pipe hole, located at (12, 2), r = 0.5, depth 0.1
    id = (x - 17) ** 2 + (y - 2) ** 2 < 0.3 ** 2
    z[id] -= 0.0

    return z


# ------------------------------------------------------------------------------
# Setup initial quantity
# ------------------------------------------------------------------------------
domain.set_quantity('elevation', topography, location = 'centroids')  # elevation is a function
domain.set_quantity('friction', 0.01)  # Constant friction
domain.set_quantity('stage', expression='elevation', location = 'centroids')  # Dry initial condition
# --------------------------

"""
We would use this method to gain the boundary indices
"""


# polygon1 = [ [10.0, 0.0], [11.0, 0.0], [11.0, 5.0], [10.0, 5.0] ]
# polygon2 = [ [10.0, 0.2], [11.0, 0.2], [11.0, 4.8], [10.0, 4.8] ]





# ------------------------------------------------------------------------------
# Setup boundaries
# ------------------------------------------------------------------------------
Bi = Dirichlet_boundary([-2.75, 0, 0])  # Inflow
Br = Reflective_boundary(domain)  # Solid reflective wall
Bo = Dirichlet_boundary([-5, 0, 0])  # Outflow

domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

# ------------------------------------------------------------------------------
# Setup inject water
# ------------------------------------------------------------------------------

region_inlet  = Region(domain, radius=1.0, center=(7., 2.))
region_outlet = Region(domain, radius=1.0, center=(17., 2.))

region_input = Region(domain, radius=1.0, center=(2., 2.))

op_input = Inlet_operator(domain, region_input, Q=0.25)

op_inlet = Inlet_operator(domain, region_inlet, Q=0.0, zero_velocity=True)
op_outlet = Inlet_operator(domain, region_outlet, Q=0.0, zero_velocity=False)  #

x = domain.centroid_coordinates[:,0]
indices = num.where(x < 10)
anuga.Set_stage(domain, stage = -2.75, indices = indices)()

from pyswmm import Simulation, Nodes, Links

sim = Simulation('./pipe_test.inp')
sim.start()
node_names = ['Inlet', 'Outlet']

link_names = ['Culvert']

nodes = [Nodes(sim)[names] for names in node_names]
links = [Links(sim)[names] for names in link_names]


culvert = Links(sim)['Culvert']
inlet   = Nodes(sim)['Inlet']
outlet  = Nodes(sim)['Outlet']

print('')

# type, area, length, orifice_coeff, free_weir_coeff, submerged_weir_coeff
inlet_opening = inlet.create_opening(4, 2.0, 1.0, 0.6, 1.6, 1.0)
inlet.coupling_area = 0.5
print('inlet_opening',inlet_opening)

# TODO: setup the outlet node
#outlet.create_opening(4, 1.0, 1.0, 0.6, 1.6, 1.0)
#outlet.coupling_area = 1.0


print("inlet_is_open?:", inlet.is_coupled)
print("outlet_is_open?:", outlet.is_coupled)

flow = 1.0
stop_release_water_time = 0 # the time for stopping releasing the water
initial_volume = domain.get_water_volume()
inlet_flow = 0.0
outlet_flow = 0.0
inlet_residual = 0.0
inlet_discrepency = 0.0

previous_culvert_volume = 0.0


domain.set_name("anuga_swmm")

for t in domain.evolve(yieldstep=1.0, finaltime=60.0):
    print("\n")
    #print(f"coupling step: {t}")
    print(80*'=')
    domain.print_timestepping_statistics()
    print(80*'=')

    
    # FIXME SR: Don't know why including inlet.depth seems to reconcile the volumes
    # Note that including inlet.coupling_area doesnt seem to work
    total_volume = domain.get_water_volume() + culvert.volume - inlet.depth
    print("total volume: ", total_volume)
    print("correct volume: ",initial_volume + t*0.25)

    print("inlet.depth",inlet.depth)
    print("domain volume: ",domain.get_water_volume())
    print('culvert volume', culvert.volume)
    
    print('Discrepancy', initial_volume + t*0.25 - total_volume)

    print('cum inlet_flow',inlet_flow)
    print('cum outlet_flow',outlet_flow)
    print('cum pipe flow', inlet_flow + outlet_flow)
    #print("inlet_discrepency",inlet_discrepency)



    previous_culvert_volume = culvert.volume

    #print("culvert depth",culvert.depth)




    #print('Outlet outflow', outlet.total_outflow)

    #print('Inlet coupling inflow', inlet.coupling_inflow)


    print(24*'-')
    print("Flows for next timestep")
    print(24*'-')

    # Setup inlet for SWMM step
    inlet.overland_depth = op_inlet.inlet.get_average_depth()
    print("inlet overland depth: ", inlet.overland_depth)
    

    volumes = sim.coupling_step(1.0)
    volumes_in_out = volumes[-1][-1]

    print(volumes)


    #print('Inlet volumes', op_inlet.domain.fractional_step_volume_integral)
    #print('Outlet volume', op_outlet.domain.fractional_step_volume_integral)

    




    #print ("node area", inlet.coupling_area)
    #print ("node inflow", inlet.coupling_inflow)
    #print ("node volume", inlet.volume)
    #print("node stuff", nodes[0].statistics)

    #print(volumes)

    print('Culvert volume', culvert.volume)
    #print('Change in culvert volume',culvert.volume-previous_culvert_volume)
    print('Culvert Flow', culvert.flow)
    print('Inlet inflow', inlet.total_inflow)
    #print('Inlet ourflow', inlet.total_outflow)
    #inlet_residual = inlet.total_inflow - inlet.total_outflow
    #print('inlet_residual',inlet_residual)
    print('Outlet inflow', outlet.total_inflow)
    
    #print(volumes_in_out)

    #print("Volume total at node Inlet" ":", volumes_in_out["Inlet"])
    #print("Volume total at node Outlet" ":", volumes_in_out["Outlet"])


    #inlet_Q  = -volumes_in_out['Inlet']
    #inlet_Q = - inlet.coupling_inflow
    inlet_Q = - inlet.total_inflow
    outlet_Q = outlet.total_inflow

    op_inlet.set_Q(inlet_Q)
    op_outlet.set_Q(outlet_Q)

    inlet_flow = inlet_flow + inlet_Q
    outlet_flow = outlet_flow + outlet_Q

    #inlet_discrepency = inlet.total_inflow + inlet_Q

    #print("inlet_discrepency",inlet_discrepency)

    previous_inlet_flow = inlet.total_inflow #volumes_in_out['Inlet']
    previous_outlet_flow = outlet.total_inflow

