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
from anuga import rectangular_cross

import anuga
import numpy as num

# ------------------------------------------------------------------------------
# Setup computational domain
# ------------------------------------------------------------------------------

length = 20.
width = 4.
dx = dy = 0.2 # .1           # Resolution: Length of subdivisions on both axes

points, vertices, boundary = rectangular_cross(int(length / dx), int(width / dy),
                                               len1=length, len2=width)
domain = Domain(points, vertices, boundary)
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


Culvert = Links(sim)['Culvert']

# type, area, length, orifice_coeff, free_weir_coeff, submerged_weir_coeff
nodes[0].create_opening(4, 1.0, 1.0, 0.6, 1.6, 1.0)
nodes[0].coupling_area = 0.25

# TODO: setup the outlet node
nodes[1].create_opening(4, 1.0, 1.0, 0.6, 1.6, 1.0)
nodes[1].coupling_area = 0.25

print('')
print("node0_is_open?:", nodes[0].is_coupled)
print("node1_is_open?:", nodes[1].is_coupled)

flow = 1.0
stop_release_water_time = 0 # the time for stopping releasing the water
initial_volume = domain.get_water_volume()
previous_inlet_flow = 0.0
previous_outlet_flow = 0.0


domain.set_name("anuga_swmm")

for t in domain.evolve(yieldstep=1.0, finaltime=60.0):
    print("\n")
    #print(f"coupling step: {t}")
    domain.print_timestepping_statistics()

    # set the overland_depth
    # TODO: set up the overland depth, modify this function


    volumes = sim.coupling_step(1.0)
    volumes_in_out = volumes[-1][-1]

    # FIXME SR: Shouldn't this go before calculating the coupling step.
    # But that messes up the volume balance!
    nodes[0].overland_depth = op_inlet.inlet.get_average_depth()
    
    #print('Inlet volumes', op_inlet.domain.fractional_step_volume_integral)
    #print('Outlet volume', op_outlet.domain.fractional_step_volume_integral)

    total_volume = domain.get_water_volume() + Culvert.volume + nodes[0].overland_depth

    print("total volume: ", total_volume)
    print("correct volume: ",initial_volume + t*0.25)
    #print("calc volume: ",initial_volume + op_inlet.domain.fractional_step_volume_integral + Culvert.volume)
    #print('Loss', t*0.25 - op_inlet.domain.fractional_step_volume_integral - nodes[0].overland_depth - Culvert.volume )
    print('Loss', initial_volume + t*0.25 - total_volume)

    print('Culvert volume', Culvert.volume)

    print('Culvert Flow', Culvert.flow)

    print("inlet overland depth: ", op_inlet.inlet.get_average_depth())

    print('Inlet inflow', nodes[0].total_inflow)
    print('Inlet ourflow', nodes[0].total_inflow)
    print('Outlet inflow', nodes[1].total_inflow)
    print('Outlet ourflow', nodes[1].total_inflow)

    #print ("node area", nodes[0].coupling_area)
    #print ("node inflow", nodes[0].coupling_inflow)
    #print ("node volume", nodes[0].volume)
    #print("node stuff", nodes[0].statistics)

    #print(volumes)


    
    #print(volumes_in_out)

    print("Volume total at node Inlet" ":", volumes_in_out["Inlet"])
    print("Volume total at node Outlet" ":", volumes_in_out["Outlet"])
    print("Outlet: ", nodes[1].total_inflow)
    op_inlet.set_Q(-1 * volumes_in_out['Inlet'])

    #op_inlet.set_Q(-nodes[0].total_inflow)
    op_outlet.set_Q(nodes[1].total_inflow)



    previous_inlet_flow = nodes[0].total_inflow #volumes_in_out['Inlet']
    previous_outlet_flow = nodes[1].total_inflow

