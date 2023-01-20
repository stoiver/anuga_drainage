"""
    This test is trying to build a Y-shape pipe system in anuga-swmm
    It has two inlet nodes in the same altitude and connect to one outlet node which is in lower altitude
    and see how water flows in the simple system.
    Using the same domain structure like previous testing "test_volume_inlet_op.py"
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
import math

# ------------------------------------------------------------------------------
# Setup computational domain
# ------------------------------------------------------------------------------

length = 20.
width = 6.
dx = dy = 0.2  # .1           # Resolution: Length of subdivisions on both axes

domain = rectangular_cross_domain(int(length / dx), int(width / dy),
                                               len1=length, len2=width)
domain.set_name('Y_shape_2inlets')  # Output name based on script name. You can add timestamp=True


# print(domain.statistics())


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
    z[id] = 3

    # wall, spilt two inlet areas
    id = (10 > x) & (2.5 < y) & (y < 3.5)
    z[id] = 3

    '''
    not consider the pit's physical depth, just note the physical location
    to setup the depth, just need to modify z[id] = <x>, x represent x meters depth
    '''

    # first pit, locate at (7, 1) radius=1
    id = (x - 7) ** 2 + (y - 1) ** 2 < 1. ** 2
    z[id] -= 0.0

    # second pit, locate at (7, 5) radius=1
    id = (x - 7) ** 2 + (y - 5) ** 2 < 1. ** 2
    z[id] -= 0.0

    # out pit, locate at (17, 3) radius=1
    id = (x - 17) ** 2 + (y - 3) ** 2 < 1 ** 2
    z[id] -= 0.0

    return z


# ------------------------------------------------------------------------------
# Setup initial quantity
# ------------------------------------------------------------------------------

domain.set_quantity('elevation', topography, location='centroids')  # elevation is a function
domain.set_quantity('friction', 0.01)  # Constant friction
domain.set_quantity('stage', expression='elevation', location='centroids')  # Dry initial condition

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

inlet1_anuga_region = Region(domain, radius=1.0, center=(7., 1.))
inlet2_anuga_region = Region(domain, radius=1.0, center=(7., 5.))
outlet_anuga_region = Region(domain, radius=1.0, center=(17., 3.))

input1_anuga_region = Region(domain, radius=1.0, center=(2., 1.))
# region_input2 = Region(domain, radius=1.0, center=(2., 5.))

inlet1_anuga_inlet_op = Inlet_operator(domain, inlet1_anuga_region, Q=0.0, zero_velocity=True)
inlet2_anuga_inlet_op = Inlet_operator(domain, inlet2_anuga_region, Q=0.0, zero_velocity=True)
outlet_anuga_inlet_op = Inlet_operator(domain, outlet_anuga_region, Q=0.0, zero_velocity=False)

input1_anuga_inlet_op = Inlet_operator(domain, input1_anuga_region, Q=1) #only inlet1 has water inflow.
# op_input2 = Inlet_operator(domain, region_input2, Q=0.1)

x = domain.centroid_coordinates[:, 0]
y = domain.centroid_coordinates[:, 1]
indices = num.where(x < 10)
# anuga.Set_stage(domain, stage=-2.75, indices=indices)()

# ------------------------------------------------------------------------------
# couple SWMM
# ------------------------------------------------------------------------------
# TODO: couple SWMM
from pyswmm import Simulation, Nodes, Links

sim = Simulation('./2inlets_short.inp')
sim.start()
print('\nsim start?? ', sim._isStarted)
node_names = ['J1', 'J2', 'Out1']
link_names = ['C1', 'C2']

nodes = [Nodes(sim)[names] for names in node_names]
links = [Links(sim)[names] for names in link_names]

culvert1 = Links(sim)['C1']
culvert2 = Links(sim)['C2']

inlet1_pyswmm_node = Nodes(sim)['J1']
inlet2_pyswmm_node = Nodes(sim)['J2']
outlet_pyswmm_node = Nodes(sim)['Out1']

# type, area, length, orifice_coeff, free_weir_coeff, submerged_weir_coeff

#initialize parameters of nodes and links
for i in range(len(node_names)):
    opening = nodes[i].create_opening(4, math.pi, 1.0, 0.6, 1.6, 1.0)
    nodes[i].coupling_area = math.pi
    # print('node opening? ', node_names[i], ' ', opening)
    print('node coupled? ', node_names[i], ' ', nodes[i].is_coupled)

# essential parameters
initial_volume = domain.get_water_volume()
print("initial volume: %.4f:" % initial_volume)
input_velocity = 1
total_input_flow = 0
total_output_flow = 0
loss = 0

# use to plot
losses = []
residual_losses = []
Q_inlet1_volumes = []
Q_inlet2_volumes = []
Q_inlet1_depth = []
Q_inlet2_depth = []
Q_acc_volumes = []
flood = 0

print('initial culvert volume: %.4f' % (culvert1.volume+culvert2.volume))

for t in domain.evolve(yieldstep=1.0, finaltime=100.0):
    print('\n')
    domain.print_timestepping_statistics()


    # Test real total volume by Chen
    total_volume_correct = t * input_velocity
    total_volume_real = domain.get_water_volume()+culvert1.volume+culvert2.volume
    loss = total_volume_real - total_volume_correct

    print('loss %.4f' % loss)
    print("correct water ", total_volume_correct)
    print("water on domain", domain.get_water_volume())
    print("expect water in pipe ", total_volume_correct-domain.get_water_volume())
    print("water in pipe: ", culvert1.volume + culvert2.volume)

    #overland depth will be calculated by the average value of current depth.
    inlet1_pyswmm_node.overland_depth = inlet1_anuga_inlet_op.inlet.get_average_depth()
    inlet2_pyswmm_node.overland_depth = inlet2_anuga_inlet_op.inlet.get_average_depth()

    #volumes parameter refers to a list of tuple of dictionary, formatted by k-v: nodes, volumes during simulation.
    volumes = sim.coupling_step(1.0)
    volumes_in_out = volumes[-1][-1]

    # FIXME: determine Q_inlet1, Q_outlet
    Q_inlet1 = -volumes[-1][1]['J1']
    Q_inlet2 = -volumes[-1][1]['J2']
    Q_outlet = outlet_pyswmm_node.total_inflow
    total_input_flow += Q_inlet1 + Q_inlet2 + Q_outlet

    print("Q_inlet1/inlet1 flow: %.4f/%.4f" % (Q_inlet1, inlet1_pyswmm_node.coupling_inflow))
    print("Q_inlet2/conduit2 flow:  %.4f/%.4f" % (Q_inlet2, culvert2.flow))
    print("Q_outlet: %.4f" % Q_outlet)
    print("inlet 1 depth/ inlet 2 depth/ outlet depth (%.4f/%.4f/%.4f)" % (inlet1_pyswmm_node.depth, inlet2_pyswmm_node.depth, outlet_pyswmm_node.depth))
    print("expected/cal/real in pipe: %.4f/%.4f/%.4f" % (total_volume_correct - domain.get_water_volume(), total_input_flow,culvert1.volume+culvert2.volume))

    max_water = 10*math.pi + math.pi*(0.05**2)*11

    print("maximum water in pipe: ", max_water)
    print(inlet1_pyswmm_node.statistics['flooding_volume']* 0.0283168466)
    print(inlet1_pyswmm_node.flooding * 0.0283168466)
    print(inlet2_pyswmm_node.flooding * 0.0283168466)

    flood += (inlet1_pyswmm_node.flooding+inlet2_pyswmm_node.flooding)*0.0283168466 #coefficiency between cubid meter and gallon.

    print("cumulative flood ", flood)

    # if inlet1_pyswmm_node.flooding > 0 and Q_inlet1 < 0:
    #     Q_inlet1 = 0
    # if inlet2_pyswmm_node.flooding > 0 and Q_inlet2 < 0:
    #     Q_inlet2 = 0

    #Q: recorded as the raw data of next stage to be calculated as average depth.
    inlet1_anuga_inlet_op.set_Q(Q_inlet1+inlet1_pyswmm_node.flooding*0.0283168466)
    inlet2_anuga_inlet_op.set_Q(Q_inlet2+inlet2_pyswmm_node.flooding*0.0283168466)
    outlet_anuga_inlet_op.set_Q(Q_outlet)

    losses.append(loss)





print('visualize the loss curve')

t = int(t)
import matplotlib.pyplot as plt
plt.subplot(2, 2, 1)
plt.title('loss')
plt.plot(range(t+1), losses)
plt.show()
