"""pyswmm example
This is a simple example based on previous work by:
Stephen Roberts, Petar Milevski, Rudy van Drie, Ole Nielsen
in December 2018
It was expanded by Zhixian Wu in September 2020 to demonstrate
some changes she made to the SWMM and PySWMM code.
"""

from pyswmm import Simulation, Nodes, Links


def run_swmm():
    # =======================================
    # setup all the nodes before starting
    # =======================================

    sim = Simulation('./pipe_test.inp')

    # =======================================
    # Start the simulation
    # =======================================
    print(sim._isStarted)

    sim.start()

    print(sim._isStarted)

    node_names = ['Inlet', 'Outlet']
    link_names = ['Culvert']

    nodes = [Nodes(sim)[names] for names in node_names]
    links = [Links(sim)[names] for names in link_names]

    # type, area, length, orifice_coeff, free_weir_coeff, submerged_weir_coeff
    opening0 = nodes[0].create_opening(4, 1.0, 1.0, 0.6, 1.6, 1.0)
    opening1 = nodes[1].create_opening(4, 1.0, 1.0, 0.6, 1.6, 1.0)

    print("\n")
    print("n0 is coupled? ", nodes[0].is_coupled)
    print("n1 is coupled? ", nodes[1].is_coupled)

    nodes[0].overland_depth = 1.0
    nodes[0].coupling_area = 1.0

    # This step_advance should be an integer multiple of the routing step
    # which is set in the ,inp file. Currently set to 1s.
    # Should be able to interrogate sim to find out what the
    # routing stepsize is. Maybe should issue a warning if
    # step_advance is set lower than the routing step size.
    # Indeed maybe step_advance should just allow advance n routing steps?



    for i in range(5):
        ids = dict()
        nodes[0].overland_depth = 1.8854
        volumes = sim.coupling_step(1.0)
        print("Step:", i)
        print("Current time:", sim.current_time)
        print("overland depth: ",nodes[0].overland_depth)
        print(volumes)

def dynamic(t):
    return -0.02 * t + 1.8854

run_swmm()