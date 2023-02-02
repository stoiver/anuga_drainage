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

#------------------------------------------------------------------------------
# FILENAMES, MODEL DOMAIN and VARIABLES
#------------------------------------------------------------------------------

basename = 'terrain'
outname = 'domain_boyd_culvert'
meshname = 'terrain.msh'

dt = 1.0
ft = 400.0

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
# BOYD PIPE CULVERT
#------------------------------------------------------------------------------

losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
ep0 = numpy.array([296660.390,6180017.186]) 
ep1 = numpy.array([296649.976,6180038.872])
    
invert_elevations=[12.40,12.20]          
culvert = Boyd_pipe_operator(domain,
    losses=losses,
    diameter=1.0,
    end_points=[ep0, ep1],
    invert_elevations=invert_elevations,
    use_momentum_jet=False,
    use_velocity_head=False,
    manning=0.013,
    logging=True,
    label='boyd_pipe', 
    verbose=False)
    
#------------------------------------------------------------------------------
# APPLY FLOW
#------------------------------------------------------------------------------

line=[[296669.258,6179974.191],[296677.321,6179976.449]]
anuga.parallel.Inlet_operator(domain, line, 1.0)

#------------------------------------------------------------------------------
# SETUP BOUNDARY CONDITIONS
#------------------------------------------------------------------------------

print ('Available boundary tags', domain.get_boundary_tags())

Br = anuga.Reflective_boundary(domain)
Bd = anuga.Dirichlet_boundary([0,0,0])

domain.set_boundary({'west': Bd, 'south': Br, 'north': Bd, 'east': Bd})
    
#------------------------------------------------------------------------------
# EVOLVE SYSTEM THROUGH TIME
#------------------------------------------------------------------------------

import time
t0 = time.time()

for t in domain.evolve(yieldstep = dt, finaltime = ft):
    domain.print_timestepping_statistics()
    domain.print_boundary_statistics(quantities='stage')  

print ('Finished')
