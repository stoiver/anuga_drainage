#------------------------------------------------------------------------------
# IMPORT NECESSARY MODULES
#------------------------------------------------------------------------------
print (' ABOUT to Start Simulation:- Importing Modules')

import anuga, numpy, time, os, glob
from anuga.operators.rate_operators import Polygonal_rate_operator
from anuga import file_function, Polygon_function, read_polygon, create_domain_from_regions, Domain, Inlet_operator
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
outname = 'domain_no_pipes'
meshname = 'model/terrain.tsh'

#------------------------------------------------------------------------------
# CREATING COMPUTATIONAL DOMAIN
#------------------------------------------------------------------------------
CatchmentDictionary = {'model/kerb/kerb1.csv':0.01, 'model/kerb/kerb2.csv':0.01}
    
bounding_polygon = anuga.read_polygon('model/domain.csv')
interior_regions = anuga.read_polygon_dir(CatchmentDictionary, 'model/kerb')

domain = create_domain_from_regions(bounding_polygon,
    boundary_tags={'south': [0], 'east': [1], 'north': [2], 'west': [3]},
    maximum_triangle_area=0.1,
    interior_regions=interior_regions,
    mesh_filename=meshname,
    use_cache=False,
    verbose=False)

domain.set_minimum_storable_height(0.015)
domain.set_name(outname) 

print (domain.statistics())

#------------------------------------------------------------------------------
# APPLY MANNING'S ROUGHNESSES
#------------------------------------------------------------------------------

domain.set_quantity('friction', 0.03)

# Set a Initial Water Level over the Domain
domain.set_quantity('stage', 0)

domain.set_quantity('elevation', filename=basename+'.csv', alpha=0.99)

#------------------------------------------------------------------------------
# SETUP BOUNDARY CONDITIONS
#------------------------------------------------------------------------------

print ('Available boundary tags', domain.get_boundary_tags())

Br = anuga.Reflective_boundary(domain)  
Bd = anuga.Dirichlet_boundary([0,0,0])

domain.set_boundary({'exterior': Bd, 'west': Bd, 'south': Bd, 'north': Bd, 'east': Bd})
 
# ------------------------------------------------------------------------------
# Setup inject water
# ------------------------------------------------------------------------------

input1_anuga_region = Region(domain, radius=1.0, center=(305694.91,6188013.94))
input1_anuga_inlet_op = Inlet_operator(domain, input1_anuga_region, Q=0.102) # i made flow exactly the same as in DRAINS example



dt = 0.25    # yield step
ft = 400    # final timestep



for t in domain.evolve(yieldstep=dt, finaltime=ft):
    print('\n')
    domain.print_timestepping_statistics()

