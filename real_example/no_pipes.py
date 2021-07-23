""" 
Strathhearn Ave flow
"""
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
outname = 'no_pipes'
meshname = 'terrain.tsh'

#------------------------------------------------------------------------------
# CREATING MESH
#------------------------------------------------------------------------------

riverWall_csv_files = glob.glob('kerb/*.csv') # Make a list of the csv files in BREAKLINES
(riverWalls, riverWall_parameters) = su.readListOfRiverWalls(riverWall_csv_files)

bounding_polygon = anuga.read_polygon('domain.csv')

create_mesh_from_regions(bounding_polygon, 
    boundary_tags={'south': [0], 'east': [1], 'north': [2], 'west': [3]},
    maximum_triangle_area=0.05,
    breaklines=riverWalls.values(),
    filename=meshname,
    use_cache=False,
    verbose=True)

#------------------------------------------------------------------------------
# SETUP COMPUTATIONAL DOMAIN
#------------------------------------------------------------------------------

domain = anuga.Domain(meshname, use_cache=True, verbose=True)
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

domain.set_boundary({'interior': Br, 'exterior': Bd, 'west': Bd, 'south': Bd, 'north': Bd, 'east': Bd})
 
# ------------------------------------------------------------------------------
# Setup inject water
# ------------------------------------------------------------------------------

input1_anuga_region = Region(domain, radius=1.0, center=(305692.98,6188014.56))
input1_anuga_inlet_op = Inlet_operator(domain, input1_anuga_region, Q=0.1) 



dt = 1    # yield step
ft = 100  # final timestep



for t in domain.evolve(yieldstep=dt, finaltime=ft):
    print('\n')
    domain.print_timestepping_statistics()

