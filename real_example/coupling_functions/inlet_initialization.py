import math
from pyswmm import Nodes
from hymo import SWMMInpFile
from anuga import Inlet_operator,Region
import numpy as np

def n_sided_inlet(n_sides, area, inlet_coordinate, rotation):
    # Computes the vertex coordinates and side length of a regular polygon with:
    # Number of sides = n_sides
    # Area = area
    if n_sides < 3:
        raise RuntimeError('A polygon should have at least 3 sides')

    one_segment = math.pi * 2 / n_sides
    side_length = math.sqrt(4.0*area*math.tan(math.pi/n_sides)/n_sides)
    
    radius = side_length/(2.0*math.sin(math.pi/n_sides))

    vertex = [
        (math.sin(one_segment * i + rotation) * radius,
        math.cos(one_segment * i + rotation) * radius)
        for i in range(n_sides)]

    vertex = [[sum(pair) for pair in zip(point, inlet_coordinate)]
            for point in vertex]
        
    return vertex, side_length

def initialize_inlets(domain, sim, inp, n_sides = 6, manhole_areas = [1], Q_in_0 = [1], rotation = 0):
    if n_sides < 3:
        raise RuntimeError('A polygon should have at least 3 sides')

    if not(isinstance(manhole_areas, int) or isinstance(manhole_areas,float) or isinstance(manhole_areas,list) or isinstance(manhole_areas,np.ndarray)):
        raise RuntimeError('Invalid ')

    inlet_operators = dict()
    elevation_list  = []
    circumferences  = []
    polygons        = []
    in_nodes        = [node for node in Nodes(sim) if node.is_junction()]

    for inlet_idx, node in enumerate(in_nodes):

        if isinstance(manhole_areas,list) or isinstance(manhole_areas,np.ndarray):
            inlet_area = manhole_areas[inlet_idx] 

        inlet_coordinates = [inp.coordinates.loc[node.nodeid].X_Coord, inp.coordinates.loc[node.nodeid].Y_Coord]
        vertices, side_length = n_sided_inlet(n_sides, inlet_area, inlet_coordinates, rotation)
        
        inlet_operators[node.nodeid] = Inlet_operator(domain, Region(domain,polygon = vertices,expand_polygon = True), Q_in_0[inlet_idx], zero_velocity=True)

        elevation_list.append(inlet_operators[node.nodeid].inlet.get_average_elevation())
        circumferences.append(n_sides*side_length)
        polygons.append(vertices)

    polygons         = np.array(polygons)
    inlet_elevations = np.array(elevation_list)
    circumferences   = np.array(circumferences)

    return inlet_operators,inlet_elevations,circumferences,vertices
