import math
from pyswmm import Nodes
from hymo import SWMMInpFile
from anuga import Inlet_operator,Region
import numpy as np

# Would a class be better here?

def initialize_inlets(domain,sim,inp,n_sides = 6,manhole_areas = [1],Q_in_0 = [1],rotation = 0):
    if n_sides < 3:
        raise NameError('A polygon should have at least 3 sides')


    inlet_operators = dict()
    elevation_list  = []
    circumferences  = []
    vertices        = []
    in_nodes = [node for node in Nodes(sim) if node.is_junction()]
    inlet_idx = 0
    one_segment = math.pi * 2 / n_sides
    for node in in_nodes:
        side_length = math.sqrt(4.0*manhole_areas[inlet_idx]*math.tan(math.pi/n_sides)/n_sides)
        radius = side_length/(2.0*math.sin(math.pi/n_sides))

# for node in Nodes(sim):

        # if node.is_junction():
        inlet_coordinates = [inp.coordinates.loc[node.nodeid].X_Coord, inp.coordinates.loc[node.nodeid].Y_Coord]

        vertex = [
            (math.sin(one_segment * i + rotation) * radius,
            math.cos(one_segment * i + rotation) * radius)
            for i in range(n_sides)]

        vertex = [[sum(pair) for pair in zip(point, inlet_coordinates)]
                for point in vertex]
        
        
        # inlet_operators[node.nodeid] = Inlet_operator(domain, Region(domain, radius=np.sqrt(1.167/np.pi), center=(inp.coordinates.loc[node.nodeid].X_Coord, inp.coordinates.loc[node.nodeid].Y_Coord),expand_polygon=True), Q_in_0[inlet_idx], zero_velocity=True)
        inlet_operators[node.nodeid] = Inlet_operator(domain, Region(domain,polygon = vertex,expand_polygon = True), Q_in_0[inlet_idx], zero_velocity=True)


        elevation_list.append(inlet_operators[node.nodeid].inlet.get_average_elevation())

        inlet_idx += 1
        circumferences.append(n_sides*side_length)
        vertices.append(vertex)

    vertices       = np.array(vertices)
    anuga_elevs    = np.array(elevation_list)
    circumferences = np.array(circumferences)

    return inlet_operators,anuga_elevs,circumferences,vertices

def n_sided_inlet(n_sides, manhole_areas = [1], rotation=0, inlet_coordinate=None):
    vertices = side_lengths = circumferences = []
    for area in manhole_areas:
        one_segment = math.pi * 2 / n_sides
        side_length = math.sqrt(4.0*area*math.tan(math.pi/n_sides)/n_sides)
        print(side_length)
        radius = side_length/(2.0*math.sin(math.pi/n_sides))
        vertex = [
            (math.sin(one_segment * i + rotation) * radius,
            math.cos(one_segment * i + rotation) * radius)
            for i in range(n_sides)]

        if inlet_coordinate:
            vertex = [[sum(pair) for pair in zip(point, inlet_coordinate)]
                    for point in vertices]
        circumference = n_sides*side_length
        
        vertices.append(vertex)
        circumferences.append(circumference)
        side_lengths.append(side_length)
        
    return vertices,circumferences,side_length