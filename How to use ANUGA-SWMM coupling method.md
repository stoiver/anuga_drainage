# How to use ANUGA-SWMM coupling method

This is an instruction about how to use the **ANUGA-SWMM** coupling method. Since ANUGA only runs on Python, we use PySWMM in the coupling method. In this tutorial, we assume you have already knew how to use **ANUGA and SWMM**. If not, please find the relevant tutorials below.

This tutorial will run on the `test-Yshape_flooding_offset.py`, alone with the swmm file `2inlets_short.inp` script.

### Relevant materials:

- ANUGA Clinic: this is the simple tutorial how to use the ANUGA, [link](https://github.com/anuga-community/anuga-clinic)
- SWMM manual: [link](https://www.epa.gov/sites/production/files/2019-02/documents/epaswmm5_1_manual_master_8-2-15.pdf)
- PySWMM manual: [link](https://pyswmm.readthedocs.io/en/stable/reference/nodes.html)

## Mechanism of the coupling method:

The main idea of the coupling method is that both of models runs simultaneously and ANUGA takes charge of **overground part** and SWMM for the **underground part**. Also, models would transmit data for each time step.

## Y-shape model and environment set up

![](https://i.loli.net/2021/05/18/C6bLfkXmn1yE9Qa.jpg)

This is an simple y-shape topology image. The main idea of this model is to let water on the LHS flows to RHS by the gravity.


1. Environment settings:
   - left/right platform elevation: -3/-5m
   - underground pip network topology:
     - Node1 <-> Node2 <-> Node3
   - water inject: inject water at the blue arrow shows in the image.
2. Water loss definition:

   ![](https://latex2image-output.s3.amazonaws.com/img-Hs83u5vM.png)
   
   where `v` represents the velocity of inject water, `V` represent the volume of water. And the code shows below. 

``` python
total_volume_correct = t * input_velocity

# In this case, the nodes in SWMM cannot store water while this could be modified, so it only calculate the volume in the culverts.
total_volume_real = domain.get_water_volume()+culvert1.volume+culvert2.volume                        
loss = total_volume_real - total_volume_correct
```

## SWMM set up:
> About how to set up SWMM see the link below:
[SWMM inp file generation tips and meaning of some important PySWMM parameters](https://github.com/20-S2-2-C-Flood-Modelling/anuga_core/wiki/SWMM-inp-file-generation-tips-and-meaning-of-some-important-PySWMM-parameters)

loading the swmm setting by the following code.

```
from pyswmm import Simulation, Nodes, Links
sim = Simulation('./2inlets_short.inp')
sim.start()
```

Assignment of variables, down below shows two ways to variable assignment:

```
node_names = ['J1', 'J2', 'Out1']
link_names = ['C1', 'C2']
nodes = [Nodes(sim)[names] for names in node_names]
links = [Links(sim)[names] for names in link_names]

culvert1 = Links(sim)['C1']
culvert2 = Links(sim)['C2']

inlet1 = Nodes(sim)['J1']
inlet2 = Nodes(sim)['J2']
outlet = Nodes(sim)['Out1']
```



Initialise opening for nodes: the parameter of create_opening shows below

```
# create_opening(type, area, length, orifice_coeff, free_weir_coeff, submerged_weir_coeff)
for i in range(len(node_names)):
    opening = nodes[i].create_opening(4, math.pi, 1.0, 0.6, 1.6, 1.0)
    nodes[i].coupling_area = math.pi
    print('node opening? ', node_names[i], ' ', opening)
    print('node coupled? ', node_names[i], ' ', nodes[i].is_coupled)
```

## Coupling method:

In the coupling method, we only transmit ANUGA overland water depth of each node to SWMM, and then SWMM would transmit the changes of water volume (as velocity) back to ANUGA in each time step.

<img src="https://latex2image-output.s3.amazonaws.com/img-E4NpzrKE.png" div align=center />


```
for t in domain.evolve(yieldstep=1.0, finaltime=300.0):
    # overland_depth to SWMM
    inlet1.overland_depth = op_inlet1.inlet.get_average_depth()
    inlet2.overland_depth = op_inlet2.inlet.get_average_depth()
    volumes = sim.coupling_step(1.0)
    
    Q_inlet1 = -volumes[-1][1]['J1']
    Q_inlet2 = -volumes[-1][1]['J2']
    Q_outlet = outlet.total_inflow
    
    # set velocity at each node
    op_inlet1.set_Q(Q_inlet1+inlet1.flooding*0.0283168466)
    op_inlet2.set_Q(Q_inlet2+inlet2.flooding*0.0283168466)
    op_outlet.set_Q(Q_outlet)
```

`volumes`is a matrix that store all volume changes in each nodes for all time step that have been calculated. Note here, since `outlet` is a special node -- outflow node, that it cannot retrieve the velocity from the `volumes` matrix instead we need to retrieve from the `outlet.total_inflow`

*Flooding*: flooding is a special mechanism in SWMM that it will automatically drop water if the volume of water in swmm exceed its maximum capacity, but measurement unit here is cubic feet, so it need to transform to cubic meter.

